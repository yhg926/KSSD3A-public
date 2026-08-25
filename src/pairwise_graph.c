#include "pairwise_graph.h"

#include <err.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "command_ani.h"
#include "kssdlib_sort.h"

static double pairwise_clean_distance(double distance)
{
	if (!isfinite(distance))
		return 1.0;
	if (distance < 0.0)
		return 0.0;
	return distance;
}

static double pairwise_mash_distance(int kmer_len, size_t x_count,
									 size_t y_count, size_t xny_count)
{
	if (kmer_len <= 0 || xny_count == 0)
		return 1.0;
	const double denom = (double)x_count + (double)y_count - (double)xny_count;
	if (denom <= 0.0)
		return 1.0;
	const double jaccard = (double)xny_count / denom;
	const double transformed = (2.0 * jaccard) / (1.0 + jaccard);
	if (transformed <= 0.0)
		return 1.0;
	return pairwise_clean_distance(-log(transformed) / (double)kmer_len);
}

static double pairwise_aaf_distance(int kmer_len, size_t x_count,
									size_t y_count, size_t xny_count)
{
	if (kmer_len <= 0 || xny_count == 0)
		return 1.0;
	const size_t min_count = x_count < y_count ? x_count : y_count;
	if (min_count == 0)
		return 1.0;
	const double containment = (double)xny_count / (double)min_count;
	if (containment <= 0.0)
		return 1.0;
	return pairwise_clean_distance(-log(containment) / (double)kmer_len);
}

const char *pairwise_metric_name(pairwise_metric_t metric)
{
	switch (metric) {
	case PAIRWISE_METRIC_CTX_MOE:
		return "ctx-moe";
	case PAIRWISE_METRIC_CTX_NAIVE:
		return "ctx-naive";
	case PAIRWISE_METRIC_P_DIST:
		return "p_dist";
	case PAIRWISE_METRIC_MASH:
		return "mash";
	case PAIRWISE_METRIC_AAF:
		return "aaf";
	}
	return "ctx-moe";
}

bool pairwise_metric_from_string(const char *arg, pairwise_metric_t *metric_out)
{
	if (!arg || !metric_out)
		return false;
	if (strcmp(arg, "ctx-moe") == 0 || strcmp(arg, "moe") == 0) {
		*metric_out = PAIRWISE_METRIC_CTX_MOE;
		return true;
	}
	if (strcmp(arg, "ctx-naive") == 0 || strcmp(arg, "naive") == 0) {
		*metric_out = PAIRWISE_METRIC_CTX_NAIVE;
		return true;
	}
	if (strcmp(arg, "p_dist") == 0 || strcmp(arg, "pdist") == 0 ||
		strcmp(arg, "p-distance") == 0 || strcmp(arg, "micro-obj") == 0 ||
		strcmp(arg, "micro_obj_rate") == 0) {
		*metric_out = PAIRWISE_METRIC_P_DIST;
		return true;
	}
	if (strcmp(arg, "mash") == 0 || strcmp(arg, "mashd") == 0 ||
		strcmp(arg, "0") == 0) {
		*metric_out = PAIRWISE_METRIC_MASH;
		return true;
	}
	if (strcmp(arg, "aaf") == 0 || strcmp(arg, "aafd") == 0 ||
		strcmp(arg, "1") == 0) {
		*metric_out = PAIRWISE_METRIC_AAF;
		return true;
	}
	return false;
}

bool pairwise_metric_is_context(pairwise_metric_t metric)
{
	return metric == PAIRWISE_METRIC_CTX_MOE ||
		   metric == PAIRWISE_METRIC_CTX_NAIVE ||
		   metric == PAIRWISE_METRIC_P_DIST;
}

const char *pairwise_dedup_strategy_name(pairwise_dedup_strategy_t strategy)
{
	switch (strategy) {
	case PAIRWISE_DEDUP_GREEDY:
		return "greedy";
	case PAIRWISE_DEDUP_COMPLETE_LINKAGE:
		return "full-linkage";
	}
	return "greedy";
}

bool pairwise_dedup_strategy_from_string(const char *arg,
										 pairwise_dedup_strategy_t *strategy_out)
{
	if (!arg || !strategy_out)
		return false;
	if (strcmp(arg, "greedy") == 0 || strcmp(arg, "rep") == 0 ||
		strcmp(arg, "representative") == 0 ||
		strcmp(arg, "rep-centric") == 0 ||
		strcmp(arg, "representative-centric") == 0) {
		*strategy_out = PAIRWISE_DEDUP_GREEDY;
		return true;
	}
	if (strcmp(arg, "full") == 0 || strcmp(arg, "full-linkage") == 0 ||
		strcmp(arg, "full_linkage") == 0 ||
		strcmp(arg, "complete") == 0 ||
		strcmp(arg, "complete-linkage") == 0 ||
		strcmp(arg, "complete_linkage") == 0 ||
		strcmp(arg, "clique") == 0) {
		*strategy_out = PAIRWISE_DEDUP_COMPLETE_LINKAGE;
		return true;
	}
	return false;
}

pairwise_metric_expr_t pairwise_metric_expr_single(pairwise_metric_t metric)
{
	pairwise_metric_expr_t expr = {
		.metrics = { metric },
		.count = 1,
		.op = PAIRWISE_METRIC_EXPR_SINGLE,
		.label = { 0 },
	};
	snprintf(expr.label, sizeof(expr.label), "%s", pairwise_metric_name(metric));
	return expr;
}

static bool pairwise_expr_append_label(pairwise_metric_expr_t *expr)
{
	const char *sep = expr->op == PAIRWISE_METRIC_EXPR_AND ? "&" :
					  expr->op == PAIRWISE_METRIC_EXPR_OR ? "|" : "";
	expr->label[0] = '\0';
	for (size_t i = 0; i < expr->count; ++i) {
		const size_t used = strlen(expr->label);
		const size_t remaining = sizeof(expr->label) - used;
		const int written =
			snprintf(expr->label + used, remaining,
					 "%s%s", i == 0 ? "" : sep,
					 pairwise_metric_name(expr->metrics[i]));
		if (written < 0 || (size_t)written >= remaining)
			return false;
	}
	return true;
}

bool pairwise_metric_expr_from_string(const char *arg,
									  pairwise_metric_expr_t *expr_out)
{
	if (!arg || !expr_out || arg[0] == '\0')
		return false;
	char compact[PAIRWISE_METRIC_EXPR_LABEL_LEN];
	size_t n = 0;
	for (const char *p = arg; *p != '\0'; ++p) {
		if (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r')
			continue;
		if (n + 1 >= sizeof(compact))
			return false;
		compact[n++] = *p;
	}
	compact[n] = '\0';
	if (n == 0)
		return false;

	bool has_and = false;
	bool has_or = false;
	for (size_t i = 0; i < n; ++i) {
		has_and = has_and || compact[i] == '&';
		has_or = has_or || compact[i] == '|';
	}
	if (has_and && has_or)
		return false;
	if (has_and || has_or) {
		const char op = has_and ? '&' : '|';
		if (compact[0] == op || compact[n - 1] == op)
			return false;
		for (size_t i = 1; i < n; ++i) {
			if (compact[i] == op && compact[i - 1] == op)
				return false;
		}
	}

	pairwise_metric_expr_t expr = {0};
	expr.op = has_and ? PAIRWISE_METRIC_EXPR_AND :
			  has_or ? PAIRWISE_METRIC_EXPR_OR :
			  PAIRWISE_METRIC_EXPR_SINGLE;
	const char delim[2] = { has_and ? '&' : has_or ? '|' : '\0', '\0' };
	char *saveptr = NULL;
	char *token = expr.op == PAIRWISE_METRIC_EXPR_SINGLE
					  ? compact
					  : strtok_r(compact, delim, &saveptr);
	while (token) {
		if (expr.count >= PAIRWISE_METRIC_EXPR_MAX)
			return false;
		if (!pairwise_metric_from_string(token, &expr.metrics[expr.count]))
			return false;
		++expr.count;
		token = expr.op == PAIRWISE_METRIC_EXPR_SINGLE
					? NULL
					: strtok_r(NULL, delim, &saveptr);
	}
	if (expr.count == 0)
		return false;
	if (expr.count == 1)
		expr.op = PAIRWISE_METRIC_EXPR_SINGLE;
	if (!pairwise_expr_append_label(&expr))
		return false;
	*expr_out = expr;
	return true;
}

bool pairwise_metric_expr_is_combined(const pairwise_metric_expr_t *expr)
{
	return expr && expr->count > 1 && expr->op != PAIRWISE_METRIC_EXPR_SINGLE;
}

bool pairwise_metric_expr_uses_context(const pairwise_metric_expr_t *expr)
{
	if (!expr)
		return false;
	for (size_t i = 0; i < expr->count; ++i) {
		if (pairwise_metric_is_context(expr->metrics[i]))
			return true;
	}
	return false;
}

pairwise_metric_t pairwise_metric_expr_primary(const pairwise_metric_expr_t *expr)
{
	return expr && expr->count > 0 ? expr->metrics[0] : PAIRWISE_METRIC_CTX_MOE;
}

const char *pairwise_metric_expr_name(const pairwise_metric_expr_t *expr)
{
	if (!expr || expr->label[0] == '\0')
		return pairwise_metric_name(pairwise_metric_expr_primary(expr));
	return expr->label;
}

void pairwise_check_compatible(const unify_sketch_t *ref, const unify_sketch_t *qry)
{
	if (!ref || !qry)
		errx(EINVAL, "%s(): NULL sketch", __func__);
	if (ref->stat_type != qry->stat_type)
		errx(EXIT_FAILURE, "%s(): ref sketch type %u != qry %u",
			 __func__, ref->stat_type, qry->stat_type);
	if (ref->hash_id != qry->hash_id)
		errx(EXIT_FAILURE, "%s(): ref hash_id %u != qry %u",
			 __func__, ref->hash_id, qry->hash_id);
}

void pairwise_prepare_lco_model(const unify_sketch_t *sketch)
{
	if (!sketch || sketch->stat_type != 2)
		return;
	dim_sketch_stat_t stat = sketch->stats.lco_stat_val;
	const_comask_init(&stat);
	ani_model_drfold = stat.drfold;
}

uint32_t pairwise_count_ctx_runs_sorted_ctxobj64(const uint64_t *a, size_t n)
{
	const uint8_t nobjbits = Bitslen.obj;
	uint32_t count = 0;
	size_t i = 0;
	while (i < n) {
		const uint64_t ctx = a[i] >> nobjbits;
		do {
			++i;
		} while (i < n && (a[i] >> nobjbits) == ctx);
		++count;
	}
	return count;
}

pairwise_eval_t pairwise_eval_arrays(pairwise_metric_t metric,
									 const dim_sketch_stat_t *stat,
									 const uint64_t *qry, size_t qry_n,
									 const uint64_t *ref, size_t ref_n,
									 uint32_t qry_ctx_count,
									 uint32_t ref_ctx_count,
									 bool ignore_ref_conflicts)
{
	pairwise_eval_t eval = {
		.distance = 1.0,
		.similarity = 0.0,
		.xny_ctx = 0,
		.n_diff_obj = 0,
		.n_diff_obj_section = 0,
		.n_mut2_ctx = 0,
		.af_qry = 0.0,
		.af_ref = 0.0,
		.max_af = 0.0,
		.valid = false,
	};
	if (!stat || !qry || !ref || qry_n == 0 || ref_n == 0)
		return eval;

	ani_features_t features;
	get_ani_features_ctx_min_over_conflicts_both_filtered(qry, qry_n, ref, ref_n,
														  ignore_ref_conflicts,
														  &features);
	eval.xny_ctx = features.XnY_ctx;
	eval.n_diff_obj = features.N_diff_obj;
	eval.n_diff_obj_section = features.N_diff_obj_section;
	eval.n_mut2_ctx = features.N_mut2_ctx;
	eval.af_qry = qry_ctx_count > 0
					  ? (double)features.XnY_ctx / (double)qry_ctx_count
					  : 0.0;
	eval.af_ref = ref_ctx_count > 0
					  ? (double)features.XnY_ctx / (double)ref_ctx_count
					  : 0.0;
	eval.max_af = eval.af_qry > eval.af_ref ? eval.af_qry : eval.af_ref;
	if (features.XnY_ctx == 0)
		return eval;

	if (metric == PAIRWISE_METRIC_MASH || metric == PAIRWISE_METRIC_AAF) {
		const int ctx_k = Bitslen.ctx > 0 ? Bitslen.ctx / 2 : stat->klen;
		if (metric == PAIRWISE_METRIC_MASH)
			eval.distance = pairwise_mash_distance(ctx_k, qry_ctx_count,
												   ref_ctx_count, features.XnY_ctx);
		else
			eval.distance = pairwise_aaf_distance(ctx_k, qry_ctx_count,
												  ref_ctx_count, features.XnY_ctx);
	} else if (metric == PAIRWISE_METRIC_CTX_NAIVE) {
		eval.distance = pairwise_clean_distance(get_naive_dist(&features));
	} else if (metric == PAIRWISE_METRIC_P_DIST) {
		eval.distance = pairwise_clean_distance(get_p_dist(&features));
	} else {
		eval.distance = pairwise_clean_distance(lm3ways_dist_from_features(&features));
	}
	eval.similarity = 1.0 - eval.distance;
	eval.valid = true;
	return eval;
}

pairwise_eval_t pairwise_eval_samples(pairwise_metric_t metric,
									  const unify_sketch_t *qry, uint32_t qn,
									  const unify_sketch_t *ref, uint32_t rn,
									  bool ignore_ref_conflicts)
{
	pairwise_check_compatible(ref, qry);
	if (qn >= (uint32_t)qry->infile_num || rn >= (uint32_t)ref->infile_num)
		errx(EINVAL, "%s(): sample index out of range", __func__);

	const uint64_t *qry_arr = qry->comb_sketch + qry->sketch_index[qn];
	const size_t qry_len = (size_t)(qry->sketch_index[qn + 1] -
								   qry->sketch_index[qn]);
	const uint64_t *ref_arr = ref->comb_sketch + ref->sketch_index[rn];
	const size_t ref_len = (size_t)(ref->sketch_index[rn + 1] -
								   ref->sketch_index[rn]);

	if (qry->stat_type == 1) {
		if (pairwise_metric_is_context(metric))
			errx(EINVAL, "%s(): context-object metrics require lco sketches",
				 __func__);
		pairwise_eval_t eval = {
			.distance = 1.0,
			.similarity = 0.0,
			.xny_ctx = 0,
			.valid = false,
		};
		const size_t overlap = count_overlaps(qry_arr, qry_len, ref_arr, ref_len);
		eval.xny_ctx = overlap > UINT32_MAX ? UINT32_MAX : (uint32_t)overlap;
		eval.af_qry = qry_len > 0 ? (double)overlap / (double)qry_len : 0.0;
		eval.af_ref = ref_len > 0 ? (double)overlap / (double)ref_len : 0.0;
		eval.max_af = eval.af_qry > eval.af_ref ? eval.af_qry : eval.af_ref;
		if (overlap == 0)
			return eval;
		if (metric == PAIRWISE_METRIC_MASH)
			eval.distance = pairwise_mash_distance(qry->kmerlen, qry_len, ref_len, overlap);
		else
			eval.distance = pairwise_aaf_distance(qry->kmerlen, qry_len, ref_len, overlap);
		eval.similarity = 1.0 - eval.distance;
		eval.valid = true;
		return eval;
	}

	pairwise_prepare_lco_model(ref);
	const uint32_t qry_ctx_count = qry->conflict
									  ? pairwise_count_ctx_runs_sorted_ctxobj64(qry_arr, qry_len)
									  : (uint32_t)qry_len;
	const uint32_t ref_ctx_count = ref->conflict
									  ? pairwise_count_ctx_runs_sorted_ctxobj64(ref_arr, ref_len)
									  : (uint32_t)ref_len;
	return pairwise_eval_arrays(metric, &ref->stats.lco_stat_val,
								qry_arr, qry_len, ref_arr, ref_len,
								qry_ctx_count, ref_ctx_count,
								ignore_ref_conflicts && ref->conflict);
}

static pairwise_eval_t pairwise_eval_expr_merge(const pairwise_metric_expr_t *expr,
												const pairwise_eval_t *evals,
												size_t n)
{
	pairwise_eval_t merged = {
		.distance = 1.0,
		.similarity = 0.0,
		.valid = false,
	};
	if (!expr || n == 0)
		return merged;
	if (expr->op == PAIRWISE_METRIC_EXPR_AND) {
		for (size_t i = 0; i < n; ++i) {
			if (!evals[i].valid)
				return merged;
			if (!merged.valid || evals[i].distance > merged.distance)
				merged = evals[i];
		}
	} else {
		for (size_t i = 0; i < n; ++i) {
			if (!evals[i].valid)
				continue;
			if (!merged.valid || evals[i].distance < merged.distance)
				merged = evals[i];
		}
	}
	if (merged.valid)
		merged.similarity = 1.0 - merged.distance;
	return merged;
}

pairwise_eval_t pairwise_eval_expr_arrays(const pairwise_metric_expr_t *expr,
										  const dim_sketch_stat_t *stat,
										  const uint64_t *qry, size_t qry_n,
										  const uint64_t *ref, size_t ref_n,
										  uint32_t qry_ctx_count,
										  uint32_t ref_ctx_count,
										  bool ignore_ref_conflicts)
{
	if (!expr || expr->count == 0)
		return pairwise_eval_arrays(PAIRWISE_METRIC_CTX_MOE, stat, qry, qry_n,
									ref, ref_n, qry_ctx_count, ref_ctx_count,
									ignore_ref_conflicts);
	if (!pairwise_metric_expr_is_combined(expr))
		return pairwise_eval_arrays(expr->metrics[0], stat, qry, qry_n,
									ref, ref_n, qry_ctx_count, ref_ctx_count,
									ignore_ref_conflicts);

	pairwise_eval_t evals[PAIRWISE_METRIC_EXPR_MAX];
	for (size_t i = 0; i < expr->count; ++i) {
		evals[i] = pairwise_eval_arrays(expr->metrics[i], stat, qry, qry_n,
										ref, ref_n, qry_ctx_count, ref_ctx_count,
										ignore_ref_conflicts);
	}
	return pairwise_eval_expr_merge(expr, evals, expr->count);
}

pairwise_eval_t pairwise_eval_expr_samples(const pairwise_metric_expr_t *expr,
										   const unify_sketch_t *qry, uint32_t qn,
										   const unify_sketch_t *ref, uint32_t rn,
										   bool ignore_ref_conflicts)
{
	if (!expr || expr->count == 0 || !pairwise_metric_expr_is_combined(expr))
		return pairwise_eval_samples(pairwise_metric_expr_primary(expr),
									 qry, qn, ref, rn, ignore_ref_conflicts);

	pairwise_check_compatible(ref, qry);
	if (qn >= (uint32_t)qry->infile_num || rn >= (uint32_t)ref->infile_num)
		errx(EINVAL, "%s(): sample index out of range", __func__);

	const uint64_t *qry_arr = qry->comb_sketch + qry->sketch_index[qn];
	const size_t qry_len = (size_t)(qry->sketch_index[qn + 1] -
								   qry->sketch_index[qn]);
	const uint64_t *ref_arr = ref->comb_sketch + ref->sketch_index[rn];
	const size_t ref_len = (size_t)(ref->sketch_index[rn + 1] -
								   ref->sketch_index[rn]);

	if (qry->stat_type == 1) {
		if (pairwise_metric_expr_uses_context(expr))
			errx(EINVAL, "%s(): context-object metrics require lco sketches",
				 __func__);
		pairwise_eval_t evals[PAIRWISE_METRIC_EXPR_MAX];
		for (size_t i = 0; i < expr->count; ++i)
			evals[i] = pairwise_eval_samples(expr->metrics[i], qry, qn, ref, rn,
											 ignore_ref_conflicts);
		return pairwise_eval_expr_merge(expr, evals, expr->count);
	}

	pairwise_prepare_lco_model(ref);
	const uint32_t qry_ctx_count = qry->conflict
									  ? pairwise_count_ctx_runs_sorted_ctxobj64(qry_arr, qry_len)
									  : (uint32_t)qry_len;
	const uint32_t ref_ctx_count = ref->conflict
									  ? pairwise_count_ctx_runs_sorted_ctxobj64(ref_arr, ref_len)
									  : (uint32_t)ref_len;
	return pairwise_eval_expr_arrays(expr, &ref->stats.lco_stat_val,
									 qry_arr, qry_len, ref_arr, ref_len,
									 qry_ctx_count, ref_ctx_count,
									 ignore_ref_conflicts && ref->conflict);
}

bool pairwise_eval_passes_edge(const pairwise_eval_t *eval, double cut,
							   uint32_t ctxcut, double max_afcut)
{
	if (!eval || !eval->valid)
		return false;
	if (eval->distance >= cut)
		return false;
	if (eval->xny_ctx < ctxcut)
		return false;
	if (eval->max_af < max_afcut)
		return false;
	return true;
}

static bool pairwise_meta_rankable(const unify_sketch_t *sketch, int idx)
{
	if (!sketch || !sketch->infile_meta || idx < 0 || idx >= sketch->infile_num)
		return false;
	const infile_meta_t *meta = &sketch->infile_meta[idx];
	return meta->meta_fmt_version == KSSD_INFILE_META_VERSION &&
		   meta->total_length_bp > 0;
}

static double pairwise_finite_or(double value, double fallback)
{
	return isfinite(value) ? value : fallback;
}

static int pairwise_positive_lower_u32_cmp(uint32_t a, uint32_t b)
{
	if (a == b)
		return 0;
	if (a == 0)
		return -1;
	if (b == 0)
		return 1;
	return a < b ? 1 : -1;
}

int pairwise_quality_compare(const unify_sketch_t *sketch, int a, int b)
{
	if (a == b)
		return 0;

	const bool a_rankable = pairwise_meta_rankable(sketch, a);
	const bool b_rankable = pairwise_meta_rankable(sketch, b);
	if (a_rankable != b_rankable)
		return a_rankable ? 1 : -1;

	if (a_rankable && b_rankable) {
		const infile_meta_t *ma = &sketch->infile_meta[a];
		const infile_meta_t *mb = &sketch->infile_meta[b];
		const bool a_fasta = ma->infile_fmt == KSSD_INFILE_FMT_FASTA;
		const bool b_fasta = mb->infile_fmt == KSSD_INFILE_FMT_FASTA;
		if (a_fasta != b_fasta)
			return a_fasta ? 1 : -1;

		const double a_asm = pairwise_finite_or((double)ma->asm_level, -1.0);
		const double b_asm = pairwise_finite_or((double)mb->asm_level, -1.0);
		if (a_asm != b_asm)
			return a_asm > b_asm ? 1 : -1;

		if (ma->total_length_bp != mb->total_length_bp)
			return ma->total_length_bp > mb->total_length_bp ? 1 : -1;

		const int record_cmp =
			pairwise_positive_lower_u32_cmp(ma->record_count, mb->record_count);
		if (record_cmp != 0)
			return record_cmp;

		if (ma->median_length_bp != mb->median_length_bp)
			return ma->median_length_bp > mb->median_length_bp ? 1 : -1;

		const double a_cv = pairwise_finite_or((double)ma->length_cv, INFINITY);
		const double b_cv = pairwise_finite_or((double)mb->length_cv, INFINITY);
		if (a_cv != b_cv)
			return a_cv < b_cv ? 1 : -1;
	}

	const uint64_t a_entries = sketch->sketch_index[a + 1] - sketch->sketch_index[a];
	const uint64_t b_entries = sketch->sketch_index[b + 1] - sketch->sketch_index[b];
	if (a_entries != b_entries)
		return a_entries > b_entries ? 1 : -1;

	return a < b ? 1 : -1;
}

void pairwise_component_result_free(pairwise_component_result_t *result)
{
	if (!result)
		return;
	free(result->parent);
	free(result->component_size);
	free(result->representative);
	free(result->remove_sample);
	free(result->rank);
	memset(result, 0, sizeof(*result));
}

int pairwise_component_find(const pairwise_component_result_t *result, int x)
{
	if (!result || !result->parent || x < 0 || x >= result->n)
		return -1;
	int root = x;
	while (result->parent[root] != root)
		root = result->parent[root];
	return root;
}

static int pairwise_component_find_mut(int *parent, int x)
{
	while (parent[x] != x) {
		parent[x] = parent[parent[x]];
		x = parent[x];
	}
	return x;
}

static void pairwise_component_union(int *parent, uint8_t *rank, int a, int b)
{
	int ra = pairwise_component_find_mut(parent, a);
	int rb = pairwise_component_find_mut(parent, b);
	if (ra == rb)
		return;
	if (rank[ra] < rank[rb])
		parent[ra] = rb;
	else if (rank[ra] > rank[rb])
		parent[rb] = ra;
	else {
		parent[rb] = ra;
		++rank[ra];
	}
}

void pairwise_component_result_init(pairwise_component_result_t *result, int n)
{
	if (!result || n < 0)
		errx(EINVAL, "%s(): invalid component arguments", __func__);
	memset(result, 0, sizeof(*result));
	result->n = n;
	if (n == 0)
		return;
	result->parent = malloc((size_t)n * sizeof(result->parent[0]));
	result->rank = calloc((size_t)n, sizeof(result->rank[0]));
	result->component_size = calloc((size_t)n, sizeof(result->component_size[0]));
	result->representative = malloc((size_t)n * sizeof(result->representative[0]));
	result->remove_sample = calloc((size_t)n, sizeof(result->remove_sample[0]));
	if (!result->parent || !result->rank || !result->component_size ||
		!result->representative || !result->remove_sample)
		err(EXIT_FAILURE, "%s(): OOM component bookkeeping", __func__);
	for (int i = 0; i < n; ++i) {
		result->parent[i] = i;
		result->representative[i] = -1;
	}
}

void pairwise_component_result_union(pairwise_component_result_t *result, int a, int b)
{
	if (!result || !result->parent || !result->rank ||
		a < 0 || b < 0 || a >= result->n || b >= result->n)
		errx(EINVAL, "%s(): invalid component union", __func__);
	pairwise_component_union(result->parent, result->rank, a, b);
}

void pairwise_component_result_finalize(pairwise_component_result_t *result,
										void *quality_ctx,
										pairwise_quality_compare_fn quality_compare)
{
	if (!result || !quality_compare)
		errx(EINVAL, "%s(): invalid component-finalize arguments", __func__);
	for (int i = 0; i < result->n; ++i) {
		const int root = pairwise_component_find_mut(result->parent, i);
		result->parent[i] = root;
		++result->component_size[root];
	}

	for (int i = 0; i < result->n; ++i) {
		const int root = pairwise_component_find_mut(result->parent, i);
		if (result->representative[root] < 0 ||
			quality_compare(quality_ctx, i, result->representative[root]) > 0)
			result->representative[root] = i;
	}

	for (int i = 0; i < result->n; ++i) {
		const int root = pairwise_component_find_mut(result->parent, i);
		if (i == root && result->component_size[root] > 1)
			++result->duplicate_clusters;
		if (i != result->representative[root]) {
			result->remove_sample[i] = true;
			++result->removed_samples;
		}
	}
}

void pairwise_edge_list_free(pairwise_edge_list_t *list)
{
	if (!list)
		return;
	free(list->edges);
	memset(list, 0, sizeof(*list));
}

void pairwise_edge_list_add(pairwise_edge_list_t *list, int a, int b)
{
	if (!list)
		errx(EINVAL, "%s(): invalid edge list", __func__);
	if (a == b)
		return;
	if (list->n == list->cap) {
		const size_t new_cap = list->cap > 0 ? list->cap * 2 : 1024;
		pairwise_edge_pair_t *new_edges =
			realloc(list->edges, new_cap * sizeof(new_edges[0]));
		if (!new_edges)
			err(EXIT_FAILURE, "%s(): OOM edge list", __func__);
		list->edges = new_edges;
		list->cap = new_cap;
	}
	list->edges[list->n++] = (pairwise_edge_pair_t){ .a = a, .b = b };
}

static bool pairwise_quality_better(void *quality_ctx,
									pairwise_quality_compare_fn quality_compare,
									int a, int b)
{
	const int cmp = quality_compare(quality_ctx, a, b);
	return cmp > 0 || (cmp == 0 && a < b);
}

static void pairwise_quality_merge_sort_rec(int *idx, int *tmp,
											int left, int right,
											void *quality_ctx,
											pairwise_quality_compare_fn quality_compare)
{
	if (right - left <= 1)
		return;
	const int mid = left + (right - left) / 2;
	pairwise_quality_merge_sort_rec(idx, tmp, left, mid,
									quality_ctx, quality_compare);
	pairwise_quality_merge_sort_rec(idx, tmp, mid, right,
									quality_ctx, quality_compare);
	int i = left;
	int j = mid;
	int k = left;
	while (i < mid && j < right) {
		if (pairwise_quality_better(quality_ctx, quality_compare, idx[i], idx[j]))
			tmp[k++] = idx[i++];
		else
			tmp[k++] = idx[j++];
	}
	while (i < mid)
		tmp[k++] = idx[i++];
	while (j < right)
		tmp[k++] = idx[j++];
	for (i = left; i < right; ++i)
		idx[i] = tmp[i];
}

static int *pairwise_quality_order(int n, void *quality_ctx,
								   pairwise_quality_compare_fn quality_compare)
{
	if (n <= 0)
		return NULL;
	int *idx = malloc((size_t)n * sizeof(idx[0]));
	int *tmp = malloc((size_t)n * sizeof(tmp[0]));
	if (!idx || !tmp)
		err(EXIT_FAILURE, "%s(): OOM quality order", __func__);
	for (int i = 0; i < n; ++i)
		idx[i] = i;
	pairwise_quality_merge_sort_rec(idx, tmp, 0, n,
									quality_ctx, quality_compare);
	free(tmp);
	return idx;
}

void pairwise_greedy_dedup_from_edges(int n,
									  const pairwise_edge_pair_t *edges,
									  size_t edge_n,
									  void *quality_ctx,
									  pairwise_quality_compare_fn quality_compare,
									  pairwise_component_result_t *result)
{
	if (!result || n < 0 || !quality_compare)
		errx(EINVAL, "%s(): invalid greedy-dedup arguments", __func__);
	pairwise_component_result_init(result, n);
	if (n == 0)
		return;

	size_t *degree = calloc((size_t)n, sizeof(degree[0]));
	size_t *offset = calloc((size_t)n + 1, sizeof(offset[0]));
	if (!degree || !offset)
		err(EXIT_FAILURE, "%s(): OOM greedy-dedup graph", __func__);
	for (size_t e = 0; e < edge_n; ++e) {
		const int a = edges[e].a;
		const int b = edges[e].b;
		if (a < 0 || b < 0 || a >= n || b >= n || a == b)
			continue;
		++degree[a];
		++degree[b];
	}
	for (int i = 0; i < n; ++i)
		offset[i + 1] = offset[i] + degree[i];
	int *adj = offset[n] > 0 ? malloc(offset[n] * sizeof(adj[0])) : NULL;
	size_t *cursor = malloc(((size_t)n + 1) * sizeof(cursor[0]));
	if ((offset[n] > 0 && !adj) || !cursor)
		err(EXIT_FAILURE, "%s(): OOM greedy-dedup adjacency", __func__);
	memcpy(cursor, offset, ((size_t)n + 1) * sizeof(cursor[0]));
	for (size_t e = 0; e < edge_n; ++e) {
		const int a = edges[e].a;
		const int b = edges[e].b;
		if (a < 0 || b < 0 || a >= n || b >= n || a == b)
			continue;
		adj[cursor[a]++] = b;
		adj[cursor[b]++] = a;
	}

	int *order = pairwise_quality_order(n, quality_ctx, quality_compare);
	for (int oi = 0; oi < n; ++oi) {
		const int rep = order[oi];
		if (result->remove_sample[rep])
			continue;
		result->parent[rep] = rep;
		result->representative[rep] = rep;
		result->component_size[rep] = 1;
		for (size_t ai = offset[rep]; ai < offset[rep + 1]; ++ai) {
			const int sample = adj[ai];
			if (sample == rep || result->remove_sample[sample])
				continue;
			result->remove_sample[sample] = true;
			result->parent[sample] = rep;
			result->representative[sample] = rep;
			++result->component_size[rep];
			++result->removed_samples;
		}
		if (result->component_size[rep] > 1)
			++result->duplicate_clusters;
	}
	for (int i = 0; i < n; ++i) {
		if (result->representative[i] < 0) {
			result->representative[i] = i;
			result->component_size[i] = 1;
		}
	}

	free(order);
	free(cursor);
	free(adj);
	free(offset);
	free(degree);
}

static size_t pairwise_square_bit_words(int n)
{
	if (n <= 0)
		return 0;
	const size_t nn = (size_t)n;
	if (nn > SIZE_MAX / nn)
		errx(EINVAL, "%s(): too many samples for full-linkage bitset", __func__);
	const size_t bits = nn * nn;
	if (bits > SIZE_MAX - 63)
		errx(EINVAL, "%s(): too many samples for full-linkage bitset", __func__);
	return (bits + 63) / 64;
}

static void pairwise_edge_bit_set(uint64_t *bits, int n, int a, int b)
{
	const size_t pos = (size_t)a * (size_t)n + (size_t)b;
	bits[pos >> 6] |= UINT64_C(1) << (pos & 63);
}

static bool pairwise_edge_bit_test(const uint64_t *bits, int n, int a, int b)
{
	if (a == b)
		return true;
	const size_t pos = (size_t)a * (size_t)n + (size_t)b;
	return (bits[pos >> 6] & (UINT64_C(1) << (pos & 63))) != 0;
}

void pairwise_complete_linkage_dedup_from_edges(int n,
												const pairwise_edge_pair_t *edges,
												size_t edge_n,
												void *quality_ctx,
												pairwise_quality_compare_fn quality_compare,
												pairwise_component_result_t *result)
{
	if (!result || n < 0 || !quality_compare)
		errx(EINVAL, "%s(): invalid full-linkage dedup arguments", __func__);
	pairwise_component_result_init(result, n);
	if (n == 0)
		return;

	const size_t bit_words = pairwise_square_bit_words(n);
	uint64_t *adj_bits = bit_words > 0 ? calloc(bit_words, sizeof(adj_bits[0])) : NULL;
	if (bit_words > 0 && !adj_bits)
		err(EXIT_FAILURE, "%s(): OOM full-linkage adjacency bitset", __func__);
	for (size_t e = 0; e < edge_n; ++e) {
		const int a = edges[e].a;
		const int b = edges[e].b;
		if (a < 0 || b < 0 || a >= n || b >= n || a == b)
			continue;
		pairwise_edge_bit_set(adj_bits, n, a, b);
		pairwise_edge_bit_set(adj_bits, n, b, a);
	}

	int *order = pairwise_quality_order(n, quality_ctx, quality_compare);
	int *members = malloc((size_t)n * sizeof(members[0]));
	bool *assigned = calloc((size_t)n, sizeof(assigned[0]));
	if (!members || !assigned)
		err(EXIT_FAILURE, "%s(): OOM full-linkage members", __func__);

	for (int oi = 0; oi < n; ++oi) {
		const int rep = order[oi];
		if (assigned[rep])
			continue;
		assigned[rep] = true;
		result->parent[rep] = rep;
		result->representative[rep] = rep;
		result->component_size[rep] = 1;
		members[0] = rep;
		int member_n = 1;

		for (int ci = 0; ci < n; ++ci) {
			const int sample = order[ci];
			if (sample == rep || assigned[sample])
				continue;
			bool linked_to_all = true;
			for (int mi = 0; mi < member_n; ++mi) {
				if (!pairwise_edge_bit_test(adj_bits, n, sample, members[mi])) {
					linked_to_all = false;
					break;
				}
			}
			if (!linked_to_all)
				continue;
			result->remove_sample[sample] = true;
			assigned[sample] = true;
			result->parent[sample] = rep;
			result->representative[sample] = rep;
			++result->component_size[rep];
			++result->removed_samples;
			members[member_n++] = sample;
		}
		if (result->component_size[rep] > 1)
			++result->duplicate_clusters;
	}
	for (int i = 0; i < n; ++i) {
		if (result->representative[i] < 0) {
			result->representative[i] = i;
			result->component_size[i] = 1;
		}
	}

	free(members);
	free(assigned);
	free(order);
	free(adj_bits);
}

void pairwise_build_components(int n, double cut, uint32_t ctxcut,
							   double max_afcut, void *cb_ctx,
							   pairwise_eval_pair_fn eval_pair,
							   pairwise_quality_compare_fn quality_compare,
							   pairwise_edge_observer_fn edge_observer,
							   void *edge_ctx,
							   pairwise_progress_fn progress_fn,
							   void *progress_ctx,
							   pairwise_component_result_t *result)
{
	if (!result || n < 0 || !eval_pair || !quality_compare)
		errx(EINVAL, "%s(): invalid component-build arguments", __func__);
	pairwise_component_result_init(result, n);
	if (n == 0)
		return;

	uint64_t completed_pairs = 0;
	const uint64_t total_pairs = (uint64_t)n * (uint64_t)(n - 1) / 2;
	for (int i = 0; i < n; ++i) {
		for (int j = i + 1; j < n; ++j) {
			const pairwise_eval_t eval = eval_pair(cb_ctx, i, j);
			if (!eval.valid || eval.distance >= cut)
				continue;
			++result->distance_edges;
			if (eval.xny_ctx < ctxcut) {
				++result->ctx_rejects;
				continue;
			}
			if (eval.max_af < max_afcut) {
				++result->max_af_rejects;
				continue;
			}
			pairwise_component_result_union(result, i, j);
			++result->duplicate_edges;
			if (edge_observer)
				edge_observer(edge_ctx, i, j, &eval);
		}
		completed_pairs += (uint64_t)(n - i - 1);
		if (progress_fn)
			progress_fn(progress_ctx, completed_pairs, total_pairs);
	}

	pairwise_component_result_finalize(result, cb_ctx, quality_compare);
}

static void pairwise_collect_dedup_edges(int n, double cut, uint32_t ctxcut,
										 double max_afcut, void *cb_ctx,
										 pairwise_eval_pair_fn eval_pair,
										 pairwise_edge_observer_fn edge_observer,
										 void *edge_ctx,
										 pairwise_progress_fn progress_fn,
										 void *progress_ctx,
										 pairwise_edge_list_t *edges,
										 uint64_t *distance_edges,
										 uint64_t *duplicate_edges,
										 uint64_t *ctx_rejects,
										 uint64_t *max_af_rejects)
{
	if (n < 0 || !eval_pair || !edges || !distance_edges ||
		!duplicate_edges || !ctx_rejects || !max_af_rejects)
		errx(EINVAL, "%s(): invalid dedup edge collection arguments", __func__);

	uint64_t completed_pairs = 0;
	const uint64_t total_pairs = (uint64_t)n * (uint64_t)(n - 1) / 2;
	*distance_edges = 0;
	*duplicate_edges = 0;
	*ctx_rejects = 0;
	*max_af_rejects = 0;
	for (int i = 0; i < n; ++i) {
		for (int j = i + 1; j < n; ++j) {
			const pairwise_eval_t eval = eval_pair(cb_ctx, i, j);
			if (!eval.valid || eval.distance >= cut)
				continue;
			++(*distance_edges);
			if (eval.xny_ctx < ctxcut) {
				++(*ctx_rejects);
				continue;
			}
			if (eval.max_af < max_afcut) {
				++(*max_af_rejects);
				continue;
			}
			pairwise_edge_list_add(edges, i, j);
			++(*duplicate_edges);
			if (edge_observer)
				edge_observer(edge_ctx, i, j, &eval);
		}
		completed_pairs += (uint64_t)(n - i - 1);
		if (progress_fn)
			progress_fn(progress_ctx, completed_pairs, total_pairs);
	}
}

void pairwise_build_greedy_dedup(int n, double cut, uint32_t ctxcut,
								 double max_afcut, void *cb_ctx,
								 pairwise_eval_pair_fn eval_pair,
								 pairwise_quality_compare_fn quality_compare,
								 pairwise_edge_observer_fn edge_observer,
								 void *edge_ctx,
								 pairwise_progress_fn progress_fn,
								 void *progress_ctx,
								 pairwise_component_result_t *result)
{
	if (!result || n < 0 || !eval_pair || !quality_compare)
		errx(EINVAL, "%s(): invalid greedy-dedup arguments", __func__);

	pairwise_edge_list_t edges = {0};
	uint64_t distance_edges = 0;
	uint64_t duplicate_edges = 0;
	uint64_t ctx_rejects = 0;
	uint64_t max_af_rejects = 0;
	pairwise_collect_dedup_edges(n, cut, ctxcut, max_afcut, cb_ctx, eval_pair,
								 edge_observer, edge_ctx,
								 progress_fn, progress_ctx, &edges,
								 &distance_edges, &duplicate_edges,
								 &ctx_rejects, &max_af_rejects);

	pairwise_greedy_dedup_from_edges(n, edges.edges, edges.n,
									 cb_ctx, quality_compare, result);
	result->distance_edges = distance_edges;
	result->duplicate_edges = duplicate_edges;
	result->ctx_rejects = ctx_rejects;
	result->max_af_rejects = max_af_rejects;
	pairwise_edge_list_free(&edges);
}

void pairwise_build_complete_linkage_dedup(int n, double cut, uint32_t ctxcut,
										   double max_afcut, void *cb_ctx,
										   pairwise_eval_pair_fn eval_pair,
										   pairwise_quality_compare_fn quality_compare,
										   pairwise_edge_observer_fn edge_observer,
										   void *edge_ctx,
										   pairwise_progress_fn progress_fn,
										   void *progress_ctx,
										   pairwise_component_result_t *result)
{
	if (!result || n < 0 || !eval_pair || !quality_compare)
		errx(EINVAL, "%s(): invalid full-linkage dedup arguments", __func__);

	pairwise_edge_list_t edges = {0};
	uint64_t distance_edges = 0;
	uint64_t duplicate_edges = 0;
	uint64_t ctx_rejects = 0;
	uint64_t max_af_rejects = 0;
	pairwise_collect_dedup_edges(n, cut, ctxcut, max_afcut, cb_ctx, eval_pair,
								 edge_observer, edge_ctx,
								 progress_fn, progress_ctx, &edges,
								 &distance_edges, &duplicate_edges,
								 &ctx_rejects, &max_af_rejects);

	pairwise_complete_linkage_dedup_from_edges(n, edges.edges, edges.n,
											   cb_ctx, quality_compare, result);
	result->distance_edges = distance_edges;
	result->duplicate_edges = duplicate_edges;
	result->ctx_rejects = ctx_rejects;
	result->max_af_rejects = max_af_rejects;
	pairwise_edge_list_free(&edges);
}

static size_t pairwise_ctxgid_upper_bound(const ctxgidobj_t *arr, size_t n,
										  uint64_t ctx, size_t low)
{
	size_t high = n;
	while (low < high) {
		const size_t mid = low + (high - low) / 2;
		if ((arr[mid].ctxgid >> GID_NBITS) <= ctx)
			low = mid + 1;
		else
			high = mid;
	}
	return low;
}

void pairwise_indexed_self_scan(const unify_sketch_t *sketch,
								const char *sorted_index_path,
								const pairwise_index_scan_options_t *opt,
								pairwise_index_edge_fn edge_fn,
								void *edge_ctx,
								pairwise_progress_fn progress_fn,
								void *progress_ctx,
								pairwise_index_scan_stats_t *stats)
{
	if (!sketch || !sorted_index_path || !opt)
		errx(EINVAL, "%s(): invalid indexed scan arguments", __func__);
	if (sketch->stat_type != 2)
		errx(EINVAL, "%s(): indexed scan requires lco sketches", __func__);
	if (sketch->conflict)
		errx(EINVAL, "%s(): indexed scan does not support conflict sketches yet", __func__);

	uint64_t candidate_pairs = 0;
	uint64_t exact_pairs = 0;
	uint64_t distance_edges = 0;
	uint64_t accepted_edges = 0;
	uint64_t ctx_rejects = 0;
	uint64_t max_af_rejects = 0;
	pairwise_metric_expr_t default_metric = pairwise_metric_expr_single(PAIRWISE_METRIC_CTX_NAIVE);
	const pairwise_metric_expr_t *metric = opt->metric ? opt->metric : &default_metric;
	const uint32_t index_min_votes = opt->index_min_votes > 0 ? opt->index_min_votes : 1;
	const uint32_t index_sample_step = opt->index_sample_step > 0 ? opt->index_sample_step : 1;
	const int threads = opt->threads > 0 ? opt->threads : 1;

	const int n = sketch->infile_num;
	const uint64_t total_entries = sketch->sketch_index[n];
	size_t index_size = 0;
	ctxgidobj_t *sorted_index = read_from_file(sorted_index_path, &index_size);
	if (index_size != (size_t)total_entries * sizeof(sorted_index[0]))
		err(EINVAL, "%s(): %s has %zu bytes, expected %zu", __func__,
			sorted_index_path, index_size,
			(size_t)total_entries * sizeof(sorted_index[0]));

	pairwise_prepare_lco_model(sketch);
	const uint64_t gidmask = (1ULL << GID_NBITS) - 1ULL;
	const uint8_t nobjbits = Bitslen.obj;
	const int fence_k = kssd_choose_k_fenceposts((size_t)total_entries, 0);
	const size_t fencepost_count = ((size_t)1 << fence_k) + 1;
	size_t *fenceposts = malloc(fencepost_count * sizeof(fenceposts[0]));
	if (!fenceposts)
		err(EXIT_FAILURE, "%s(): OOM sorted-index fenceposts", __func__);
	if (kssd_build_fenceposts_ctxgid(sorted_index, (size_t)total_entries,
									 fence_k, fenceposts) != 0)
		errx(EXIT_FAILURE, "%s(): failed to build sorted-index fenceposts", __func__);

	uint64_t completed_rows = 0;
#pragma omp parallel num_threads(threads) reduction(+:candidate_pairs,exact_pairs,distance_edges,accepted_edges,ctx_rejects,max_af_rejects)
	{
		uint32_t *votes = calloc((size_t)n, sizeof(votes[0]));
		uint32_t *touched = malloc((size_t)n * sizeof(touched[0]));
		if (!votes || !touched)
			err(EXIT_FAILURE, "%s(): OOM indexed scan work buffers", __func__);

#pragma omp for schedule(dynamic, 1)
		for (int qn = 0; qn < n; ++qn) {
			const uint64_t q_begin = sketch->sketch_index[qn];
			const uint64_t q_end = sketch->sketch_index[qn + 1];
			const uint64_t *qry = sketch->comb_sketch + q_begin;
			const size_t qry_n = (size_t)(q_end - q_begin);
			if (qry_n == 0) {
				if (progress_fn) {
					uint64_t done;
#pragma omp atomic capture
					done = ++completed_rows;
#pragma omp critical(pairwise_index_progress_callback)
					progress_fn(progress_ctx, done, (uint64_t)n);
				}
				continue;
			}

			const size_t lookup_n = (qry_n + index_sample_step - 1) / index_sample_step;
			uint64_t *sampled_qry = NULL;
			const uint64_t *lookup_qry = qry;
			if (index_sample_step > 1) {
				sampled_qry = malloc(lookup_n * sizeof(sampled_qry[0]));
				if (!sampled_qry)
					err(EXIT_FAILURE, "%s(): OOM sampled query contexts", __func__);
				for (size_t si = 0; si < lookup_n; ++si)
					sampled_qry[si] = qry[si * index_sample_step];
				lookup_qry = sampled_qry;
			}

			size_t *first = kssd_find_first_occurrences_fenceposts(lookup_qry, lookup_n,
																   sorted_index,
																   (size_t)total_entries,
																   fenceposts,
																   fence_k,
																   nobjbits);
			if (!first)
				err(EXIT_FAILURE, "%s(): OOM index lookup buffer", __func__);
			size_t touched_n = 0;

			for (size_t qi = 0; qi < lookup_n; ++qi) {
				if (first[qi] == SIZE_MAX)
					continue;
				const uint64_t qry_ctx = lookup_qry[qi] >> nobjbits;
				const size_t ctx_end =
					pairwise_ctxgid_upper_bound(sorted_index, (size_t)total_entries,
												qry_ctx, first[qi]);
				if (opt->index_max_ctx_freq > 0 &&
					ctx_end - first[qi] > opt->index_max_ctx_freq)
					continue;
				for (size_t d = first[qi]; d < ctx_end; ++d) {
					const uint32_t gid = (uint32_t)(sorted_index[d].ctxgid & gidmask);
					if (gid >= (uint32_t)qn)
						break;
					if (votes[gid] == 0)
						touched[touched_n++] = gid;
					if (votes[gid] < UINT32_MAX)
						++votes[gid];
				}
			}
			free(first);
			free(sampled_qry);

			for (size_t ti = 0; ti < touched_n; ++ti) {
				const uint32_t rn = touched[ti];
				++candidate_pairs;
				if (votes[rn] < index_min_votes) {
					votes[rn] = 0;
					continue;
				}
				++exact_pairs;
				const pairwise_eval_t eval =
					pairwise_eval_expr_samples(metric, sketch, (uint32_t)qn,
											   sketch, rn, false);
				if (eval.valid && eval.distance < opt->cut) {
					++distance_edges;
					if (eval.xny_ctx < opt->ctxcut) {
						++ctx_rejects;
					} else if (eval.max_af < opt->max_afcut) {
						++max_af_rejects;
					} else {
						++accepted_edges;
						if (edge_fn) {
#pragma omp critical(pairwise_index_edge_callback)
							edge_fn(edge_ctx, qn, (int)rn, &eval);
						}
					}
				}
				votes[rn] = 0;
			}
			if (progress_fn) {
				uint64_t done;
#pragma omp atomic capture
				done = ++completed_rows;
#pragma omp critical(pairwise_index_progress_callback)
				progress_fn(progress_ctx, done, (uint64_t)n);
			}
		}
		free(votes);
		free(touched);
	}

	free(fenceposts);
	free_read_from_file(sorted_index, index_size);
	if (stats) {
		stats->candidate_pairs = candidate_pairs;
		stats->exact_pairs = exact_pairs;
		stats->distance_edges = distance_edges;
		stats->accepted_edges = accepted_edges;
		stats->ctx_rejects = ctx_rejects;
		stats->max_af_rejects = max_af_rejects;
	}
}
