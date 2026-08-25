#include "command_matrix.h"
#include "global_basic.h"
#include "command_progress.h"
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <err.h>
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <libgen.h>
#include <dirent.h>
#include <omp.h>
#include <stdatomic.h>
#include <ctype.h>
#include "../klib/khash.h"
#include "../klib/kstring.h"
#include "kssdlib_sort.h"
#include "command_ani.h"

extern const char sorted_comb_ctxgid64obj32[];
// core functions
double get_mashD(uint32_t K, uint32_t X, uint32_t Y, uint32_t XnY)
{
	return (-log(2 * JCD(X, Y, XnY) / (1 + JCD(X, Y, XnY))) / (K));
}
double get_aafD(uint32_t K, uint32_t X, uint32_t Y, uint32_t XnY)
{
	return (-log(CTM(X, Y, XnY)) / (K));
}
typedef double (*Dist)(uint32_t, uint32_t, uint32_t, uint32_t);

#define MATRIX_PROGRESS_AUTO_MIN 1000ULL

static bool matrix_progress_should_enable(const matrix_opt_t *opt, uint64_t total)
{
	if (opt->outf[0] == '\0' || opt->progress_mode == MATRIX_PROGRESS_OFF)
		return false;
	if (opt->progress_mode == MATRIX_PROGRESS_ON)
		return true;
	return total >= MATRIX_PROGRESS_AUTO_MIN;
}

static kssd_progress_t matrix_progress_start(const matrix_opt_t *opt,
											 const char *label,
											 const char *unit,
											 uint64_t total)
{
	return kssd_progress_start(matrix_progress_should_enable(opt, total),
							   "matrix", label, unit, total);
}

static FILE *matrix_open_output(const char *path)
{
	FILE *output = path[0] == '\0' ? stdout : fopen(path, "w");
	if (output == NULL)
		err(errno, "%s(): %s", __func__, path);
	return output;
}

static void matrix_make_phylip_id(int one_based_index, char *buf, size_t buf_size);
static void matrix_write_full_matrix_idmap(const matrix_opt_t *opt,
										   const unify_sketch_t *sketch);

static matrix_format_t matrix_effective_format(const matrix_opt_t *opt, bool rectangular)
{
	if (opt->format != MATRIX_FORMAT_AUTO)
		return opt->format;
	return rectangular ? MATRIX_FORMAT_FULL : MATRIX_FORMAT_TRIANGLE;
}

static unsigned matrix_parse_flags_for_format(matrix_format_t fmt)
{
	return (fmt == MATRIX_FORMAT_CLUSTERS || fmt == MATRIX_FORMAT_DEDUP_PLAN)
			   ? SKETCH_PARSE_INFILE_META
			   : SKETCH_PARSE_NONE;
}

static double matrix_eval_value(const matrix_opt_t *opt, const pairwise_eval_t *eval)
{
	return eval->valid ? eval->distance : opt->e;
}

static void matrix_print_edge_header(FILE *out)
{
	fprintf(out, "Qry\tRef\tDistance\tSimilarity\tMetric\tXnY_ctx\tQry_align_fraction\tRef_align_fraction\tmax_AF\tN_diff_obj\tN_diff_obj_section\tN_mut2_ctx\n");
}

static void matrix_print_edge(FILE *out, const char *qry_name, const char *ref_name,
							  const pairwise_eval_t *eval,
							  const pairwise_metric_expr_t *metric)
{
	fprintf(out, "%s\t%s\t%.12g\t%.12g\t%s\t%u\t%.12g\t%.12g\t%.12g\t%u\t%u\t%u\n",
			qry_name, ref_name, eval->distance, eval->similarity,
			pairwise_metric_expr_name(metric), eval->xny_ctx, eval->af_qry,
			eval->af_ref, eval->max_af, eval->n_diff_obj,
			eval->n_diff_obj_section, eval->n_mut2_ctx);
}

static void matrix_kput_fixed6(kstring_t *s, double value)
{
	if (!isfinite(value) || value > 900000000000.0 || value < -900000000000.0) {
		ksprintf(s, "%lf", value);
		return;
	}
	if (value < 0.0) {
		kputc('-', s);
		value = -value;
	}
	const uint64_t scaled = (uint64_t)(value * 1000000.0 + 0.5);
	const uint64_t whole = scaled / 1000000ULL;
	uint32_t frac = (uint32_t)(scaled % 1000000ULL);
	kputl((long)whole, s);
	kputc('.', s);
	char frac_buf[6];
	for (int i = 5; i >= 0; --i) {
		frac_buf[i] = (char)('0' + (frac % 10));
		frac /= 10;
	}
	kputsn(frac_buf, 6, s);
}

static uint32_t matrix_scale_fixed6(double value)
{
	if (!isfinite(value))
		value = 1.0;
	if (value < 0.0)
		value = 0.0;
	if (value > 4294.0)
		value = 4294.0;
	return (uint32_t)(value * 1000000.0 + 0.5);
}

static void matrix_kput_scaled6(kstring_t *s, uint32_t scaled)
{
	const uint32_t whole = scaled / 1000000U;
	uint32_t frac = scaled % 1000000U;
	kputuw(whole, s);
	kputc('.', s);
	char frac_buf[6];
	for (int i = 5; i >= 0; --i) {
		frac_buf[i] = (char)('0' + (frac % 10));
		frac /= 10;
	}
	kputsn(frac_buf, 6, s);
}

static inline size_t matrix_lower_offset(int row, int col)
{
	return (size_t)row * (size_t)(row - 1) / 2u + (size_t)col;
}

static bool matrix_meta_valid(const unify_sketch_t *sketch, int idx)
{
	if (!sketch || !sketch->infile_meta || idx < 0 || idx >= sketch->infile_num)
		return false;
	const infile_meta_t *meta = &sketch->infile_meta[idx];
	return meta->meta_fmt_version == KSSD_INFILE_META_VERSION &&
		   meta->total_length_bp > 0;
}

static const char *matrix_infile_fmt_name(int8_t fmt)
{
	switch (fmt) {
	case KSSD_INFILE_FMT_FASTA:
		return "FASTA";
	case KSSD_INFILE_FMT_FASTQ:
		return "FASTQ";
	default:
		return "NA";
	}
}

static uint64_t matrix_sketch_entries(const unify_sketch_t *sketch, int idx)
{
	return (uint64_t)(sketch->sketch_index[idx + 1] - sketch->sketch_index[idx]);
}

static const char *matrix_rep_rule(const unify_sketch_t *sketch, int sample, int rep)
{
	if (sample == rep)
		return "self";

	const bool sample_meta_valid = matrix_meta_valid(sketch, sample);
	const bool rep_meta_valid = matrix_meta_valid(sketch, rep);
	if (sample_meta_valid != rep_meta_valid)
		return rep_meta_valid ? "valid_metadata" : "sample_valid_metadata";

	if (sample_meta_valid && rep_meta_valid) {
		const infile_meta_t *sample_meta = &sketch->infile_meta[sample];
		const infile_meta_t *rep_meta = &sketch->infile_meta[rep];
		const bool sample_fasta = sample_meta->infile_fmt == KSSD_INFILE_FMT_FASTA;
		const bool rep_fasta = rep_meta->infile_fmt == KSSD_INFILE_FMT_FASTA;
		if (sample_fasta != rep_fasta)
			return rep_fasta ? "fasta_input" : "sample_fasta_input";
		if (sample_meta->asm_level != rep_meta->asm_level)
			return rep_meta->asm_level > sample_meta->asm_level
					   ? "higher_asm_level"
					   : "sample_higher_asm_level";
		if (sample_meta->total_length_bp != rep_meta->total_length_bp)
			return rep_meta->total_length_bp > sample_meta->total_length_bp
					   ? "longer_total_length"
					   : "sample_longer_total_length";
		if (sample_meta->record_count != rep_meta->record_count) {
			if (sample_meta->record_count == 0)
				return "sample_missing_record_count";
			if (rep_meta->record_count == 0)
				return "rep_missing_record_count";
			return rep_meta->record_count < sample_meta->record_count
					   ? "fewer_records"
					   : "sample_fewer_records";
		}
		if (sample_meta->median_length_bp != rep_meta->median_length_bp)
			return rep_meta->median_length_bp > sample_meta->median_length_bp
					   ? "longer_median_length"
					   : "sample_longer_median_length";
		if (sample_meta->length_cv != rep_meta->length_cv)
			return rep_meta->length_cv < sample_meta->length_cv
					   ? "lower_length_cv"
					   : "sample_lower_length_cv";
	}

	const uint64_t sample_entries = matrix_sketch_entries(sketch, sample);
	const uint64_t rep_entries = matrix_sketch_entries(sketch, rep);
	if (sample_entries != rep_entries)
		return rep_entries > sample_entries ? "larger_sketch_entries" : "sample_larger_sketch_entries";
	return rep < sample ? "earlier_sample_order" : "sample_earlier_sample_order";
}

static void matrix_print_meta_field(FILE *output, const unify_sketch_t *sketch, int idx)
{
	if (!matrix_meta_valid(sketch, idx)) {
		fprintf(output, "\tNA\tNA\tNA\tNA\tNA\tNA");
		return;
	}
	const infile_meta_t *meta = &sketch->infile_meta[idx];
	fprintf(output, "\t%s\t%.6g\t%" PRIu64 "\t%u\t%u\t%.6g",
			matrix_infile_fmt_name(meta->infile_fmt), (double)meta->asm_level,
			meta->total_length_bp, meta->record_count,
			meta->median_length_bp, (double)meta->length_cv);
}

static void matrix_write_glist(const matrix_opt_t *opt, const unify_sketch_t *sketch)
{
	if (opt->gl[0] == '\0')
		return;
	FILE *glout = fopen(opt->gl, "w");
	if (!glout)
		err(errno, "%s(): Failed to open file:%s\n", __func__, opt->gl);
	for (int i = 0; i < sketch->infile_num; i++)
		fprintf(glout, "%d\t%s\n", i + 1, sketch->gname[i]);
	fclose(glout);
}

static void matrix_compute_dense(const matrix_opt_t *opt, const unify_sketch_t *ref,
								 const unify_sketch_t *qry, bool same_sketch)
{
	if (opt->matrix_out_format == MATRIX_KEEP_MATRIX_PHYLIP &&
		(!same_sketch || ref->infile_num != qry->infile_num))
		errx(EINVAL, "%s(): --matrix-format phylip requires one square sketch", __func__);
	FILE *output = matrix_open_output(opt->outf);
	if (opt->matrix_out_format == MATRIX_KEEP_MATRIX_PHYLIP)
		matrix_write_full_matrix_idmap(opt, ref);
	kssd_progress_t progress =
		matrix_progress_start(opt, "full rows", "rows",
							  (uint64_t)(opt->matrix_out_format == MATRIX_KEEP_MATRIX_PHYLIP
											 ? ref->infile_num
											 : qry->infile_num));
	if (opt->matrix_out_format == MATRIX_KEEP_MATRIX_PHYLIP) {
		fprintf(output, "%d\n", ref->infile_num);
		for (int rn = 0; rn < ref->infile_num; rn++) {
			char id[16];
			matrix_make_phylip_id(rn + 1, id, sizeof(id));
			fprintf(output, "%-10s", id);
			for (int qn = 0; qn < qry->infile_num; qn++) {
				double value = opt->diagonal_value;
				if (!(same_sketch && rn == qn)) {
					pairwise_eval_t eval =
						pairwise_eval_expr_samples(&opt->metric, qry, (uint32_t)qn,
												   ref, (uint32_t)rn, false);
					value = matrix_eval_value(opt, &eval);
				}
				fprintf(output, " %lf", value);
			}
			fprintf(output, "\n");
			kssd_progress_update(&progress, (uint64_t)rn + 1, false);
		}
	} else {
		for (int rn = 0; rn < ref->infile_num; rn++)
			fprintf(output, "\t%s", ref->gname[rn]);
		fprintf(output, "\n");

		for (int qn = 0; qn < qry->infile_num; qn++) {
			fprintf(output, "%s", qry->gname[qn]);
			for (int rn = 0; rn < ref->infile_num; rn++) {
				double value = opt->diagonal_value;
				if (!(same_sketch && rn == qn)) {
					pairwise_eval_t eval =
						pairwise_eval_expr_samples(&opt->metric, qry, (uint32_t)qn,
												   ref, (uint32_t)rn, false);
					value = matrix_eval_value(opt, &eval);
				}
				fprintf(output, "\t%lf", value);
			}
			fprintf(output, "\n");
			kssd_progress_update(&progress, (uint64_t)qn + 1, false);
		}
	}
	kssd_progress_done(&progress);
	if (output != stdout)
		fclose(output);
}

static int *matrix_collect_kept_indices(const pairwise_component_result_t *comp,
										int *kept_n_out)
{
	if (!comp || !kept_n_out)
		errx(EINVAL, "%s(): invalid kept-index arguments", __func__);
	int kept_n = 0;
	for (int i = 0; i < comp->n; ++i) {
		if (!comp->remove_sample || !comp->remove_sample[i])
			++kept_n;
	}
	int *kept = kept_n > 0 ? malloc((size_t)kept_n * sizeof(kept[0])) : NULL;
	if (kept_n > 0 && !kept)
		err(EXIT_FAILURE, "%s(): OOM kept sample list", __func__);
	int k = 0;
	for (int i = 0; i < comp->n; ++i) {
		if (!comp->remove_sample || !comp->remove_sample[i])
			kept[k++] = i;
	}
	*kept_n_out = kept_n;
	return kept;
}

static void matrix_make_phylip_id(int one_based_index, char *buf, size_t buf_size)
{
	if (one_based_index <= 0 || one_based_index > 999999999)
		errx(EINVAL, "%s(): PHYLIP short id index out of range", __func__);
	if (snprintf(buf, buf_size, "S%09d", one_based_index) >= (int)buf_size)
		errx(EINVAL, "%s(): PHYLIP sample id overflow", __func__);
}

static void matrix_write_keep_matrix_idmap(const matrix_opt_t *opt,
										   const unify_sketch_t *sketch,
										   const int *kept,
										   int kept_n)
{
	if (opt->keep_matrix_idmap_outf[0] == '\0')
		return;
	FILE *idmap = matrix_open_output(opt->keep_matrix_idmap_outf);
	fprintf(idmap, "id\tsample\n");
	for (int i = 0; i < kept_n; ++i) {
		char id[16];
		matrix_make_phylip_id(i + 1, id, sizeof(id));
		fprintf(idmap, "%s\t%s\n", id, sketch->gname[kept[i]]);
	}
	if (idmap != stdout)
		fclose(idmap);
}

static void matrix_write_full_matrix_idmap(const matrix_opt_t *opt,
										   const unify_sketch_t *sketch)
{
	if (opt->matrix_idmap_outf[0] == '\0')
		return;
	FILE *idmap = matrix_open_output(opt->matrix_idmap_outf);
	fprintf(idmap, "id\tsample\n");
	for (int i = 0; i < sketch->infile_num; ++i) {
		char id[16];
		matrix_make_phylip_id(i + 1, id, sizeof(id));
		fprintf(idmap, "%s\t%s\n", id, sketch->gname[i]);
	}
	if (idmap != stdout)
		fclose(idmap);
}

static bool matrix_metric_can_use_context_index(const matrix_opt_t *opt)
{
	if (!opt || opt->metric.count == 0)
		return false;
	for (size_t i = 0; i < opt->metric.count; ++i) {
		switch (opt->metric.metrics[i]) {
		case PAIRWISE_METRIC_CTX_NAIVE:
		case PAIRWISE_METRIC_CTX_MOE:
		case PAIRWISE_METRIC_MASH:
		case PAIRWISE_METRIC_AAF:
			break;
		default:
			return false;
		}
	}
	return true;
}

static double matrix_context_metric_distance(pairwise_metric_t metric,
											 ani_features_t *features,
											 uint32_t qry_ctx,
											 uint32_t ref_ctx)
{
	switch (metric) {
	case PAIRWISE_METRIC_CTX_NAIVE:
		return get_naive_dist(features);
	case PAIRWISE_METRIC_CTX_MOE:
		return lm3ways_dist_from_features(features);
	case PAIRWISE_METRIC_MASH:
		return get_mashD(Bitslen.ctx / 2, ref_ctx, qry_ctx, features->XnY_ctx);
	case PAIRWISE_METRIC_AAF:
		return get_aafD(Bitslen.ctx / 2, ref_ctx, qry_ctx, features->XnY_ctx);
	}
	return 1.0;
}

static double matrix_context_feature_value(const matrix_opt_t *opt,
										   const ani_features_t *features,
										   uint32_t qry_ctx,
										   uint32_t ref_ctx);
static int matrix_choose_context_block_size(int ref_n, int qry_n);
static int matrix_choose_dense_context_block_size(int ref_n, int qry_n);
static uint32_t *matrix_context_counts_for_sketch(const unify_sketch_t *sketch,
												  bool ignoreconflict);
static ctxgidobj_t *matrix_load_sorted_index(const unify_sketch_t *sketch,
											 const char *sketch_dir,
											 size_t *index_bytes,
											 bool *from_file);
static void matrix_free_sorted_index(ctxgidobj_t *index, size_t index_bytes,
									 bool from_file);
static void matrix_fill_query_ctx_counts(uint32_t *ctx_counts,
										 const unify_sketch_t *sketch,
										 int offset_gid,
										 int block_size);

static bool matrix_write_indexed_self_matrix(const matrix_opt_t *opt,
											 const unify_sketch_t *sketch,
											 const char *sketch_dir,
											 matrix_format_t format)
{
	if (!matrix_metric_can_use_context_index(opt) || sketch->stat_type != 2)
		return false;
	if (format != MATRIX_FORMAT_FULL && format != MATRIX_FORMAT_TRIANGLE)
		return false;

	const int n = sketch->infile_num;
	size_t index_bytes = 0;
	bool index_from_file = false;
	ctxgidobj_t *sorted_index =
		matrix_load_sorted_index(sketch, sketch_dir, &index_bytes, &index_from_file);
	const size_t ref_sksize = (size_t)sketch->sketch_index[n];

	uint32_t *ctx_counts = matrix_context_counts_for_sketch(sketch, false);
	const int block_size = matrix_choose_dense_context_block_size(n, n);
	if (block_size <= 0)
		errx(EINVAL, "%s(): invalid matrix block size", __func__);
	ctx_mut2_t *ctx = malloc((size_t)n * (size_t)block_size * sizeof(ctx[0]));
	obj_section_t *obj = malloc((size_t)n * (size_t)block_size * sizeof(obj[0]));
	if (!ctx || !obj)
		err(EXIT_FAILURE, "%s(): OOM indexed dense matrix block", __func__);

	ani_opt_t ani_opt = {0};
	ani_opt.p = opt->p > 0 ? opt->p : 1;
	ani_opt.fmt = 1;
	ani_opt.ntop = -1;
	ani_opt.e = (int)opt->e;
	ani_opt.s = 4;

	FILE *output = matrix_open_output(opt->outf);
	kstring_t row = {0, 0, NULL};
	if (format == MATRIX_FORMAT_FULL) {
		const size_t lower_n = (size_t)n * (size_t)(n - 1) / 2u;
		uint32_t *lower = lower_n > 0 ? malloc(lower_n * sizeof(lower[0])) : NULL;
		if (lower_n > 0 && !lower)
			err(EXIT_FAILURE, "%s(): OOM indexed full lower matrix", __func__);

		kssd_progress_t compute_progress =
			matrix_progress_start(opt, "indexed full lower pairs", "pairs",
								  (uint64_t)lower_n);
		for (int offset = 0; offset < n; offset += block_size) {
			const int this_block = (offset + block_size > n) ? n - offset : block_size;
			matrix_fill_query_ctx_counts(ctx_counts, sketch, offset, this_block);
			memset(ctx, 0, (size_t)n * (size_t)this_block * sizeof(ctx[0]));
			memset(obj, 0, (size_t)n * (size_t)this_block * sizeof(obj[0]));
			count_ctx_obj_frm_comb_sketch_section_lower(
				ctx, obj, sorted_index, ref_sksize, n, offset, this_block,
				sketch->comb_sketch + sketch->sketch_index[offset],
				sketch->sketch_index + offset, &ani_opt);

			for (int i = 0; i < this_block; ++i) {
				const int qgid = offset + i;
				for (int rgid = 0; rgid < qgid; ++rgid) {
					ani_features_t features = {
						.XnY_ctx = MCTX(n, i, rgid).num_ctx,
						.N_diff_obj = MOBJ(n, i, rgid).diff_obj,
						.N_diff_obj_section = MOBJ(n, i, rgid).diff_obj_section,
						.N_mut2_ctx = MCTX(n, i, rgid).num_mut2_ctx,
					};
					const double value = matrix_context_feature_value(
						opt, &features, ctx_counts[qgid], ctx_counts[rgid]);
					lower[matrix_lower_offset(qgid, rgid)] = matrix_scale_fixed6(value);
				}
				kssd_progress_update(&compute_progress,
									   (uint64_t)matrix_lower_offset(qgid + 1, 0),
									   false);
			}
		}
		kssd_progress_done(&compute_progress);

		if (opt->matrix_out_format == MATRIX_KEEP_MATRIX_PHYLIP) {
			matrix_write_full_matrix_idmap(opt, sketch);
			ksprintf(&row, "%d\n", n);
		} else {
			for (int q = 0; q < n; ++q) {
				kputc('\t', &row);
				kputs(sketch->gname[q], &row);
			}
			kputc('\n', &row);
		}
		fwrite(row.s, 1, row.l, output);

		const uint32_t diag_scaled = matrix_scale_fixed6(opt->diagonal_value);
		kssd_progress_t write_progress =
			matrix_progress_start(opt, "indexed full write rows", "rows",
								  (uint64_t)n);
		for (int qgid = 0; qgid < n; ++qgid) {
			row.l = 0;
			if (opt->matrix_out_format == MATRIX_KEEP_MATRIX_PHYLIP) {
				char id[16];
				matrix_make_phylip_id(qgid + 1, id, sizeof(id));
				ksprintf(&row, "%-10s", id);
			} else {
				kputs(sketch->gname[qgid], &row);
			}
			for (int rgid = 0; rgid < n; ++rgid) {
				uint32_t scaled = diag_scaled;
				if (qgid > rgid)
					scaled = lower[matrix_lower_offset(qgid, rgid)];
				else if (qgid < rgid)
					scaled = lower[matrix_lower_offset(rgid, qgid)];
				kputc(opt->matrix_out_format == MATRIX_KEEP_MATRIX_PHYLIP ? ' ' : '\t',
					  &row);
				matrix_kput_scaled6(&row, scaled);
			}
			kputc('\n', &row);
			fwrite(row.s, 1, row.l, output);
			kssd_progress_update(&write_progress, (uint64_t)qgid + 1, false);
		}
		kssd_progress_done(&write_progress);

		if (output != stdout)
			fclose(output);
		free(lower);
		free(row.s);
		free(ctx);
		free(obj);
		free(ctx_counts);
		matrix_free_sorted_index(sorted_index, index_bytes, index_from_file);
		return true;
	}
	const uint64_t lower_n = (uint64_t)n * (uint64_t)(n - 1) / 2u;
	kssd_progress_t progress =
		matrix_progress_start(opt, "indexed triangle lower pairs", "pairs", lower_n);
	for (int offset = 0; offset < n; offset += block_size) {
		const int this_block = (offset + block_size > n) ? n - offset : block_size;
		matrix_fill_query_ctx_counts(ctx_counts, sketch, offset, this_block);
		memset(ctx, 0, (size_t)n * (size_t)this_block * sizeof(ctx[0]));
		memset(obj, 0, (size_t)n * (size_t)this_block * sizeof(obj[0]));
		count_ctx_obj_frm_comb_sketch_section_lower(
			ctx, obj, sorted_index, ref_sksize, n, offset, this_block,
			sketch->comb_sketch + sketch->sketch_index[offset],
			sketch->sketch_index + offset, &ani_opt);

		for (int i = 0; i < this_block; ++i) {
			const int qgid = offset + i;
			row.l = 0;
			kputs(sketch->gname[qgid], &row);
			const int col_end = format == MATRIX_FORMAT_TRIANGLE ? qgid : n;
			for (int rgid = 0; rgid < col_end; ++rgid) {
				double value = opt->diagonal_value;
				if (qgid != rgid) {
					ani_features_t features = {
						.XnY_ctx = MCTX(n, i, rgid).num_ctx,
						.N_diff_obj = MOBJ(n, i, rgid).diff_obj,
						.N_diff_obj_section = MOBJ(n, i, rgid).diff_obj_section,
						.N_mut2_ctx = MCTX(n, i, rgid).num_mut2_ctx,
					};
					value = matrix_context_feature_value(
						opt, &features, ctx_counts[qgid], ctx_counts[rgid]);
				}
				kputc('\t', &row);
				matrix_kput_fixed6(&row, value);
			}
			if (format == MATRIX_FORMAT_TRIANGLE && opt->d) {
				kputc('\t', &row);
				matrix_kput_fixed6(&row, opt->diagonal_value);
			}
			kputc('\n', &row);
			fwrite(row.s, 1, row.l, output);
			kssd_progress_update(&progress,
								   (uint64_t)matrix_lower_offset(qgid + 1, 0),
								   false);
		}
	}
	kssd_progress_done(&progress);

	if (output != stdout)
		fclose(output);
	free(row.s);
	free(ctx);
	free(obj);
	free(ctx_counts);
	matrix_free_sorted_index(sorted_index, index_bytes, index_from_file);
	return true;
}

static double matrix_context_feature_value(const matrix_opt_t *opt,
										   const ani_features_t *features,
										   uint32_t qry_ctx,
										   uint32_t ref_ctx)
{
	if (!features || features->XnY_ctx == 0 || qry_ctx == 0 || ref_ctx == 0)
		return opt->e;
	bool have_value = false;
	double selected = opt->metric.op == PAIRWISE_METRIC_EXPR_AND ? 0.0 : 1.0;
	for (size_t i = 0; i < opt->metric.count; ++i) {
		ani_features_t tmp = *features;
		double distance =
			matrix_context_metric_distance(opt->metric.metrics[i], &tmp, qry_ctx, ref_ctx);
		if (!isfinite(distance))
			return opt->e;
		if (distance < 0.0)
			distance = 0.0;
		if (distance > 1.0)
			distance = 1.0;
		if (!have_value) {
			selected = distance;
			have_value = true;
		} else if (opt->metric.op == PAIRWISE_METRIC_EXPR_AND) {
			if (distance > selected)
				selected = distance;
		} else {
			if (distance < selected)
				selected = distance;
		}
	}
	if (!have_value)
		return opt->e;
	double distance = selected;
	if (!isfinite(distance))
		return opt->e;
	if (distance < 0.0)
		distance = 0.0;
	if (distance > 1.0)
		distance = 1.0;
	return distance;
}

static int matrix_choose_context_block_size(int ref_n, int qry_n)
{
	if (ref_n <= 0 || qry_n <= 0)
		return 0;
	const size_t bytes_per_query =
		(size_t)ref_n * (sizeof(ctx_mut2_t) + sizeof(obj_section_t));
	if (bytes_per_query == 0)
		return 0;
	size_t budget = GetAvailableMemory() / 4;
	if (budget < bytes_per_query * 64)
		budget = bytes_per_query * 64;
	size_t block = budget / bytes_per_query;
	if (block < 64)
		block = 64;
	if (block > 4096)
		block = 4096;
	if (block > (size_t)qry_n)
		block = (size_t)qry_n;
	return (int)block;
}

static int matrix_choose_dense_context_block_size(int ref_n, int qry_n)
{
	int block = matrix_choose_context_block_size(ref_n, qry_n);
	if (block > 64)
		block = 64;
	const char *env = getenv("KSSD3A_MATRIX_BLOCK_SIZE");
	if (env && env[0] != '\0') {
		char *end = NULL;
		errno = 0;
		long requested = strtol(env, &end, 10);
		if (errno == 0 && end != env && *end == '\0' && requested > 0) {
			if (requested > qry_n)
				requested = qry_n;
			block = (int)requested;
		}
	}
	return block;
}

static uint32_t *matrix_context_counts_for_sketch(const unify_sketch_t *sketch,
												  bool ignoreconflict)
{
	uint32_t *counts = calloc((size_t)sketch->infile_num, sizeof(counts[0]));
	if (!counts)
		err(EXIT_FAILURE, "%s(): OOM context counts", __func__);
	for (int i = 0; i < sketch->infile_num; ++i) {
		const uint64_t begin = sketch->sketch_index[i];
		const uint64_t end = sketch->sketch_index[i + 1];
		const uint64_t *arr = sketch->comb_sketch + begin;
		const size_t len = (size_t)(end - begin);
		counts[i] = sketch->conflict
						? pairwise_count_ctx_runs_sorted_ctxobj64(arr, len)
						: (uint32_t)len;
		if (ignoreconflict && sketch->conflict) {
			uint32_t nonconflict = 0;
			const uint8_t nobjbits = Bitslen.obj;
			for (size_t p = 0; p < len;) {
				const uint64_t ctx = arr[p] >> nobjbits;
				const size_t first = p;
				do {
					++p;
				} while (p < len && (arr[p] >> nobjbits) == ctx);
				if (p - first == 1)
					++nonconflict;
			}
			counts[i] = nonconflict;
		}
	}
	return counts;
}

static ctxgidobj_t *matrix_load_sorted_index(const unify_sketch_t *sketch,
											 const char *sketch_dir,
											 size_t *index_bytes,
											 bool *from_file)
{
	*index_bytes = 0;
	*from_file = false;
	const size_t total_entries = (size_t)sketch->sketch_index[sketch->infile_num];
	if (sketch_dir && sketch_dir[0] != '\0' &&
		file_exists_in_folder((char *)sketch_dir, (char *)sorted_comb_ctxgid64obj32)) {
		char *index_path = test_get_fullpath(sketch_dir, sorted_comb_ctxgid64obj32);
		ctxgidobj_t *index = read_from_file(index_path, index_bytes);
		free(index_path);
		if (*index_bytes != total_entries * sizeof(index[0]))
			errx(EINVAL, "%s(): sorted index size mismatch", __func__);
		*from_file = true;
		return index;
	}
	*index_bytes = total_entries * sizeof(ctxgidobj_t);
	return comb_sortedsketch64_2sortedcomb_ctxgid64obj32((unify_sketch_t *)sketch);
}

static void matrix_free_sorted_index(ctxgidobj_t *index, size_t index_bytes,
									 bool from_file)
{
	if (!index)
		return;
	if (from_file)
		free_read_from_file(index, index_bytes);
	else
		free(index);
}

static void matrix_fill_query_ctx_counts(uint32_t *ctx_counts,
										 const unify_sketch_t *sketch,
										 int offset_gid,
										 int block_size)
{
	for (int i = 0; i < block_size; ++i) {
		const int gid = offset_gid + i;
		const uint64_t begin = sketch->sketch_index[gid];
		const uint64_t end = sketch->sketch_index[gid + 1];
		const uint64_t *arr = sketch->comb_sketch + begin;
		const size_t len = (size_t)(end - begin);
		ctx_counts[gid] = sketch->conflict
							   ? pairwise_count_ctx_runs_sorted_ctxobj64(arr, len)
							   : (uint32_t)len;
	}
}

static bool matrix_write_keep_matrix_indexed(const matrix_opt_t *opt,
											 const unify_sketch_t *sketch,
											 const int *kept,
											 int kept_n,
											 const char *sketch_dir)
{
	if (!matrix_metric_can_use_context_index(opt) || sketch->stat_type != 2)
		return false;

	const int n = sketch->infile_num;
	int *rank_by_gid = n > 0 ? calloc((size_t)n, sizeof(rank_by_gid[0])) : NULL;
	bool *keep_flag = n > 0 ? calloc((size_t)n, sizeof(keep_flag[0])) : NULL;
	if ((n > 0 && !rank_by_gid) || (n > 0 && !keep_flag))
		err(EXIT_FAILURE, "%s(): OOM keep rank", __func__);
	for (int i = 0; i < kept_n; ++i) {
		rank_by_gid[kept[i]] = i + 1;
		keep_flag[kept[i]] = true;
	}

	size_t index_bytes = 0;
	bool index_from_file = false;
	ctxgidobj_t *sorted_index =
		matrix_load_sorted_index(sketch, sketch_dir, &index_bytes, &index_from_file);
	const size_t ref_sksize = (size_t)sketch->sketch_index[n];

	uint32_t *ctx_counts = matrix_context_counts_for_sketch(sketch, false);
	const int block_size = matrix_choose_context_block_size(n, n);
	if (block_size <= 0)
		errx(EINVAL, "%s(): invalid matrix block size", __func__);
	ctx_mut2_t *ctx = malloc((size_t)n * (size_t)block_size * sizeof(ctx[0]));
	obj_section_t *obj = malloc((size_t)n * (size_t)block_size * sizeof(obj[0]));
	if (!ctx || !obj)
		err(EXIT_FAILURE, "%s(): OOM indexed kept matrix block", __func__);

	ani_opt_t ani_opt = {0};
	ani_opt.p = opt->p > 0 ? opt->p : 1;
	ani_opt.fmt = 1;
	ani_opt.ntop = -1;
	ani_opt.e = (int)opt->e;
	ani_opt.s = 4;

	FILE *output = matrix_open_output(opt->keep_matrix_outf);
	if (opt->keep_matrix_format == MATRIX_KEEP_MATRIX_PHYLIP) {
		fprintf(output, "%d\n", kept_n);
	} else {
		for (int q = 0; q < kept_n; ++q)
			fprintf(output, "\t%s", sketch->gname[kept[q]]);
		fprintf(output, "\n");
	}

	int rows_done = 0;
	kssd_progress_t progress =
		matrix_progress_start(opt, "indexed kept-matrix rows", "rows",
							  (uint64_t)kept_n);
	for (int offset = 0; offset < n; offset += block_size) {
		const int this_block = (offset + block_size > n) ? n - offset : block_size;
		matrix_fill_query_ctx_counts(ctx_counts, sketch, offset, this_block);
		memset(ctx, 0, (size_t)n * (size_t)this_block * sizeof(ctx[0]));
		memset(obj, 0, (size_t)n * (size_t)this_block * sizeof(obj[0]));
		count_ctx_obj_frm_comb_sketch_section(
			ctx, obj, sorted_index, ref_sksize, n, ctx_counts,
			ctx_counts + offset, this_block,
			sketch->comb_sketch + sketch->sketch_index[offset],
			sketch->sketch_index + offset, NULL, NULL, &ani_opt);

		for (int i = 0; i < this_block; ++i) {
			const int qgid = offset + i;
			if (!keep_flag[qgid])
				continue;
			if (opt->keep_matrix_format == MATRIX_KEEP_MATRIX_PHYLIP) {
				char id[16];
				matrix_make_phylip_id(rank_by_gid[qgid], id, sizeof(id));
				fprintf(output, "%-10s", id);
			} else {
				fprintf(output, "%s", sketch->gname[qgid]);
			}
			for (int k = 0; k < kept_n; ++k) {
				const int rgid = kept[k];
				double value = opt->diagonal_value;
				if (qgid != rgid) {
					ani_features_t features = {
						.XnY_ctx = MCTX(n, i, rgid).num_ctx,
						.N_diff_obj = MOBJ(n, i, rgid).diff_obj,
						.N_diff_obj_section = MOBJ(n, i, rgid).diff_obj_section,
						.N_mut2_ctx = MCTX(n, i, rgid).num_mut2_ctx,
					};
					value = matrix_context_feature_value(
						opt, &features, ctx_counts[qgid], ctx_counts[rgid]);
				}
				if (opt->keep_matrix_format == MATRIX_KEEP_MATRIX_PHYLIP)
					fprintf(output, " %.12g", value);
				else
					fprintf(output, "\t%lf", value);
			}
			fprintf(output, "\n");
			++rows_done;
			kssd_progress_update(&progress, (uint64_t)rows_done, false);
		}
	}
	kssd_progress_done(&progress);

	if (output != stdout)
		fclose(output);
	free(ctx);
	free(obj);
	free(ctx_counts);
	matrix_free_sorted_index(sorted_index, index_bytes, index_from_file);
	free(rank_by_gid);
	free(keep_flag);
	return true;
}

static void matrix_write_keep_matrix(const matrix_opt_t *opt,
									 const unify_sketch_t *sketch,
									 const pairwise_component_result_t *comp,
									 const char *sketch_dir)
{
	if (opt->keep_matrix_outf[0] == '\0')
		return;

	int kept_n = 0;
	int *kept = matrix_collect_kept_indices(comp, &kept_n);
	matrix_write_keep_matrix_idmap(opt, sketch, kept, kept_n);
	if (kept_n > 0 &&
		matrix_write_keep_matrix_indexed(opt, sketch, kept, kept_n, sketch_dir)) {
		free(kept);
		return;
	}
	FILE *output = matrix_open_output(opt->keep_matrix_outf);
	if (opt->keep_matrix_format == MATRIX_KEEP_MATRIX_PHYLIP) {
		fprintf(output, "%d\n", kept_n);
	} else {
		for (int q = 0; q < kept_n; ++q)
			fprintf(output, "\t%s", sketch->gname[kept[q]]);
		fprintf(output, "\n");
	}

	double *row = kept_n > 0 ? malloc((size_t)kept_n * sizeof(row[0])) : NULL;
	if (kept_n > 0 && !row)
		err(EXIT_FAILURE, "%s(): OOM kept matrix row", __func__);
	kssd_progress_t progress =
		matrix_progress_start(opt, "kept-matrix rows", "rows",
							  (uint64_t)kept_n);
	for (int r = 0; r < kept_n; ++r) {
		const int rn = kept[r];
#pragma omp parallel for num_threads((opt->p)) schedule(dynamic, 1)
		for (int q = 0; q < kept_n; ++q) {
			const int qn = kept[q];
			if (rn == qn) {
				row[q] = opt->diagonal_value;
			} else {
				pairwise_eval_t eval =
					pairwise_eval_expr_samples(&opt->metric, sketch, (uint32_t)qn,
											   sketch, (uint32_t)rn, false);
				row[q] = matrix_eval_value(opt, &eval);
			}
		}
		if (opt->keep_matrix_format == MATRIX_KEEP_MATRIX_PHYLIP) {
			char id[16];
			matrix_make_phylip_id(r + 1, id, sizeof(id));
			fprintf(output, "%-10s", id);
			for (int q = 0; q < kept_n; ++q)
				fprintf(output, " %.12g", row[q]);
		} else {
			fprintf(output, "%s", sketch->gname[rn]);
			for (int q = 0; q < kept_n; ++q)
				fprintf(output, "\t%lf", row[q]);
		}
		fprintf(output, "\n");
		kssd_progress_update(&progress, (uint64_t)r + 1, false);
	}
	kssd_progress_done(&progress);
	free(row);
	if (output != stdout)
		fclose(output);
	free(kept);
}

static void matrix_compute_triangle_report(const matrix_opt_t *opt,
										   const unify_sketch_t *sketch)
{
	FILE *output = matrix_open_output(opt->outf);
	const uint64_t n = (uint64_t)sketch->infile_num;
	const uint64_t total_pairs = n * (n - 1) / 2;
	uint64_t completed_pairs = 0;
	kssd_progress_t progress =
		matrix_progress_start(opt, "triangle pairs", "pairs", total_pairs);
	for (int qn = 0; qn < sketch->infile_num; qn++) {
		fprintf(output, "%s", sketch->gname[qn]);
		for (int rn = 0; rn < qn; rn++) {
			pairwise_eval_t eval =
				pairwise_eval_expr_samples(&opt->metric, sketch, (uint32_t)qn,
										   sketch, (uint32_t)rn, false);
			fprintf(output, "\t%lf", matrix_eval_value(opt, &eval));
		}
		if (opt->d)
			fprintf(output, "\t%lf", opt->diagonal_value);
		fprintf(output, "\n");
		completed_pairs += (uint64_t)qn;
		kssd_progress_update(&progress, completed_pairs, false);
	}
	kssd_progress_done(&progress);
	if (output != stdout)
		fclose(output);
}

static void matrix_compute_edges(const matrix_opt_t *opt, const unify_sketch_t *ref,
								 const unify_sketch_t *qry, bool same_sketch)
{
	FILE *output = matrix_open_output(opt->outf);
	kssd_progress_t progress =
		matrix_progress_start(opt, "edge rows", "rows",
							  (uint64_t)ref->infile_num);
	matrix_print_edge_header(output);
	for (int rn = 0; rn < ref->infile_num; rn++) {
		const int q_start = same_sketch ? rn + 1 : 0;
		for (int qn = q_start; qn < qry->infile_num; qn++) {
			pairwise_eval_t eval =
				pairwise_eval_expr_samples(&opt->metric, qry, (uint32_t)qn,
										   ref, (uint32_t)rn, false);
			if (pairwise_eval_passes_edge(&eval, opt->cut, opt->ctxcut, opt->max_afcut))
				matrix_print_edge(output, qry->gname[qn], ref->gname[rn],
								  &eval, &opt->metric);
		}
		kssd_progress_update(&progress, (uint64_t)rn + 1, false);
	}
	kssd_progress_done(&progress);
	if (output != stdout)
		fclose(output);
}

typedef struct matrix_component_ctx
{
	const matrix_opt_t *opt;
	const unify_sketch_t *sketch;
	FILE *edge_out;
} matrix_component_ctx_t;

static pairwise_eval_t matrix_component_eval_pair(void *ctx, int a, int b)
{
	const matrix_component_ctx_t *mctx = ctx;
	return pairwise_eval_expr_samples(&mctx->opt->metric, mctx->sketch, (uint32_t)a,
									  mctx->sketch, (uint32_t)b, false);
}

static int matrix_component_quality_compare(void *ctx, int a, int b)
{
	const matrix_component_ctx_t *mctx = ctx;
	return pairwise_quality_compare(mctx->sketch, a, b);
}

static void matrix_component_edge_observer(void *ctx, int a, int b,
										   const pairwise_eval_t *eval)
{
	const matrix_component_ctx_t *mctx = ctx;
	matrix_print_edge(mctx->edge_out, mctx->sketch->gname[a],
					  mctx->sketch->gname[b], eval, &mctx->opt->metric);
}

static bool matrix_can_use_sorted_index(const matrix_opt_t *opt,
										const unify_sketch_t *sketch,
										const char *sketch_dir)
{
	(void)opt;
	return sketch->stat_type == 2 &&
		   !sketch->conflict &&
		   file_exists_in_folder((char *)sketch_dir, (char *)sorted_comb_ctxgid64obj32);
}

static void matrix_write_component_report(const matrix_opt_t *opt,
										  const unify_sketch_t *sketch,
										  bool dedup_plan,
										  const pairwise_component_result_t *comp,
										  const char *sketch_dir)
{
	const int n = sketch->infile_num;
	int *cluster_id = n > 0 ? calloc((size_t)n, sizeof(cluster_id[0])) : NULL;
	if (n > 0 && !cluster_id)
		err(EXIT_FAILURE, "%s(): OOM cluster ids", __func__);

	int next_cluster = 1;
	for (int i = 0; i < n; ++i) {
		const int root = pairwise_component_find(comp, i);
		if (cluster_id[root] == 0)
			cluster_id[root] = next_cluster++;
	}

	FILE *output = matrix_open_output(opt->outf);
	FILE *keep_out = NULL;
	FILE *remove_out = NULL;
	if (dedup_plan && opt->keep_outf[0] != '\0')
		keep_out = matrix_open_output(opt->keep_outf);
	if (dedup_plan && opt->remove_outf[0] != '\0')
		remove_out = matrix_open_output(opt->remove_outf);
	if (dedup_plan)
		fprintf(output, "sample\tcomponent_id\taction\trepresentative\treason\trep_rule\tcomponent_size\tsample_sketch_entries\trep_sketch_entries\tsample_fmt\tsample_asm_level\tsample_total_length_bp\tsample_record_count\tsample_median_length_bp\tsample_length_cv\trep_fmt\trep_asm_level\trep_total_length_bp\trep_record_count\trep_median_length_bp\trep_length_cv\n");
	else
		fprintf(output, "sample\tcluster_id\tcomponent_size\trepresentative\n");

	for (int i = 0; i < n; ++i) {
		const int root = pairwise_component_find(comp, i);
		const int rep = comp->representative[root];
		if (dedup_plan) {
			const bool singleton = comp->component_size[root] == 1;
			const bool keep = i == rep;
			const char *reason = singleton ? "singleton" : (keep ? "best_representative" : "duplicate");
			if (keep && keep_out)
				fprintf(keep_out, "%s\n", sketch->gname[i]);
			if (!keep && remove_out)
				fprintf(remove_out, "%s\n", sketch->gname[i]);
			fprintf(output, "%s\tcluster_%06d\t%s\t%s\t%s\t%s\t%d\t%" PRIu64 "\t%" PRIu64,
					sketch->gname[i], cluster_id[root],
					keep ? "keep" : "remove", sketch->gname[rep], reason,
					matrix_rep_rule(sketch, i, rep), comp->component_size[root],
					matrix_sketch_entries(sketch, i),
					matrix_sketch_entries(sketch, rep));
			matrix_print_meta_field(output, sketch, i);
			matrix_print_meta_field(output, sketch, rep);
			fprintf(output, "\n");
		} else {
			fprintf(output, "%s\tcluster_%06d\t%d\t%s\n",
					sketch->gname[i], cluster_id[root],
					comp->component_size[root], sketch->gname[rep]);
		}
	}
	if (output != stdout)
		fclose(output);
	if (keep_out && keep_out != stdout)
		fclose(keep_out);
	if (remove_out && remove_out != stdout)
		fclose(remove_out);
	matrix_write_keep_matrix(opt, sketch, comp, sketch_dir);
	free(cluster_id);
}

typedef struct matrix_indexed_edge_ctx
{
	const matrix_opt_t *opt;
	const unify_sketch_t *sketch;
	pairwise_component_result_t *comp;
	pairwise_edge_list_t *edges;
	FILE *edge_out;
} matrix_indexed_edge_ctx_t;

static void matrix_indexed_edge_observer(void *ctx, int qry, int ref,
										 const pairwise_eval_t *eval)
{
	matrix_indexed_edge_ctx_t *mctx = ctx;
	if (mctx->comp)
		pairwise_component_result_union(mctx->comp, qry, ref);
	if (mctx->edges)
		pairwise_edge_list_add(mctx->edges, qry, ref);
	if (mctx->edge_out)
		matrix_print_edge(mctx->edge_out, mctx->sketch->gname[qry],
						  mctx->sketch->gname[ref], eval, &mctx->opt->metric);
}

static void matrix_compute_indexed_sparse(const matrix_opt_t *opt,
										  const char *sketch_dir,
										  const unify_sketch_t *sketch,
										  bool components,
										  bool dedup_plan)
{
	const int n = sketch->infile_num;
	const char *index_path = test_get_fullpath(sketch_dir, sorted_comb_ctxgid64obj32);

	FILE *edge_out = NULL;
	if (components) {
		if (opt->edge_outf[0] != '\0') {
			edge_out = matrix_open_output(opt->edge_outf);
			matrix_print_edge_header(edge_out);
		}
	} else {
		edge_out = matrix_open_output(opt->outf);
		matrix_print_edge_header(edge_out);
	}

	pairwise_component_result_t comp = {0};
	pairwise_edge_list_t dedup_edges = {0};
	if (components && !dedup_plan)
		pairwise_component_result_init(&comp, n);

	matrix_indexed_edge_ctx_t edge_ctx = {
		.opt = opt,
		.sketch = sketch,
		.comp = components && !dedup_plan ? &comp : NULL,
		.edges = components && dedup_plan ? &dedup_edges : NULL,
		.edge_out = edge_out,
	};
	pairwise_index_scan_options_t scan_opt = {
		.metric = &opt->metric,
		.cut = opt->cut,
		.ctxcut = opt->ctxcut,
		.max_afcut = opt->max_afcut,
		.index_max_ctx_freq = opt->index_max_ctx_freq,
		.index_min_votes = opt->index_min_votes,
		.index_sample_step = opt->index_sample_step,
		.threads = opt->p,
	};
	pairwise_index_scan_stats_t scan_stats = {0};
	kssd_progress_t progress =
		matrix_progress_start(opt,
							  dedup_plan ? "indexed dedup rows" :
							  components ? "indexed component rows" : "indexed edge rows",
							  "rows", (uint64_t)n);
	pairwise_indexed_self_scan(sketch, index_path, &scan_opt,
							   matrix_indexed_edge_observer, &edge_ctx,
							   kssd_progress_cb, &progress,
							   &scan_stats);
	kssd_progress_done(&progress);

	if (edge_out && edge_out != stdout)
		fclose(edge_out);

	if (components) {
		comp.distance_edges = scan_stats.distance_edges;
		comp.duplicate_edges = scan_stats.accepted_edges;
		comp.ctx_rejects = scan_stats.ctx_rejects;
		comp.max_af_rejects = scan_stats.max_af_rejects;
		matrix_component_ctx_t comp_ctx = {
			.opt = opt,
			.sketch = sketch,
			.edge_out = NULL,
		};
		if (dedup_plan) {
			if (opt->dedup_strategy == PAIRWISE_DEDUP_COMPLETE_LINKAGE) {
				pairwise_complete_linkage_dedup_from_edges(n, dedup_edges.edges,
														   dedup_edges.n,
														   &comp_ctx,
														   matrix_component_quality_compare,
														   &comp);
			} else {
				pairwise_greedy_dedup_from_edges(n, dedup_edges.edges,
												 dedup_edges.n,
												 &comp_ctx,
												 matrix_component_quality_compare,
												 &comp);
			}
			comp.distance_edges = scan_stats.distance_edges;
			comp.duplicate_edges = scan_stats.accepted_edges;
			comp.ctx_rejects = scan_stats.ctx_rejects;
			comp.max_af_rejects = scan_stats.max_af_rejects;
		} else {
			pairwise_component_result_finalize(&comp, &comp_ctx,
											   matrix_component_quality_compare);
		}
		matrix_write_component_report(opt, sketch, dedup_plan, &comp, sketch_dir);
		pairwise_component_result_free(&comp);
		pairwise_edge_list_free(&dedup_edges);
	}
}

static void matrix_compute_components(const matrix_opt_t *opt,
									  const unify_sketch_t *sketch,
									  bool dedup_plan)
{
	const int n = sketch->infile_num;
	FILE *edge_out = NULL;
	if (opt->edge_outf[0] != '\0') {
		edge_out = matrix_open_output(opt->edge_outf);
		matrix_print_edge_header(edge_out);
	}

	matrix_component_ctx_t ctx = {
		.opt = opt,
		.sketch = sketch,
		.edge_out = edge_out,
	};
	pairwise_component_result_t comp = {0};
	const uint64_t total_pairs = (uint64_t)n * (uint64_t)(n - 1) / 2;
	kssd_progress_t progress =
		matrix_progress_start(opt,
							  dedup_plan ? "dedup pairs" : "cluster pairs",
							  "pairs", total_pairs);
	if (dedup_plan) {
		if (opt->dedup_strategy == PAIRWISE_DEDUP_COMPLETE_LINKAGE) {
			pairwise_build_complete_linkage_dedup(n, opt->cut, opt->ctxcut,
												  opt->max_afcut, &ctx,
												  matrix_component_eval_pair,
												  matrix_component_quality_compare,
												  edge_out ? matrix_component_edge_observer : NULL,
												  &ctx, kssd_progress_cb, &progress,
												  &comp);
		} else {
			pairwise_build_greedy_dedup(n, opt->cut, opt->ctxcut, opt->max_afcut,
										&ctx, matrix_component_eval_pair,
										matrix_component_quality_compare,
										edge_out ? matrix_component_edge_observer : NULL,
										&ctx, kssd_progress_cb, &progress,
										&comp);
		}
	} else {
		pairwise_build_components(n, opt->cut, opt->ctxcut, opt->max_afcut,
								  &ctx, matrix_component_eval_pair,
								  matrix_component_quality_compare,
								  edge_out ? matrix_component_edge_observer : NULL,
								  &ctx, kssd_progress_cb, &progress,
								  &comp);
	}
	kssd_progress_done(&progress);
	if (edge_out && edge_out != stdout)
		fclose(edge_out);

	matrix_write_component_report(opt, sketch, dedup_plan, &comp, NULL);
	pairwise_component_result_free(&comp);
}

int compute_triangle(matrix_opt_t *matrix_opt)
{
	const matrix_format_t fmt = matrix_effective_format(matrix_opt, false);
	unify_sketch_t *result = generic_sketch_parse(
		matrix_opt->qrydir, matrix_parse_flags_for_format(fmt));
	pairwise_prepare_lco_model(result);
	if (fmt == MATRIX_FORMAT_FULL) {
		if (!matrix_write_indexed_self_matrix(matrix_opt, result, matrix_opt->qrydir, fmt))
			matrix_compute_dense(matrix_opt, result, result, true);
	}
	else if (fmt == MATRIX_FORMAT_EDGES) {
		if (matrix_can_use_sorted_index(matrix_opt, result, matrix_opt->qrydir))
			matrix_compute_indexed_sparse(matrix_opt, matrix_opt->qrydir, result, false, false);
		else
			matrix_compute_edges(matrix_opt, result, result, true);
	} else if (fmt == MATRIX_FORMAT_CLUSTERS) {
		if (matrix_can_use_sorted_index(matrix_opt, result, matrix_opt->qrydir))
			matrix_compute_indexed_sparse(matrix_opt, matrix_opt->qrydir, result, true, false);
		else
			matrix_compute_components(matrix_opt, result, false);
	}
	else if (fmt == MATRIX_FORMAT_DEDUP_PLAN) {
		if (!matrix_opt->max_afcut_set)
			matrix_opt->max_afcut = 0.8;
		if (matrix_can_use_sorted_index(matrix_opt, result, matrix_opt->qrydir))
			matrix_compute_indexed_sparse(matrix_opt, matrix_opt->qrydir, result, true, true);
		else
			matrix_compute_components(matrix_opt, result, true);
	} else
	{
		if (!matrix_write_indexed_self_matrix(matrix_opt, result, matrix_opt->qrydir, fmt))
			matrix_compute_triangle_report(matrix_opt, result);
	}
	matrix_write_glist(matrix_opt, result);
	free_unify_sketch(result);
	return 1;
}

int compute_matrix(matrix_opt_t *matrix_opt)
{
	unify_sketch_t *ref_result = generic_sketch_parse(matrix_opt->refdir, SKETCH_PARSE_NONE);
	unify_sketch_t *qry_result = generic_sketch_parse(matrix_opt->qrydir, SKETCH_PARSE_NONE);
	pairwise_check_compatible(ref_result, qry_result);
	pairwise_prepare_lco_model(ref_result);
	const bool same_sketch = strcmp(matrix_opt->refdir, matrix_opt->qrydir) == 0;
	const matrix_format_t fmt = matrix_effective_format(matrix_opt, true);
	if (fmt == MATRIX_FORMAT_EDGES)
		matrix_compute_edges(matrix_opt, ref_result, qry_result, same_sketch);
	else if (fmt == MATRIX_FORMAT_FULL)
		matrix_compute_dense(matrix_opt, ref_result, qry_result, same_sketch);
	else
		errx(EINVAL, "%s(): --format clusters, dedup-plan, and triangle require one sketch", __func__);
	free_unify_sketch(ref_result);
	free_unify_sketch(qry_result);
	return 1;
}

#define CONFLICT_OBJ (0)
KHASH_MAP_INIT_INT64(u64, uint64_t)
int compute_ani_matrix(matrix_opt_t *matrix_opt)
{ // ref is the sketch(es) to be hashed.

	unify_sketch_t *ref_result = generic_sketch_parse(matrix_opt->refdir, SKETCH_PARSE_NONE), *qry_result = generic_sketch_parse(matrix_opt->qrydir, SKETCH_PARSE_NONE);
	if (ref_result->stat_type != qry_result->stat_type)
		err(EXIT_FAILURE, "%s(): ref sketch type %u != qry %u", __func__, ref_result->stat_type, qry_result->stat_type);
	else if (ref_result->hash_id != qry_result->hash_id)
		err(EXIT_FAILURE, "%s(): ref hash_id %u != qry %u", __func__, ref_result->hash_id, qry_result->hash_id);
	dim_sketch_stat_t *lco_stat_readin = (dim_sketch_stat_t *)ref_result->mem_stat;
	int obj_len = lco_stat_readin->klen - 2 * lco_stat_readin->hclen;
	if (obj_len == 0)
		err(EXIT_FAILURE, "%s() abort!: sketching mode has 0bp object", __func__);

	uint64_t tmp_var = UINT64_MAX >> (64 - 2 * lco_stat_readin->hclen);
	uint64_t ctxmask = (tmp_var << (2 * (lco_stat_readin->klen - lco_stat_readin->hclen - lco_stat_readin->holen))) | (tmp_var << (2 * (lco_stat_readin->holen)));
	//  int kmerlen = ref_result->kmerlen;
	// print header
	FILE *output = matrix_opt->outf[0] == '\0' ? stdout : fopen(matrix_opt->outf, "w");
	if (output == NULL)
		err(errno, "%s(): %s", __func__, matrix_opt->outf);
	kssd_progress_t progress =
		matrix_progress_start(matrix_opt, "legacy ANI rows", "rows",
							  (uint64_t)ref_result->infile_num);
	for (int qn = 0; qn < qry_result->infile_num; qn++)
		fprintf(output, "\t%s", qry_result->gname[qn]);
	fprintf(output, "\n");

	co_distance_t *ctx_diff_obj_cnt = malloc(qry_result->infile_num * sizeof(co_distance_t));
	for (int rn = 0; rn < ref_result->infile_num; rn++)
	{
		//     int X_size = ref_result->sketch_index[rn+1] - ref_result->sketch_index[rn];
		khash_t(u64) *h = kh_init(u64);
		int ret;
		// hash rn-th ref genome
		for (uint64_t ri = ref_result->sketch_index[rn]; ri < ref_result->sketch_index[rn + 1]; ri++)
		{
			khiter_t k = kh_put(u64, h, ref_result->comb_sketch[ri] & ctxmask, &ret);
			if (ret == 1)
				kh_value(h, k) = ref_result->comb_sketch[ri];
			else if (ret == 0)
				kh_value(h, k) = CONFLICT_OBJ; // ret == 0: mark confict obj using 0; assuming no duplicted k-mers in a sketch
		}
		memset(ctx_diff_obj_cnt, 0, qry_result->infile_num * sizeof(co_distance_t));
#pragma omp parallel for num_threads((matrix_opt->p))
		for (int qn = 0; qn < qry_result->infile_num; qn++)
		{
			for (uint64_t qi = qry_result->sketch_index[qn]; qi < qry_result->sketch_index[qn + 1]; qi++)
			{
				khiter_t it = kh_get(u64, h, qry_result->comb_sketch[qi] & ctxmask);
				if (it != kh_end(h) && kh_value(h, it) != CONFLICT_OBJ)
				{
					ctx_diff_obj_cnt[qn].ctx_ct++;
					if (kh_value(h, it) != qry_result->comb_sketch[qi])
						ctx_diff_obj_cnt[qn].diff_obj++;
				}
			} // for qi
		} // for qn
		fprintf(output, "%s", ref_result->gname[rn]);
		for (int qn = 0; qn < qry_result->infile_num; qn++)
		{
			// int Y_size =  qry_result->sketch_index[qn+1] - qry_result->sketch_index[qn];
			double ani = 0;
			if (ctx_diff_obj_cnt[qn].ctx_ct > 100)
			{
				double dist = (double)ctx_diff_obj_cnt[qn].diff_obj / ctx_diff_obj_cnt[qn].ctx_ct;
				ani = pow((1 - dist), (1.0 / obj_len));
			}
			//	double dist = ctx_diff_obj_cnt[qn].ctx_ct == 0 ? matrix_opt->e : (double) ctx_diff_obj_cnt[qn].diff_obj / ctx_diff_obj_cnt[qn].ctx_ct ;
			fprintf(output, "\t%lf", ani);
		}
		fprintf(output, "\n");
		kh_destroy(u64, h);
		kssd_progress_update(&progress, (uint64_t)rn + 1, false);
	} // loop rn end
	kssd_progress_done(&progress);
	free_unify_sketch(ref_result);
	free_unify_sketch(qry_result);
	free(ctx_diff_obj_cnt);
	if (output != stdout)
		fclose(output);
	return 1;
}
