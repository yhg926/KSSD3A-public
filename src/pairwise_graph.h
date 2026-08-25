#ifndef PAIRWISE_GRAPH_H
#define PAIRWISE_GRAPH_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "global_basic.h"

typedef enum pairwise_metric
{
	PAIRWISE_METRIC_CTX_MOE = 0,
	PAIRWISE_METRIC_CTX_NAIVE = 1,
	PAIRWISE_METRIC_P_DIST = 2,
	PAIRWISE_METRIC_MASH = 3,
	PAIRWISE_METRIC_AAF = 4
} pairwise_metric_t;

typedef enum pairwise_dedup_strategy
{
	PAIRWISE_DEDUP_GREEDY = 0,
	PAIRWISE_DEDUP_COMPLETE_LINKAGE = 1
} pairwise_dedup_strategy_t;

#define PAIRWISE_METRIC_EXPR_MAX 4
#define PAIRWISE_METRIC_EXPR_LABEL_LEN 96

typedef enum pairwise_metric_expr_op
{
	PAIRWISE_METRIC_EXPR_SINGLE = 0,
	PAIRWISE_METRIC_EXPR_AND,
	PAIRWISE_METRIC_EXPR_OR
} pairwise_metric_expr_op_t;

typedef struct pairwise_metric_expr
{
	pairwise_metric_t metrics[PAIRWISE_METRIC_EXPR_MAX];
	size_t count;
	pairwise_metric_expr_op_t op;
	char label[PAIRWISE_METRIC_EXPR_LABEL_LEN];
} pairwise_metric_expr_t;

typedef struct pairwise_eval
{
	double distance;
	double similarity;
	uint32_t xny_ctx;
	uint32_t n_diff_obj;
	uint32_t n_diff_obj_section;
	uint32_t n_mut2_ctx;
	double af_qry;
	double af_ref;
	double max_af;
	bool valid;
} pairwise_eval_t;

typedef pairwise_eval_t (*pairwise_eval_pair_fn)(void *ctx, int a, int b);
typedef int (*pairwise_quality_compare_fn)(void *ctx, int a, int b);
typedef void (*pairwise_edge_observer_fn)(void *ctx, int a, int b,
										  const pairwise_eval_t *eval);
typedef void (*pairwise_index_edge_fn)(void *ctx, int qry, int ref,
									   const pairwise_eval_t *eval);
typedef void (*pairwise_progress_fn)(void *ctx, uint64_t done, uint64_t total);

typedef struct pairwise_index_scan_options
{
	const pairwise_metric_expr_t *metric;
	double cut;
	uint32_t ctxcut;
	double max_afcut;
	uint32_t index_max_ctx_freq;
	uint32_t index_min_votes;
	uint32_t index_sample_step;
	int threads;
} pairwise_index_scan_options_t;

typedef struct pairwise_index_scan_stats
{
	uint64_t candidate_pairs;
	uint64_t exact_pairs;
	uint64_t distance_edges;
	uint64_t accepted_edges;
	uint64_t ctx_rejects;
	uint64_t max_af_rejects;
} pairwise_index_scan_stats_t;

typedef struct pairwise_component_result
{
	int n;
	int *parent;
	int *component_size;
	int *representative;
	bool *remove_sample;
	uint8_t *rank;
	uint64_t distance_edges;
	uint64_t duplicate_edges;
	uint64_t ctx_rejects;
	uint64_t max_af_rejects;
	int duplicate_clusters;
	int removed_samples;
} pairwise_component_result_t;

typedef struct pairwise_edge_pair
{
	int a;
	int b;
} pairwise_edge_pair_t;

typedef struct pairwise_edge_list
{
	pairwise_edge_pair_t *edges;
	size_t n;
	size_t cap;
} pairwise_edge_list_t;

const char *pairwise_metric_name(pairwise_metric_t metric);
bool pairwise_metric_from_string(const char *arg, pairwise_metric_t *metric_out);
bool pairwise_metric_is_context(pairwise_metric_t metric);
const char *pairwise_dedup_strategy_name(pairwise_dedup_strategy_t strategy);
bool pairwise_dedup_strategy_from_string(const char *arg,
										 pairwise_dedup_strategy_t *strategy_out);
pairwise_metric_expr_t pairwise_metric_expr_single(pairwise_metric_t metric);
bool pairwise_metric_expr_from_string(const char *arg,
									  pairwise_metric_expr_t *expr_out);
bool pairwise_metric_expr_is_combined(const pairwise_metric_expr_t *expr);
bool pairwise_metric_expr_uses_context(const pairwise_metric_expr_t *expr);
pairwise_metric_t pairwise_metric_expr_primary(const pairwise_metric_expr_t *expr);
const char *pairwise_metric_expr_name(const pairwise_metric_expr_t *expr);

void pairwise_check_compatible(const unify_sketch_t *ref, const unify_sketch_t *qry);
void pairwise_prepare_lco_model(const unify_sketch_t *sketch);

uint32_t pairwise_count_ctx_runs_sorted_ctxobj64(const uint64_t *a, size_t n);

pairwise_eval_t pairwise_eval_arrays(pairwise_metric_t metric,
									 const dim_sketch_stat_t *stat,
									 const uint64_t *qry, size_t qry_n,
									 const uint64_t *ref, size_t ref_n,
									 uint32_t qry_ctx_count,
									 uint32_t ref_ctx_count,
									 bool ignore_ref_conflicts);

pairwise_eval_t pairwise_eval_samples(pairwise_metric_t metric,
									  const unify_sketch_t *qry, uint32_t qn,
									  const unify_sketch_t *ref, uint32_t rn,
									  bool ignore_ref_conflicts);

pairwise_eval_t pairwise_eval_expr_arrays(const pairwise_metric_expr_t *expr,
										  const dim_sketch_stat_t *stat,
										  const uint64_t *qry, size_t qry_n,
										  const uint64_t *ref, size_t ref_n,
										  uint32_t qry_ctx_count,
										  uint32_t ref_ctx_count,
										  bool ignore_ref_conflicts);

pairwise_eval_t pairwise_eval_expr_samples(const pairwise_metric_expr_t *expr,
										   const unify_sketch_t *qry, uint32_t qn,
										   const unify_sketch_t *ref, uint32_t rn,
										   bool ignore_ref_conflicts);

bool pairwise_eval_passes_edge(const pairwise_eval_t *eval, double cut,
							   uint32_t ctxcut, double max_afcut);

int pairwise_quality_compare(const unify_sketch_t *sketch, int a, int b);

void pairwise_component_result_free(pairwise_component_result_t *result);
int pairwise_component_find(const pairwise_component_result_t *result, int x);
void pairwise_component_result_init(pairwise_component_result_t *result, int n);
void pairwise_component_result_union(pairwise_component_result_t *result, int a, int b);
void pairwise_component_result_finalize(pairwise_component_result_t *result,
										void *quality_ctx,
										pairwise_quality_compare_fn quality_compare);
void pairwise_edge_list_free(pairwise_edge_list_t *list);
void pairwise_edge_list_add(pairwise_edge_list_t *list, int a, int b);
void pairwise_greedy_dedup_from_edges(int n,
									  const pairwise_edge_pair_t *edges,
									  size_t edge_n,
									  void *quality_ctx,
									  pairwise_quality_compare_fn quality_compare,
									  pairwise_component_result_t *result);
void pairwise_complete_linkage_dedup_from_edges(int n,
												const pairwise_edge_pair_t *edges,
												size_t edge_n,
												void *quality_ctx,
												pairwise_quality_compare_fn quality_compare,
												pairwise_component_result_t *result);
void pairwise_build_components(int n, double cut, uint32_t ctxcut,
							   double max_afcut, void *cb_ctx,
							   pairwise_eval_pair_fn eval_pair,
							   pairwise_quality_compare_fn quality_compare,
							   pairwise_edge_observer_fn edge_observer,
							   void *edge_ctx,
							   pairwise_progress_fn progress_fn,
							   void *progress_ctx,
							   pairwise_component_result_t *result);
void pairwise_build_greedy_dedup(int n, double cut, uint32_t ctxcut,
								 double max_afcut, void *cb_ctx,
								 pairwise_eval_pair_fn eval_pair,
								 pairwise_quality_compare_fn quality_compare,
								 pairwise_edge_observer_fn edge_observer,
								 void *edge_ctx,
								 pairwise_progress_fn progress_fn,
								 void *progress_ctx,
								 pairwise_component_result_t *result);
void pairwise_build_complete_linkage_dedup(int n, double cut, uint32_t ctxcut,
										   double max_afcut, void *cb_ctx,
										   pairwise_eval_pair_fn eval_pair,
										   pairwise_quality_compare_fn quality_compare,
										   pairwise_edge_observer_fn edge_observer,
										   void *edge_ctx,
										   pairwise_progress_fn progress_fn,
										   void *progress_ctx,
										   pairwise_component_result_t *result);
void pairwise_indexed_self_scan(const unify_sketch_t *sketch,
								const char *sorted_index_path,
								const pairwise_index_scan_options_t *opt,
								pairwise_index_edge_fn edge_fn,
								void *edge_ctx,
								pairwise_progress_fn progress_fn,
								void *progress_ctx,
								pairwise_index_scan_stats_t *stats);

#endif
