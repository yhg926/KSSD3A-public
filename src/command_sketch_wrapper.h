#ifndef SKETCH_WRAPPER_H
#define SKETCH_WRAPPER_H

#include "global_basic.h"
#include "pairwise_graph.h"
#include <stdbool.h>
#include <argp.h>

typedef pairwise_metric_expr_t sketch_dedup_metric_t;

typedef struct sketch_opt
{
	// sketch stats
	int hclen;	// half context length, 1..16
	int holen;	// half outer object length, 0..8
	int iolen;	// half outer object length, 0..8
	int drfold; // dimension reduction fold 2^n , 0..32
	int kmerocrs;
	double npercentile; // dynamic k-mer occurrence threshold percentile, 0 disables
	int ncap; // maximum effective lower k-mer occurrence threshold, 0 disables
	bool reads_qc; // infer a read k-mer count range from the frequency distribution
	int p; // threads counts
	bool abundance;
	bool asone;	// treat input genomes as parts of final genome.
	bool conflict; // keep conflict context-object or not 
	bool anno; // write FASTA/FASTQ header annotations to lcofiles.anno
	bool compute_meta; // write per-input metadata to lcofiles.infilemeta
	bool position; // write per-context-object sequence positions to comblco.position
	bool drop_position; // drop comblco.position during filtering maintenance modes
	bool merge_comblco;
	bool append_comblco;
	bool remove_comblco;
	bool keep_comblco;
	bool dedup_comblco;
	bool dedup_raw_build_from_inputs;
	bool append_copy_mode;
	bool remove_copy_mode;
	bool keep_copy_mode;
	bool dedup_copy_mode;
	bool sketch_qc;
	double dedup_cutoff;
	double dedup_max_afcut;
	uint32_t dedup_ctxcut;
	sketch_dedup_metric_t dedup_metric;
	pairwise_dedup_strategy_t dedup_strategy;
	bool dedup_index;
	uint32_t dedup_index_max_ctx_freq;
	uint32_t dedup_index_min_votes;
	uint32_t dedup_index_sample_step;
	int print_mode; // 1 samples, 2 sketch entries, 3 sorted index, 4 positions
	bool split_mfa;
	bool coden_ctxobj_pattern; 	 
	char index[PATHLEN];
	char unique_index[PATHLEN];
	char *remove_list;
	char *remove_source;
	char *keep_list;
	char *keep_source;
	char *dedup_source;
	char *fpath;
	char *outdir;  // results dir
	char *pipecmd; // pipe command
	int num_remaining_args;
	char **remaining_args;
} sketch_opt_t;

int cmd_sketch(struct argp_state *state);
/*
extern void compute_sketch(sketch_opt_t*, infile_tab_t*);
extern void gen_inverted_index4comblco(const char* sketchdir);
extern int merge_comblco (sketch_opt_t * sketch_opt_val);
extern uint32_t get_sketching_id(uint32_t hclen, uint32_t holen,uint32_t iolen,uint32_t drfold,uint32_t FILTER);
*/
#endif
