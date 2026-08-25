#ifndef COMMAND_MATRIX
#define COMMAND_MATRIX
#include <argp.h>
#include <math.h>
#include <tgmath.h>
#include "global_basic.h"
#include "pairwise_graph.h"

// distance macro
#define MIN(X,Y) ( (X) < (Y) ? (X): (Y) )
#define U(X,Y,XnY) ( (X) + (Y) - (XnY))  
#define JCD(X,Y,XnY) ( (double)(XnY)/(U(X,Y,XnY)))
#define CTM(X,Y,XnY) ((double)(XnY)/MIN(X,Y))

#define MASHD (K,X,Y,XnY) (-log( 2*JCD(X,Y,XnY) / (1 + JCD(X,Y,XnY) )) / (K) )
#define AAFD(K,X,Y,XnY) (-log(CTM(X,Y,XnY)) / (K))

typedef enum matrix_format
{
	MATRIX_FORMAT_AUTO = 0,
	MATRIX_FORMAT_FULL,
	MATRIX_FORMAT_TRIANGLE,
	MATRIX_FORMAT_EDGES,
	MATRIX_FORMAT_CLUSTERS,
	MATRIX_FORMAT_DEDUP_PLAN
} matrix_format_t;

typedef enum matrix_keep_matrix_format
{
	MATRIX_KEEP_MATRIX_TSV = 0,
	MATRIX_KEEP_MATRIX_PHYLIP
} matrix_keep_matrix_format_t;

typedef enum matrix_progress_mode
{
	MATRIX_PROGRESS_AUTO = 0,
	MATRIX_PROGRESS_ON,
	MATRIX_PROGRESS_OFF
} matrix_progress_mode_t;

typedef struct matrix_opt
{
	pairwise_metric_expr_t metric;
	matrix_format_t format;
	matrix_keep_matrix_format_t matrix_out_format;
	matrix_keep_matrix_format_t keep_matrix_format;
	matrix_progress_mode_t progress_mode;
	pairwise_dedup_strategy_t dedup_strategy;
	double c; // compatibility alias for cut
	double cut;
	bool cut_set;
	double max_afcut;
	bool max_afcut_set;
	uint32_t ctxcut;
	uint32_t index_max_ctx_freq;
	uint32_t index_min_votes;
	uint32_t index_sample_step;
	int p; //threads
	bool d; //diagnal
	double diagonal_value;
	bool ani;
	double e;
  char qrydir[PATHLEN];
	char refdir[PATHLEN];
  char outf[PATHLEN];
	char edge_outf[PATHLEN];
	char keep_outf[PATHLEN];
	char remove_outf[PATHLEN];
	char matrix_idmap_outf[PATHLEN];
	char keep_matrix_outf[PATHLEN];
	char keep_matrix_idmap_outf[PATHLEN];
	char gl[PATHLEN]; // genome list with selection code 	
  int num_remaining_args;
  char **remaining_args;

} matrix_opt_t;

double get_mashD (uint32_t K, uint32_t X, uint32_t Y, uint32_t XnY);

double get_aafD (uint32_t K, uint32_t X, uint32_t Y, uint32_t XnY);
typedef double (*Dist) (uint32_t,uint32_t,uint32_t,uint32_t);
int cmd_matrix(struct argp_state* state);
int compute_triangle(matrix_opt_t *);
int compute_matrix(matrix_opt_t *);
int compute_ani_matrix(matrix_opt_t *matrix_opt);
#endif
