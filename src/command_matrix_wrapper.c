#include "command_matrix.h"
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <argp.h>
#include <argz.h>
#include <err.h>
#include <errno.h>
#include <math.h>
#include <libgen.h>
#include <dirent.h>
#include <omp.h>
#include <stdatomic.h>
#include <ctype.h>
#include <limits.h>

/*** argp wrapper ***/
struct arg_matrix
{
  struct arg_global* global;

  char* name;
  int dedup_strategy_seen;
};

enum
{
  MATRIX_KEY_FORMAT = 1000,
  MATRIX_KEY_CUT,
  MATRIX_KEY_MAX_AFCUT,
  MATRIX_KEY_CTXCUT,
  MATRIX_KEY_EDGE_OUT,
  MATRIX_KEY_DIAGONAL_VALUE,
  MATRIX_KEY_KEEP_OUT,
  MATRIX_KEY_REMOVE_OUT,
  MATRIX_KEY_MATRIX_OUT_FORMAT,
  MATRIX_KEY_MATRIX_IDMAP,
  MATRIX_KEY_KEEP_MATRIX_OUT,
  MATRIX_KEY_KEEP_MATRIX_FORMAT,
  MATRIX_KEY_KEEP_MATRIX_IDMAP,
  MATRIX_KEY_INDEX_MAX_CTX_FREQ,
  MATRIX_KEY_INDEX_MIN_VOTES,
  MATRIX_KEY_INDEX_SAMPLE_STEP,
  MATRIX_KEY_PROGRESS,
  MATRIX_KEY_SPARSE,
  MATRIX_KEY_DEDUP_STRATEGY
};

static struct argp_option opt_matrix[] =
{
	{"ref",'r',"<DIR>", 0, "Reference sketches. Omit for triangle mode.",9},
	{"query",'q',"<DIR>", 0, "Query sketches.",1},
	{"metric",'m',"<METRIC|EXPR>", 0, "Pairwise metric. Sparse formats also accept quoted A&B or A|B expressions. Numeric 0/1 alias mash/aaf. [ctx-naive]",2},
	{"format", MATRIX_KEY_FORMAT, "<FORMAT>", 0, "Report format: full, triangle, edges/sparse, clusters, or dedup-plan. [full for -r/-q; triangle for one sketch]", 2},
	{"sparse", MATRIX_KEY_SPARSE, 0, 0, "Alias for --format sparse/edges; requires --cut.", 2},
	{"cut", MATRIX_KEY_CUT, "<FLOAT>", 0, "Distance cutoff for sparse edge, cluster, and dedup-plan reports.", 3},
	{"dedup-strategy", MATRIX_KEY_DEDUP_STRATEGY, "<greedy|full-linkage>", 0, "Dedup-plan grouping strategy: greedy representative-neighbor or full-linkage clique. [greedy]", 3},
	{"control",'c',"<FLOAT>", OPTION_HIDDEN, "Compatibility alias for --cut.",3},
	{"max-afcut", MATRIX_KEY_MAX_AFCUT, "<FLOAT>", 0, "Minimum max(Qry_align_fraction, Ref_align_fraction) for sparse edges. [0; 0.8 for dedup-plan]", 3},
	{"ctxcut", MATRIX_KEY_CTXCUT, "<INT>", 0, "Minimum shared context count for sparse edges. [0]", 3},
	{"index-max-ctx-freq", MATRIX_KEY_INDEX_MAX_CTX_FREQ, "<INT>", 0, "Indexed sparse candidate mode: skip context groups above this frequency; 0 keeps all contexts. [256]", 3},
	{"index-min-votes", MATRIX_KEY_INDEX_MIN_VOTES, "<INT>", 0, "Indexed sparse candidate mode: exact-score pairs with at least this many candidate context votes. [1]", 3},
	{"index-sample-step", MATRIX_KEY_INDEX_SAMPLE_STEP, "<INT>", 0, "Indexed sparse candidate mode: use every Nth query context for nomination. [1]", 3},
	{"progress", MATRIX_KEY_PROGRESS, "<auto|on|off>", 0, "Report matrix progress to stderr only when main output uses -o. [auto]", 3},
	{"glist",'g',"<FILE>",0,"Sample output file for KSSD set grouping.",4},
	{"outfile",'o',"<FILE>",0,"Matrix output file. [STDOUT]",5},
	{"matrix-format", MATRIX_KEY_MATRIX_OUT_FORMAT, "<tsv|phylip>", 0, "Output format for --format full. [tsv]", 5},
	{"matrix-idmap", MATRIX_KEY_MATRIX_IDMAP, "<FILE>", 0, "With --matrix-format phylip, write short-ID to sample-name map.", 5},
	{"edge-out", MATRIX_KEY_EDGE_OUT, "<FILE>", 0, "Optional sparse edge TSV for cluster/dedup-plan reports.", 5},
	{"keep-out", MATRIX_KEY_KEEP_OUT, "<FILE>", 0, "Write one kept sample per line for dedup-plan reports.", 5},
	{"remove-out", MATRIX_KEY_REMOVE_OUT, "<FILE>", 0, "Write one removed sample per line for dedup-plan reports.", 5},
	{"keep-matrix-out", MATRIX_KEY_KEEP_MATRIX_OUT, "<FILE>", 0, "With dedup-plan, write a dense distance matrix for kept representatives.", 5},
	{"keep-matrix-format", MATRIX_KEY_KEEP_MATRIX_FORMAT, "<tsv|phylip>", 0, "Output format for --keep-matrix-out. [tsv]", 5},
	{"keep-matrix-idmap", MATRIX_KEY_KEEP_MATRIX_IDMAP, "<FILE>", 0, "With --keep-matrix-format phylip, write short-ID to sample-name map.", 5},
	{"threads",'p',"<INT>", 0, "Number of threads to use.",6},
	{"diagonal",'d',0, 0, "Print diagonal values in triangle mode. [0]",7},
	{"diagonal-value", MATRIX_KEY_DIAGONAL_VALUE, "<FLOAT>", 0, "Set diagonal value and enable diagonal output.",7},
	{"exception",'e',"<FLOAT>", 0, "Distance value to use when XnY == 0. [1]",8},
	{"ani",'n',0, OPTION_HIDDEN, "Use context-object ANI mode.",9},
  { 0 }
};

static char doc_matrix[] =
  "\n"
  "Report pairwise sketch distance matrices, sparse edge lists, clusters, or dedup plans."
  "\v"
  "Use one sketch argument, -q, or -r/-q for a rectangular query-row by reference-column report.\n"
  "\n"
  "Examples:\n"
	  "  kssd3a matrix --format triangle -q sketches -o triangle.tsv\n"
	  "  kssd3a matrix -r ref_sketches -q qry_sketches --format full -o matrix.tsv\n"
	  "  kssd3a matrix --sparse --cut 0.05 sketches > sparse_edges.tsv\n"
	  "  kssd3a matrix --format edges -m 'ctx-naive&aaf' --cut 0.05 sketches\n"
	  "  kssd3a matrix --format clusters --cut 0.05 sketches > clusters.tsv\n"
  "  kssd3a matrix --format dedup-plan --cut 0.001 --keep-out keep.txt \\\n"
  "    --keep-matrix-out keep_matrix.tsv sketches\n"
  "  kssd3a matrix --format dedup-plan --cut 0.001 --keep-out keep.txt \\\n"
  "    --keep-matrix-out keep_matrix.phy --keep-matrix-format phylip \\\n"
  "    --keep-matrix-idmap keep_matrix.idmap.tsv sketches"
  ;

matrix_opt_t matrix_opt ={
	.metric = {
		.metrics = { PAIRWISE_METRIC_CTX_NAIVE },
		.count = 1,
		.op = PAIRWISE_METRIC_EXPR_SINGLE,
		.label = "ctx-naive",
	},
	.format = MATRIX_FORMAT_AUTO,
	.matrix_out_format = MATRIX_KEEP_MATRIX_TSV,
	.keep_matrix_format = MATRIX_KEEP_MATRIX_TSV,
	.progress_mode = MATRIX_PROGRESS_AUTO,
	.dedup_strategy = PAIRWISE_DEDUP_GREEDY,
	.c = 0.0,
	.cut = 0.0,
	.cut_set = false,
	.max_afcut = 0.0,
	.max_afcut_set = false,
	.ctxcut = 0,
	.index_max_ctx_freq = 256,
	.index_min_votes = 1,
	.index_sample_step = 1,
	.p = 1,
	.d = 0, //diagonal
	.diagonal_value = 0.0,
	.ani = 0,
	.e = 1.0,
	.refdir[0] = '\0',
	.qrydir[0] = '\0',
	.outf[0] = '\0',
	.edge_outf[0] = '\0',
	.keep_outf[0] = '\0',
	.remove_outf[0] = '\0',
	.matrix_idmap_outf[0] = '\0',
	.keep_matrix_outf[0] = '\0',
	.keep_matrix_idmap_outf[0] = '\0',
	.gl[0]= '\0',
	.num_remaining_args = 0, //int num_remaining_args; no option arguments num.
	.remaining_args = NULL //char **remaining_args; no option arguments array.
};

static int parse_int_range(struct argp_state *state, const char *option_name,
                           const char *arg, int min_value, int max_value)
{
  char *end = NULL;
  errno = 0;
  long value = strtol(arg, &end, 10);
  if (errno != 0 || end == arg || *end != '\0' || value < min_value || value > max_value)
    argp_error(state, "%s requires an integer in range %d..%d", option_name, min_value, max_value);
  return (int)value;
}

static double parse_nonnegative_double(struct argp_state *state, const char *option_name,
                                       const char *arg)
{
  char *end = NULL;
  errno = 0;
  double value = strtod(arg, &end);
  if (errno != 0 || end == arg || *end != '\0' || !isfinite(value) || value < 0.0)
    argp_error(state, "%s requires a non-negative number", option_name);
  return value;
}

static double parse_double_range(struct argp_state *state, const char *option_name,
                                 const char *arg, double min_value, double max_value)
{
  char *end = NULL;
  errno = 0;
  double value = strtod(arg, &end);
  if (errno != 0 || end == arg || *end != '\0' || !isfinite(value) ||
      value < min_value || value > max_value)
    argp_error(state, "%s requires a number in range %.12g..%.12g",
               option_name, min_value, max_value);
  return value;
}

static matrix_format_t parse_matrix_format(struct argp_state *state, const char *arg)
{
  if (strcmp(arg, "full") == 0 || strcmp(arg, "matrix") == 0)
    return MATRIX_FORMAT_FULL;
  if (strcmp(arg, "triangle") == 0 || strcmp(arg, "tri") == 0)
    return MATRIX_FORMAT_TRIANGLE;
  if (strcmp(arg, "edges") == 0 || strcmp(arg, "edge") == 0 || strcmp(arg, "sparse") == 0)
    return MATRIX_FORMAT_EDGES;
  if (strcmp(arg, "clusters") == 0 || strcmp(arg, "cluster") == 0 ||
      strcmp(arg, "components") == 0)
    return MATRIX_FORMAT_CLUSTERS;
  if (strcmp(arg, "dedup-plan") == 0 || strcmp(arg, "dedup") == 0)
    return MATRIX_FORMAT_DEDUP_PLAN;
  argp_error(state, "--format must be one of full, triangle, edges/sparse, clusters, or dedup-plan");
  return MATRIX_FORMAT_AUTO;
}

static matrix_keep_matrix_format_t parse_keep_matrix_format(struct argp_state *state,
                                                            const char *arg)
{
  if (strcmp(arg, "tsv") == 0 || strcmp(arg, "table") == 0)
    return MATRIX_KEEP_MATRIX_TSV;
  if (strcmp(arg, "phylip") == 0 || strcmp(arg, "phy") == 0 ||
      strcmp(arg, "fillip") == 0)
    return MATRIX_KEEP_MATRIX_PHYLIP;
  argp_error(state, "--keep-matrix-format must be one of tsv or phylip");
  return MATRIX_KEEP_MATRIX_TSV;
}

static matrix_progress_mode_t parse_matrix_progress(struct argp_state *state,
                                                    const char *arg)
{
  if (strcmp(arg, "auto") == 0)
    return MATRIX_PROGRESS_AUTO;
  if (strcmp(arg, "on") == 0 || strcmp(arg, "1") == 0 || strcmp(arg, "yes") == 0)
    return MATRIX_PROGRESS_ON;
  if (strcmp(arg, "off") == 0 || strcmp(arg, "0") == 0 || strcmp(arg, "no") == 0)
    return MATRIX_PROGRESS_OFF;
  argp_error(state, "--progress must be one of auto, on, or off");
  return MATRIX_PROGRESS_AUTO;
}

static void copy_path_arg(struct argp_state *state, const char *option_name,
                          char *dest, size_t dest_size, const char *arg)
{
  if (strlen(arg) >= dest_size)
    argp_error(state, "%s path is too long; maximum supported length is %zu bytes", option_name, dest_size - 1);
  snprintf(dest, dest_size, "%s", arg);
}

static error_t parse_matrix(int key, char* arg, struct argp_state* state) {
  struct arg_matrix* matrix = state->input;
	
  switch(key)
  {		
		case 'm':
		{
				if (!pairwise_metric_expr_from_string(arg, &matrix_opt.metric))
					argp_error(state, "-m/--metric must be one metric or an expression like 'ctx-naive&aaf' or 'ctx-moe|mash'");
				break;
		}
    case MATRIX_KEY_FORMAT:
    {
      matrix_opt.format = parse_matrix_format(state, arg);
      break;
    }
    case MATRIX_KEY_SPARSE:
    {
      matrix_opt.format = MATRIX_FORMAT_EDGES;
      break;
    }
		case 'c':
    case MATRIX_KEY_CUT:
		{
			matrix_opt.cut = parse_nonnegative_double(state, "--cut", arg);
			matrix_opt.c = matrix_opt.cut;
      matrix_opt.cut_set = true;
			break;
		}
    case MATRIX_KEY_MAX_AFCUT:
    {
      matrix_opt.max_afcut = parse_double_range(state, "--max-afcut", arg, 0.0, 1.0);
      matrix_opt.max_afcut_set = true;
      break;
    }
    case MATRIX_KEY_CTXCUT:
    {
      matrix_opt.ctxcut = (uint32_t)parse_int_range(state, "--ctxcut", arg, 0, INT32_MAX);
      break;
    }
    case MATRIX_KEY_INDEX_MAX_CTX_FREQ:
    {
      matrix_opt.index_max_ctx_freq = (uint32_t)parse_int_range(state, "--index-max-ctx-freq", arg, 0, INT32_MAX);
      break;
    }
    case MATRIX_KEY_INDEX_MIN_VOTES:
    {
      matrix_opt.index_min_votes = (uint32_t)parse_int_range(state, "--index-min-votes", arg, 1, INT32_MAX);
      break;
    }
    case MATRIX_KEY_INDEX_SAMPLE_STEP:
    {
      matrix_opt.index_sample_step = (uint32_t)parse_int_range(state, "--index-sample-step", arg, 1, INT32_MAX);
      break;
    }
    case MATRIX_KEY_PROGRESS:
    {
      matrix_opt.progress_mode = parse_matrix_progress(state, arg);
      break;
    }
    case MATRIX_KEY_DEDUP_STRATEGY:
    {
      matrix->dedup_strategy_seen = 1;
      if (!pairwise_dedup_strategy_from_string(arg, &matrix_opt.dedup_strategy))
        argp_error(state, "--dedup-strategy must be one of greedy, full-linkage, complete-linkage, or clique");
      break;
    }
		case 'e':
		{
      matrix_opt.e = parse_nonnegative_double(state, "-e/--exception", arg);
			break;
		}
		case 'p':
		{
			matrix_opt.p = parse_int_range(state, "-p/--threads", arg, 1, 65536);
			break;
		}
		case 'q':
		{
			copy_path_arg(state, "-q/--query", matrix_opt.qrydir, sizeof(matrix_opt.qrydir), arg);
			break;
		}
		case 'r':
    {
      copy_path_arg(state, "-r/--ref", matrix_opt.refdir, sizeof(matrix_opt.refdir), arg);
      break;
    }
		case 'o':
		{	
			copy_path_arg(state, "-o/--outfile", matrix_opt.outf, sizeof(matrix_opt.outf), arg);
			break;
		}
		case 'g':
		{
			copy_path_arg(state, "-g/--glist", matrix_opt.gl, sizeof(matrix_opt.gl), arg);
			break;
		}
    case MATRIX_KEY_EDGE_OUT:
    {
      copy_path_arg(state, "--edge-out", matrix_opt.edge_outf, sizeof(matrix_opt.edge_outf), arg);
      break;
    }
    case MATRIX_KEY_KEEP_OUT:
    {
      copy_path_arg(state, "--keep-out", matrix_opt.keep_outf, sizeof(matrix_opt.keep_outf), arg);
      break;
    }
    case MATRIX_KEY_REMOVE_OUT:
    {
      copy_path_arg(state, "--remove-out", matrix_opt.remove_outf, sizeof(matrix_opt.remove_outf), arg);
      break;
    }
    case MATRIX_KEY_MATRIX_OUT_FORMAT:
    {
      matrix_opt.matrix_out_format = parse_keep_matrix_format(state, arg);
      break;
    }
    case MATRIX_KEY_MATRIX_IDMAP:
    {
      copy_path_arg(state, "--matrix-idmap", matrix_opt.matrix_idmap_outf, sizeof(matrix_opt.matrix_idmap_outf), arg);
      break;
    }
    case MATRIX_KEY_KEEP_MATRIX_OUT:
    {
      copy_path_arg(state, "--keep-matrix-out", matrix_opt.keep_matrix_outf, sizeof(matrix_opt.keep_matrix_outf), arg);
      break;
    }
    case MATRIX_KEY_KEEP_MATRIX_FORMAT:
    {
      matrix_opt.keep_matrix_format = parse_keep_matrix_format(state, arg);
      break;
    }
    case MATRIX_KEY_KEEP_MATRIX_IDMAP:
    {
      copy_path_arg(state, "--keep-matrix-idmap", matrix_opt.keep_matrix_idmap_outf, sizeof(matrix_opt.keep_matrix_idmap_outf), arg);
      break;
    }
		case 'd':
		{
			matrix_opt.d = 1 ;
			break;
		}
    case MATRIX_KEY_DIAGONAL_VALUE:
    {
      matrix_opt.diagonal_value = parse_nonnegative_double(state, "--diagonal-value", arg);
      matrix_opt.d = 1;
      break;
    }
    case 'n':
    {
      matrix_opt.ani = 1 ;
      break;
    }
		case ARGP_KEY_ARGS:
		{
			matrix_opt.num_remaining_args = state->argc - state->next;
			matrix_opt.remaining_args  = state->argv + state->next;
			break;
		}
    case ARGP_KEY_END:
    {	
			if (matrix_opt.qrydir[0] == '\0' && matrix_opt.num_remaining_args == 1)
				copy_path_arg(state, "sketch", matrix_opt.qrydir, sizeof(matrix_opt.qrydir),
				              matrix_opt.remaining_args[0]);
			else if (matrix_opt.num_remaining_args > 0)
				argp_error(state, "matrix accepts at most one positional sketch argument; use -r/-q for rectangular reports");
			if(matrix_opt.qrydir[0] == '\0')
				argp_error(state, "missing sketch input; provide one sketch argument or -q/--query");
			if (matrix_opt.refdir[0] != '\0' && matrix_opt.format == MATRIX_FORMAT_TRIANGLE)
				argp_error(state, "--format triangle is only valid for one sketch");
			const bool rectangular = matrix_opt.refdir[0] != '\0';
			const matrix_format_t effective_format =
				matrix_opt.format != MATRIX_FORMAT_AUTO
					? matrix_opt.format
					: (rectangular ? MATRIX_FORMAT_FULL : MATRIX_FORMAT_TRIANGLE);
			if (pairwise_metric_expr_is_combined(&matrix_opt.metric) &&
			    (effective_format == MATRIX_FORMAT_FULL ||
			     effective_format == MATRIX_FORMAT_TRIANGLE))
				argp_error(state, "combined --metric expressions are only valid for --sparse/--format edges, clusters, or dedup-plan");
			if ((matrix_opt.format == MATRIX_FORMAT_EDGES ||
			     matrix_opt.format == MATRIX_FORMAT_CLUSTERS ||
			     matrix_opt.format == MATRIX_FORMAT_DEDUP_PLAN) &&
			    !matrix_opt.cut_set)
				argp_error(state, "--sparse/--format edges, clusters, or dedup-plan requires --cut");
			if ((matrix_opt.format == MATRIX_FORMAT_CLUSTERS ||
			     matrix_opt.format == MATRIX_FORMAT_DEDUP_PLAN) &&
			    matrix_opt.refdir[0] != '\0')
				argp_error(state, "--format clusters/dedup-plan requires one sketch, not -r/-q");
			if ((matrix_opt.keep_outf[0] != '\0' || matrix_opt.remove_outf[0] != '\0' ||
			     matrix_opt.keep_matrix_outf[0] != '\0' ||
			     matrix_opt.keep_matrix_idmap_outf[0] != '\0') &&
			    effective_format != MATRIX_FORMAT_DEDUP_PLAN)
				argp_error(state, "--keep-out/--remove-out/--keep-matrix-out/--keep-matrix-idmap require --format dedup-plan");
			if (matrix->dedup_strategy_seen && effective_format != MATRIX_FORMAT_DEDUP_PLAN)
				argp_error(state, "--dedup-strategy requires --format dedup-plan");
			if (matrix_opt.matrix_out_format != MATRIX_KEEP_MATRIX_TSV &&
			    effective_format != MATRIX_FORMAT_FULL)
				argp_error(state, "--matrix-format phylip requires --format full");
			if (matrix_opt.matrix_out_format == MATRIX_KEEP_MATRIX_PHYLIP &&
			    matrix_opt.refdir[0] != '\0')
				argp_error(state, "--matrix-format phylip requires one square sketch, not -r/-q");
			if (matrix_opt.matrix_idmap_outf[0] != '\0' &&
			    matrix_opt.matrix_out_format != MATRIX_KEEP_MATRIX_PHYLIP)
				argp_error(state, "--matrix-idmap requires --matrix-format phylip");
			if (matrix_opt.keep_matrix_format != MATRIX_KEEP_MATRIX_TSV &&
			    matrix_opt.keep_matrix_outf[0] == '\0')
				argp_error(state, "--keep-matrix-format requires --keep-matrix-out");
			if (matrix_opt.keep_matrix_idmap_outf[0] != '\0' &&
			    matrix_opt.keep_matrix_outf[0] == '\0')
				argp_error(state, "--keep-matrix-idmap requires --keep-matrix-out");
			if (matrix_opt.keep_matrix_idmap_outf[0] != '\0' &&
			    matrix_opt.keep_matrix_format != MATRIX_KEEP_MATRIX_PHYLIP)
				argp_error(state, "--keep-matrix-idmap requires --keep-matrix-format phylip");
			if (matrix_opt.keep_outf[0] != '\0' && matrix_opt.remove_outf[0] != '\0' &&
			    strcmp(matrix_opt.keep_outf, matrix_opt.remove_outf) == 0)
				argp_error(state, "--keep-out and --remove-out must be different files");
			if (matrix_opt.keep_outf[0] != '\0' && matrix_opt.outf[0] != '\0' &&
			    strcmp(matrix_opt.keep_outf, matrix_opt.outf) == 0)
				argp_error(state, "--keep-out must differ from --outfile");
			if (matrix_opt.remove_outf[0] != '\0' && matrix_opt.outf[0] != '\0' &&
			    strcmp(matrix_opt.remove_outf, matrix_opt.outf) == 0)
				argp_error(state, "--remove-out must differ from --outfile");
			if (matrix_opt.keep_matrix_outf[0] != '\0' && matrix_opt.outf[0] != '\0' &&
			    strcmp(matrix_opt.keep_matrix_outf, matrix_opt.outf) == 0)
				argp_error(state, "--keep-matrix-out must differ from --outfile");
			if (matrix_opt.matrix_idmap_outf[0] != '\0' && matrix_opt.outf[0] != '\0' &&
			    strcmp(matrix_opt.matrix_idmap_outf, matrix_opt.outf) == 0)
				argp_error(state, "--matrix-idmap must differ from --outfile");
			if (matrix_opt.edge_outf[0] != '\0' && matrix_opt.outf[0] != '\0' &&
			    strcmp(matrix_opt.edge_outf, matrix_opt.outf) == 0)
				argp_error(state, "--edge-out must differ from --outfile");
			if (matrix_opt.keep_outf[0] != '\0' && matrix_opt.edge_outf[0] != '\0' &&
			    strcmp(matrix_opt.keep_outf, matrix_opt.edge_outf) == 0)
				argp_error(state, "--keep-out must differ from --edge-out");
			if (matrix_opt.remove_outf[0] != '\0' && matrix_opt.edge_outf[0] != '\0' &&
			    strcmp(matrix_opt.remove_outf, matrix_opt.edge_outf) == 0)
				argp_error(state, "--remove-out must differ from --edge-out");
			if (matrix_opt.keep_matrix_outf[0] != '\0' && matrix_opt.edge_outf[0] != '\0' &&
			    strcmp(matrix_opt.keep_matrix_outf, matrix_opt.edge_outf) == 0)
				argp_error(state, "--keep-matrix-out must differ from --edge-out");
			if (matrix_opt.keep_matrix_outf[0] != '\0' && matrix_opt.keep_outf[0] != '\0' &&
			    strcmp(matrix_opt.keep_matrix_outf, matrix_opt.keep_outf) == 0)
				argp_error(state, "--keep-matrix-out must differ from --keep-out");
			if (matrix_opt.keep_matrix_outf[0] != '\0' && matrix_opt.remove_outf[0] != '\0' &&
			    strcmp(matrix_opt.keep_matrix_outf, matrix_opt.remove_outf) == 0)
				argp_error(state, "--keep-matrix-out must differ from --remove-out");
			if (matrix_opt.keep_matrix_idmap_outf[0] != '\0' && matrix_opt.outf[0] != '\0' &&
			    strcmp(matrix_opt.keep_matrix_idmap_outf, matrix_opt.outf) == 0)
				argp_error(state, "--keep-matrix-idmap must differ from --outfile");
			if (matrix_opt.keep_matrix_idmap_outf[0] != '\0' && matrix_opt.edge_outf[0] != '\0' &&
			    strcmp(matrix_opt.keep_matrix_idmap_outf, matrix_opt.edge_outf) == 0)
				argp_error(state, "--keep-matrix-idmap must differ from --edge-out");
			if (matrix_opt.keep_matrix_idmap_outf[0] != '\0' && matrix_opt.keep_outf[0] != '\0' &&
			    strcmp(matrix_opt.keep_matrix_idmap_outf, matrix_opt.keep_outf) == 0)
				argp_error(state, "--keep-matrix-idmap must differ from --keep-out");
			if (matrix_opt.keep_matrix_idmap_outf[0] != '\0' && matrix_opt.remove_outf[0] != '\0' &&
			    strcmp(matrix_opt.keep_matrix_idmap_outf, matrix_opt.remove_outf) == 0)
				argp_error(state, "--keep-matrix-idmap must differ from --remove-out");
			if (matrix_opt.keep_matrix_idmap_outf[0] != '\0' && matrix_opt.keep_matrix_outf[0] != '\0' &&
			    strcmp(matrix_opt.keep_matrix_idmap_outf, matrix_opt.keep_matrix_outf) == 0)
				argp_error(state, "--keep-matrix-idmap must differ from --keep-matrix-out");
			if (matrix_opt.matrix_idmap_outf[0] != '\0' && matrix_opt.edge_outf[0] != '\0' &&
			    strcmp(matrix_opt.matrix_idmap_outf, matrix_opt.edge_outf) == 0)
				argp_error(state, "--matrix-idmap must differ from --edge-out");
			if (matrix_opt.matrix_idmap_outf[0] != '\0' && matrix_opt.keep_outf[0] != '\0' &&
			    strcmp(matrix_opt.matrix_idmap_outf, matrix_opt.keep_outf) == 0)
				argp_error(state, "--matrix-idmap must differ from --keep-out");
			if (matrix_opt.matrix_idmap_outf[0] != '\0' && matrix_opt.remove_outf[0] != '\0' &&
			    strcmp(matrix_opt.matrix_idmap_outf, matrix_opt.remove_outf) == 0)
				argp_error(state, "--matrix-idmap must differ from --remove-out");
			if (matrix_opt.matrix_idmap_outf[0] != '\0' && matrix_opt.keep_matrix_outf[0] != '\0' &&
			    strcmp(matrix_opt.matrix_idmap_outf, matrix_opt.keep_matrix_outf) == 0)
				argp_error(state, "--matrix-idmap must differ from --keep-matrix-out");
			if (matrix_opt.matrix_idmap_outf[0] != '\0' && matrix_opt.keep_matrix_idmap_outf[0] != '\0' &&
			    strcmp(matrix_opt.matrix_idmap_outf, matrix_opt.keep_matrix_idmap_outf) == 0)
				argp_error(state, "--matrix-idmap must differ from --keep-matrix-idmap");
			break;
/*
			if(state->argc<2)
			{
      	printf("\v");
				argp_state_help(state,stdout,ARGP_HELP_SHORT_USAGE);
				printf("\v");
      	argp_state_help(state,stdout,ARGP_HELP_LONG);
      	printf("\v");
      	return EINVAL;
			}
*/
    }
    default:
      return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

static struct argp argp_matrix =
{
  opt_matrix,
  parse_matrix,
	0,//  "[arguments ...]",
  doc_matrix
};

int cmd_matrix(struct argp_state* state)
{
  struct arg_matrix matrix = { 0, };
  int    argc = state->argc - state->next + 1;
  char** argv = &state->argv[state->next - 1];
  matrix.global = state->input;
  argp_parse(&argp_matrix, argc, argv, ARGP_IN_ORDER, &argc, &matrix);

  state->next += argc - 1;
  if (matrix_opt.qrydir[0] != '\0') {
  	if (matrix_opt.refdir[0] == '\0')
	  	return compute_triangle(&matrix_opt);
		else if(matrix_opt.ani)
			return compute_ani_matrix(&matrix_opt);
		else
			return compute_matrix(&matrix_opt);
	}
  return 1;
}
