#include "command_sketch_wrapper.h"
#include "sketch_inspect.h"
#include "sketch_rearrange.h"
#include "global_basic.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <argp.h>
#include <argz.h>
#include <sys/stat.h>
#include <dirent.h>
#include <err.h>
#include <errno.h>
#include <unistd.h>
#include <stdbool.h>
#include <math.h>
#include <stdint.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

/*** argp wrapper ***/
struct arg_sketch
{
  struct arg_global *global;

  char *name;
  int manual_pattern_seen;
  int coden_pattern_seen;
  int outdir_seen;
  int dedup_metric_seen;
  int dedup_max_afcut_seen;
  int dedup_ctxcut_seen;
  int dedup_strategy_seen;
  int dedup_index_seen;
};

enum
{
  SKETCH_GROUP_IO = 1,
  SKETCH_GROUP_PRESET,
  SKETCH_GROUP_MANUAL,
  SKETCH_GROUP_FILTERS,
  SKETCH_GROUP_QC,
  SKETCH_GROUP_CONTENT,
  SKETCH_GROUP_LAYOUT,
  SKETCH_GROUP_INSPECT,
  SKETCH_GROUP_MODES,
  SKETCH_GROUP_MODIFY
};

enum
{
  SKETCH_PRINT_SAMPLES = 897,
  SKETCH_PRINT_SKETCH = 898,
  SKETCH_PRINT_INDEX = 899,
  SKETCH_PRINT_POSITIONS = 900,
  SKETCH_KEEP_SAMPLES = 901,
  SKETCH_DEDUP_SAMPLES = 902,
  SKETCH_DEDUP_METRIC = 903,
  SKETCH_DEDUP_MAX_AFCUT = 904,
  SKETCH_DEDUP_CTXCUT = 905,
  SKETCH_DEDUP_INDEX = 906,
  SKETCH_DEDUP_INDEX_MAX_CTX_FREQ = 907,
  SKETCH_DEDUP_INDEX_MIN_VOTES = 908,
  SKETCH_DEDUP_INDEX_SAMPLE_STEP = 909,
  SKETCH_DROP_POSITION = 910,
  SKETCH_DEDUP_STRATEGY = 911
};

static struct argp_option opt_sketch[] =
    {
        {0, 0, 0, 0, "Input, output, and execution:", SKETCH_GROUP_IO},
        {"list", 'l', "FILE", 0, "File containing input sequence paths.", SKETCH_GROUP_IO},
        {"outdir", 'o', "<path>", 0, "Output directory for result files.", SKETCH_GROUP_IO},
        {"pipecmd", 'P', "<cmd>", 0, "Stream each input through command; '{}' replaces the path.", SKETCH_GROUP_IO},
        {"threads", 'p', "<INT>", 0, "Number of threads to use. [1]", SKETCH_GROUP_IO},

        {0, 0, 0, 0, "Sketch preset and downsampling:", SKETCH_GROUP_PRESET},
        {"use_coden_ctxobj", 'T', 0, 0, "Use coden context-object structure pattern. [default]", SKETCH_GROUP_PRESET},
        {"DimRdcFold", 'f', "<INT>", 0, "K-mer space downsampling rate 1/2^f. [8]", SKETCH_GROUP_PRESET},

        {0, 0, 0, 0, "Manual context-object sizing:", SKETCH_GROUP_MANUAL},
        {"ctxlen", 'C', "<INT>", 0, "Half context length. [11]", SKETCH_GROUP_MANUAL},
        {"outerobjlen", 'O', "<INT>", 0, "Half outer object length. [0]", SKETCH_GROUP_MANUAL},
        {"innerobjlen", 'I', "<INT>", 0, "Inner object length. [0]", SKETCH_GROUP_MANUAL},

        {0, 0, 0, 0, "Read count thresholds:", SKETCH_GROUP_FILTERS},
        {"LstKmerOcrs", 'n', "INT", 0, "Least k-mer occurrence required in FASTQ input. [1]", SKETCH_GROUP_FILTERS},
        {"npercentile", 889, "<0..1>", 0, "Use the weighted lower-tail count percentile among k-mers passing -n; each k-mer contributes its count as weight.", SKETCH_GROUP_FILTERS},
        {"ncap", 891, "INT", OPTION_HIDDEN, "Cap the effective lower k-mer occurrence cutoff after -n, --npercentile, and --readsQC are combined; 0 disables. [0]", SKETCH_GROUP_FILTERS},

        {0, 0, 0, 0, "Read and sketch QC:", SKETCH_GROUP_QC},
        {"readsQC", 890, 0, 0, "Infer a count-weighted-median-anchored, negative-binomial-like k-mer count range for read QC and keep only that range.", SKETCH_GROUP_QC},
        {"sketchQC", 892, 0, 0, "Apply stored lcofiles.qc ranges to an existing abundance sketch.", SKETCH_GROUP_QC},

        {0, 0, 0, 0, "Sketch content:", SKETCH_GROUP_CONTENT},
        {"abundance", 'A', 0, 0, "Enable abundance-estimation mode.", SKETCH_GROUP_CONTENT},
        {"anno", 893, 0, 0, "Write FASTA/FASTQ header annotations to lcofiles.anno.", SKETCH_GROUP_CONTENT},
        {"nocomputemeta", 894, 0, 0, "Do not write per-input metadata sidecar lcofiles.infilemeta.", SKETCH_GROUP_CONTENT},
        {"position", 895, 0, 0, "Write zero-based sequence-stream positions to comblco.position.", SKETCH_GROUP_CONTENT},
        {"conflict", 666, 0, 0, "Keep conflicting context-objects for raw-read sketching.", SKETCH_GROUP_CONTENT},

        {0, 0, 0, 0, "Sample layout:", SKETCH_GROUP_LAYOUT},
        {"asone", 'a', 0, 0, "Treat input genomes as parts of one final genome.", SKETCH_GROUP_LAYOUT},
        {"splitmfa", 888, 0, 0, "Treat a multi-FASTA file as many genomes.", SKETCH_GROUP_LAYOUT},

        {0, 0, 0, 0, "Inspection modes:", SKETCH_GROUP_INSPECT},
        {"psmp", SKETCH_PRINT_SAMPLES, 0, 0, "Print sample names and sketch-entry counts.", SKETCH_GROUP_INSPECT},
        {"psketch", SKETCH_PRINT_SKETCH, 0, 0, "Print sketch content.", SKETCH_GROUP_INSPECT},
        {"pindex", SKETCH_PRINT_INDEX, 0, 0, "Print context/genome/object index content.", SKETCH_GROUP_INSPECT},
        {"ppos", SKETCH_PRINT_POSITIONS, 0, 0, "Print sketch positions as sample, sketch entry, and zero-based position.", SKETCH_GROUP_INSPECT},

        {0, 0, 0, 0, "Maintenance modes:", SKETCH_GROUP_MODES},
        {"index", 'i', "<FILE>", 0, "Build an inverted index for a combined sketch.", SKETCH_GROUP_MODES},
        {"merge", 777, 0, OPTION_HIDDEN, "Deprecated alias for --append copy mode.", SKETCH_GROUP_MODIFY},
        {"append", 778, 0, 0, "Append sketches; with -o writes a copy, without -o modifies the first sketch.", SKETCH_GROUP_MODIFY},
        {"remove", 896, "<FILE>", 0, "Remove listed samples; with -o writes a copy, without -o modifies the first sketch.", SKETCH_GROUP_MODIFY},
        {"keep", SKETCH_KEEP_SAMPLES, "<FILE>", 0, "Keep only listed samples; with -o writes a copy, without -o modifies the first sketch.", SKETCH_GROUP_MODIFY},
        {"drop-position", SKETCH_DROP_POSITION, 0, 0, "Drop comblco.position while writing --keep, --remove, or --dedup output.", SKETCH_GROUP_MODIFY},
        {"dedup", SKETCH_DEDUP_SAMPLES, "<DIST>", 0, "Deduplicate samples with distance < DIST; with -o can build from FASTA/FASTQ or write a sketch copy.", SKETCH_GROUP_MODIFY},
        {"dedup-strategy", SKETCH_DEDUP_STRATEGY, "<greedy|full-linkage>", 0, "Dedup grouping strategy: greedy representative-neighbor or full-linkage clique. [greedy]", SKETCH_GROUP_MODIFY},
        {"dedup-max-afcut", SKETCH_DEDUP_MAX_AFCUT, "<FLOAT>", 0, "Require max pairwise context alignment fraction for --dedup. [0.8]", SKETCH_GROUP_MODIFY},
        {"dedup-ctxcut", SKETCH_DEDUP_CTXCUT, "<INT>", 0, "Require at least this many shared contexts for --dedup. [0]", SKETCH_GROUP_MODIFY},
        {"dedup-index", SKETCH_DEDUP_INDEX, 0, 0, "Use an existing inverted index to nominate --dedup candidate pairs.", SKETCH_GROUP_MODIFY},
        {"dedup-index-max-ctx-freq", SKETCH_DEDUP_INDEX_MAX_CTX_FREQ, "<INT>", 0, "Indexed --dedup: skip context groups above this frequency. [256]", SKETCH_GROUP_MODIFY},
        {"dedup-index-min-votes", SKETCH_DEDUP_INDEX_MIN_VOTES, "<INT>", 0, "Indexed --dedup: exact-score pairs with at least this many sampled votes. [1]", SKETCH_GROUP_MODIFY},
        {"dedup-index-sample-step", SKETCH_DEDUP_INDEX_SAMPLE_STEP, "<INT>", 0, "Indexed --dedup: use every Nth query context for candidate nomination. [1]", SKETCH_GROUP_MODIFY},
        {"metric", SKETCH_DEDUP_METRIC, "<METRIC|EXPR>", 0, "Distance metric for --dedup; quoted A&B requires both, A|B allows either. [ctx-moe]", SKETCH_GROUP_MODIFY},
        {0}
    };
static char doc_sketch[] =
    "\n"
    "Create sketches from FASTA/FASTQ sequences, modify sketches, or build indexes."
    "\v"
    "Use '-' as one input to read FASTA/FASTQ from stdin.\n"
    "Use --pipecmd CMD to stream each input through a command; '{}' is replaced by the input path, otherwise the path is appended.\n"
    "\n"
    "Examples:\n"
    "  kssd3a sketch -o ref_sketches refs/*.fasta\n"
    "  kssd3a sketch -f8 -p8 -o qry_sketches queries/*.fq.gz\n"
    "  samtools fastq reads.bam | kssd3a sketch --conflict -o reads_sketch -\n"
    "  kssd3a sketch --pipecmd 'samtools fastq {}' --conflict \\\n"
    "    -o reads_sketch reads.bam\n"
    "  kssd3a sketch --psmp ref_sketches\n"
    "  kssd3a sketch --ppos positioned_sketches\n"
    "  kssd3a sketch -f8 -C9 -O3 -I3 -p8 -o custom_sketches genomes/*.fasta\n"
    "  kssd3a sketch -i ref_sketches\n"
    "  kssd3a sketch --append -o merged_sketches base_sketch add_sketch\n"
    "  kssd3a sketch --append base_sketch add_sketch\n"
    "  kssd3a sketch --remove remove_names.txt -o filtered_sketch ref_sketches\n"
    "  kssd3a sketch --remove remove_names.txt ref_sketches\n"
    "  kssd3a sketch --keep keep_names.txt -o kept_sketch ref_sketches\n"
    "  kssd3a sketch --keep keep_names.txt ref_sketches\n"
    "  kssd3a sketch --dedup 0.001 --metric ctx-moe -o dedup_ref ref_sketches\n"
    "  kssd3a sketch --dedup 0.0001 --metric ctx-naive --dedup-max-afcut 0.8 ref_sketches\n"
    "  kssd3a sketch --dedup 0.001 ref_sketches\n"
    "  kssd3a sketch --dedup 0.001 -o dedup_sketch genomes/*.fasta";

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

static double parse_double_range(struct argp_state *state, const char *option_name,
                                 const char *arg, double min_value, double max_value)
{
  char *end = NULL;
  errno = 0;
  double value = strtod(arg, &end);
  if (errno != 0 || end == arg || *end != '\0' || !isfinite(value) || value < min_value || value > max_value)
    argp_error(state, "%s requires a number in range %.6g..%.6g", option_name, min_value, max_value);
  return value;
}

static double parse_nonnegative_double(struct argp_state *state, const char *option_name,
                                       const char *arg)
{
  char *end = NULL;
  errno = 0;
  double value = strtod(arg, &end);
  if (errno != 0 || end == arg || *end != '\0' || !isfinite(value) || value < 0.0)
    argp_error(state, "%s requires a finite non-negative number", option_name);
  return value;
}

static sketch_dedup_metric_t parse_sketch_dedup_metric(struct argp_state *state,
                                                       const char *arg)
{
  pairwise_metric_expr_t metric = pairwise_metric_expr_single(PAIRWISE_METRIC_CTX_MOE);
  if (pairwise_metric_expr_from_string(arg, &metric))
    return metric;
  argp_error(state, "--metric must be one metric such as ctx-naive, p_dist, mash, or aaf, or an expression like 'ctx-naive&aaf' or 'ctx-moe|mash'");
  return pairwise_metric_expr_single(PAIRWISE_METRIC_CTX_MOE);
}

static void copy_path_arg(struct argp_state *state, const char *option_name,
                          char *dest, size_t dest_size, const char *arg)
{
  if (strlen(arg) >= dest_size)
    argp_error(state, "%s path is too long; maximum supported length is %zu bytes", option_name, dest_size - 1);
  snprintf(dest, dest_size, "%s", arg);
}

sketch_opt_t sketch_opt = {
    .hclen = 11, //
    .holen = 0,
    .iolen = 0,
    .drfold = 8,
    .kmerocrs = 1,
    .npercentile = 0.0,
    .ncap = 0,
    .reads_qc = false,
    .asone = 0,    // treat input genomes as parts of final genome.
    .p = 1,         // threads num: p
    .abundance = 0, // no abundance
    .conflict = 0, // no conflict context-objet.
    .anno = 0,
    .compute_meta = 1,
    .position = 0,
    .drop_position = false,
    .merge_comblco = 0,
    .append_comblco = 0,
    .remove_comblco = 0,
    .keep_comblco = 0,
    .dedup_comblco = 0,
    .dedup_raw_build_from_inputs = 0,
    .append_copy_mode = 0,
    .remove_copy_mode = 0,
    .keep_copy_mode = 0,
    .dedup_copy_mode = 0,
    .sketch_qc = 0,
    .dedup_cutoff = 0.0,
    .dedup_max_afcut = 0.8,
    .dedup_ctxcut = 0,
    .dedup_metric = {
        .metrics = { PAIRWISE_METRIC_CTX_MOE },
        .count = 1,
        .op = PAIRWISE_METRIC_EXPR_SINGLE,
        .label = "ctx-moe",
    },
    .dedup_strategy = PAIRWISE_DEDUP_GREEDY,
    .dedup_index = false,
    .dedup_index_max_ctx_freq = 256,
    .dedup_index_min_votes = 1,
    .dedup_index_sample_step = 1,
    .print_mode = 0,
    .split_mfa = 0,
    .coden_ctxobj_pattern = true, // coden ctxobj pattern by default
                                   //  .fpath[0] ='\0',
    .outdir = "./",
    .index[0] = '\0',
    .remove_list = NULL,
    .remove_source = NULL,
    .keep_list = NULL,
    .keep_source = NULL,
    .dedup_source = NULL,
    //	.pipecmd[0] = '\0', // no pipe command
    .num_remaining_args = 0, // int num_remaining_args; no option arguments num.
    .remaining_args = NULL   // char **remaining_args; no option arguments array.
};

static error_t parse_sketch(int key, char *arg, struct argp_state *state)
{
  struct arg_sketch *sketch = state->input;
  switch (key)
  {
  case 'C':
  {
    int val = parse_int_range(state, "-C/--ctxlen", arg, 4, 16);
    sketch->manual_pattern_seen = 1;
    sketch_opt.coden_ctxobj_pattern = false;
    sketch_opt.hclen = val;
    break;
  }
  case 'O':
  {
    int val = parse_int_range(state, "-O/--outerobjlen", arg, 0, 8);
    sketch->manual_pattern_seen = 1;
    sketch_opt.coden_ctxobj_pattern = false;
    sketch_opt.holen = val;
    break;
  }
  case 'I':
  {
    int val = parse_int_range(state, "-I/--innerobjlen", arg, 0, 8);
    sketch->manual_pattern_seen = 1;
    sketch_opt.coden_ctxobj_pattern = false;
    sketch_opt.iolen = val;
    break;
  }
  case 'a':
  {
    sketch_opt.asone = 1;
    break;
  }
  case 'f':
  {
    int val = parse_int_range(state, "-f/--DimRdcFold", arg, 0, 24);
    sketch_opt.drfold = val;
    break;
  }
  case 'n':
  {
    int val = parse_int_range(state, "-n/--LstKmerOcrs", arg, 1, 65536);
    sketch_opt.kmerocrs = val;
    break;
  }
  case 889:
  {
    sketch_opt.npercentile = parse_double_range(state, "--npercentile", arg, 0.0, 1.0);
    break;
  }
  case 890:
  {
    sketch_opt.reads_qc = true;
    if (sketch_opt.kmerocrs < 2)
      sketch_opt.kmerocrs = 2;
    break;
  }
  case 891:
  {
    sketch_opt.ncap = parse_int_range(state, "--ncap", arg, 0, 65536);
    break;
  }
  case 'i':
  {
    copy_path_arg(state, "-i/--index", sketch_opt.index, sizeof(sketch_opt.index), arg);
    break;
  }
  case 'p':
  {
#ifdef _OPENMP
    sketch_opt.p = parse_int_range(state, "-p/--threads", arg, 1, 65536);
#else
    int ignored_threads = parse_int_range(state, "-p/--threads", arg, 1, 65536);
    warnx("This version of kssd was built without OpenMP and "
          "thus does not support multi threading. Ignoring -p %d",
          ignored_threads);
#endif
    break;
  }
  case 'A':
  {
    sketch_opt.abundance = 1;
    break;
  }
  case 893:
  {
    sketch_opt.anno = 1;
    break;
  }
  case 894:
  {
    sketch_opt.compute_meta = 0;
    break;
  }
  case 895:
  {
    sketch_opt.position = 1;
    break;
  }
  case 'T':
  {
    sketch->coden_pattern_seen = 1;
    sketch_opt.coden_ctxobj_pattern = 1;
    break;
  }
  case 'P':
  {
    sketch_opt.pipecmd = malloc(strlen(arg) + 1);
    if (sketch_opt.pipecmd == NULL)
      err(EXIT_FAILURE, "%s(): failed to allocate pipe command", __func__);
    strcpy(sketch_opt.pipecmd, arg);
    break;
  }
  case 'l':
  {
    sketch_opt.fpath = malloc(strlen(arg) + 1);
    if (sketch_opt.fpath == NULL)
      err(EXIT_FAILURE, "%s(): failed to allocate input list path", __func__);
    strcpy(sketch_opt.fpath, arg);
    break;
  }
  case 'o':
  {
    sketch_opt.outdir = malloc(strlen(arg) + 10);
    strcpy(sketch_opt.outdir, arg);
    sketch->outdir_seen = 1;
    break;
  }
  case 666:
  {
    sketch_opt.conflict = 1; // keep conflict context-objet(for raw reads sketching).
    break;
  }
  case 777:
  {
    sketch_opt.append_comblco = 1;
    sketch_opt.append_copy_mode = 1;
    break;
  }
  case 778:
  {
    sketch_opt.append_comblco = 1;
    break;
  }
  case 896:
  {
    sketch_opt.remove_comblco = 1;
    sketch_opt.remove_list = malloc(strlen(arg) + 1);
    if (sketch_opt.remove_list == NULL)
      err(EXIT_FAILURE, "%s(): failed to allocate remove list path", __func__);
    strcpy(sketch_opt.remove_list, arg);
    break;
  }
  case SKETCH_KEEP_SAMPLES:
  {
    sketch_opt.keep_comblco = 1;
    sketch_opt.keep_list = malloc(strlen(arg) + 1);
    if (sketch_opt.keep_list == NULL)
      err(EXIT_FAILURE, "%s(): failed to allocate keep list path", __func__);
    strcpy(sketch_opt.keep_list, arg);
    break;
  }
  case SKETCH_DEDUP_SAMPLES:
  {
    sketch_opt.dedup_comblco = 1;
    sketch_opt.dedup_cutoff = parse_nonnegative_double(state, "--dedup", arg);
    break;
  }
  case SKETCH_DEDUP_METRIC:
  {
    sketch->dedup_metric_seen = 1;
    sketch_opt.dedup_metric = parse_sketch_dedup_metric(state, arg);
    break;
  }
  case SKETCH_DEDUP_STRATEGY:
  {
    sketch->dedup_strategy_seen = 1;
    if (!pairwise_dedup_strategy_from_string(arg, &sketch_opt.dedup_strategy))
      argp_error(state, "--dedup-strategy must be one of greedy, full-linkage, complete-linkage, or clique");
    break;
  }
  case SKETCH_DEDUP_MAX_AFCUT:
  {
    sketch->dedup_max_afcut_seen = 1;
    sketch_opt.dedup_max_afcut = parse_double_range(state, "--dedup-max-afcut", arg, 0.0, 1.0);
    break;
  }
  case SKETCH_DEDUP_CTXCUT:
  {
    sketch->dedup_ctxcut_seen = 1;
    sketch_opt.dedup_ctxcut = (uint32_t)parse_int_range(state, "--dedup-ctxcut", arg, 0, INT32_MAX);
    break;
  }
  case SKETCH_DEDUP_INDEX:
  {
    sketch->dedup_index_seen = 1;
    sketch_opt.dedup_index = true;
    break;
  }
  case SKETCH_DEDUP_INDEX_MAX_CTX_FREQ:
  {
    sketch->dedup_index_seen = 1;
    sketch_opt.dedup_index = true;
    sketch_opt.dedup_index_max_ctx_freq =
        (uint32_t)parse_int_range(state, "--dedup-index-max-ctx-freq", arg, 0, INT32_MAX);
    break;
  }
  case SKETCH_DEDUP_INDEX_MIN_VOTES:
  {
    sketch->dedup_index_seen = 1;
    sketch_opt.dedup_index = true;
    sketch_opt.dedup_index_min_votes =
        (uint32_t)parse_int_range(state, "--dedup-index-min-votes", arg, 1, INT32_MAX);
    break;
  }
  case SKETCH_DEDUP_INDEX_SAMPLE_STEP:
  {
    sketch->dedup_index_seen = 1;
    sketch_opt.dedup_index = true;
    sketch_opt.dedup_index_sample_step =
        (uint32_t)parse_int_range(state, "--dedup-index-sample-step", arg, 1, INT32_MAX);
    break;
  }
  case SKETCH_DROP_POSITION:
  {
    sketch_opt.drop_position = true;
    break;
  }
  case 888:
  {
    sketch_opt.split_mfa = 1;
    break;
  }
  case 892:
  {
    sketch_opt.sketch_qc = 1;
    break;
  }
  case SKETCH_PRINT_SAMPLES:
  {
    sketch_opt.print_mode = 1;
    break;
  }
  case SKETCH_PRINT_SKETCH:
  {
    sketch_opt.print_mode = 2;
    break;
  }
  case SKETCH_PRINT_INDEX:
  {
    sketch_opt.print_mode = 3;
    break;
  }
  case SKETCH_PRINT_POSITIONS:
  {
    sketch_opt.print_mode = 4;
    break;
  }
  case ARGP_KEY_ARGS:
  {
    sketch_opt.num_remaining_args = state->argc - state->next;
    sketch_opt.remaining_args = state->argv + state->next;
    break;
  }
  case ARGP_KEY_END:
  {
    int mode_count = (sketch_opt.merge_comblco ? 1 : 0)
                     + (sketch_opt.append_comblco ? 1 : 0)
                     + (sketch_opt.remove_comblco ? 1 : 0)
                     + (sketch_opt.keep_comblco ? 1 : 0)
                     + (sketch_opt.dedup_comblco ? 1 : 0)
                     + (sketch_opt.sketch_qc ? 1 : 0)
                     + (sketch_opt.print_mode ? 1 : 0)
                     + (sketch_opt.index[0] != '\0' ? 1 : 0);
    if (mode_count > 1)
      argp_error(state, "Use only one of --merge, --append, --remove, --keep, --dedup, --sketchQC, --psmp/--psketch/--pindex/--ppos, or -i/--index.");
    if (sketch->dedup_metric_seen && !sketch_opt.dedup_comblco)
      argp_error(state, "--metric is currently only used with --dedup.");
    if (sketch->dedup_strategy_seen && !sketch_opt.dedup_comblco)
      argp_error(state, "--dedup-strategy is currently only used with --dedup.");
    if (!sketch_opt.dedup_comblco &&
        (sketch->dedup_max_afcut_seen || sketch->dedup_ctxcut_seen))
      argp_error(state, "--dedup-max-afcut and --dedup-ctxcut are currently only used with --dedup.");
    if (sketch->dedup_index_seen && !sketch_opt.dedup_comblco)
      argp_error(state, "--dedup-index options are currently only used with --dedup.");
    if (sketch_opt.drop_position &&
        !(sketch_opt.remove_comblco || sketch_opt.keep_comblco || sketch_opt.dedup_comblco))
      argp_error(state, "--drop-position is only used with --remove, --keep, or --dedup.");
    if (sketch_opt.append_comblco)
    {
      if (sketch->outdir_seen)
        sketch_opt.append_copy_mode = 1;
      if (sketch_opt.append_copy_mode)
      {
        if (!sketch->outdir_seen)
          argp_error(state, "--merge/--append copy mode requires -o/--outdir to name the output sketch directory.");
        if (sketch_opt.num_remaining_args < 2)
          argp_error(state, "--append copy mode requires a base sketch followed by at least one sketch to append.");
      }
      else
      {
        if (sketch_opt.num_remaining_args < 2)
          argp_error(state, "--append in-place mode requires a target sketch followed by at least one source sketch.");
        sketch_opt.outdir = sketch_opt.remaining_args[0];
        sketch_opt.remaining_args++;
        sketch_opt.num_remaining_args--;
      }
    }
    if (sketch_opt.remove_comblco)
    {
      if (sketch_opt.remove_list == NULL || sketch_opt.remove_list[0] == '\0')
        argp_error(state, "--remove requires a newline-delimited sample-name list file.");
      if (sketch->outdir_seen)
      {
        sketch_opt.remove_copy_mode = 1;
        if (sketch_opt.num_remaining_args != 1)
          argp_error(state, "--remove copy mode requires exactly one source sketch directory after the list file.");
        sketch_opt.remove_source = sketch_opt.remaining_args[0];
      }
      else
      {
        if (sketch_opt.num_remaining_args != 1)
          argp_error(state, "--remove in-place mode requires exactly one target sketch directory after the list file.");
        sketch_opt.outdir = sketch_opt.remaining_args[0];
        sketch_opt.remove_source = sketch_opt.outdir;
      }
    }
    if (sketch_opt.keep_comblco)
    {
      if (sketch_opt.keep_list == NULL || sketch_opt.keep_list[0] == '\0')
        argp_error(state, "--keep requires a newline-delimited sample-name list file.");
      if (sketch->outdir_seen)
      {
        sketch_opt.keep_copy_mode = 1;
        if (sketch_opt.num_remaining_args != 1)
          argp_error(state, "--keep copy mode requires exactly one source sketch directory after the list file.");
        sketch_opt.keep_source = sketch_opt.remaining_args[0];
      }
      else
      {
        if (sketch_opt.num_remaining_args != 1)
          argp_error(state, "--keep in-place mode requires exactly one target sketch directory after the list file.");
        sketch_opt.outdir = sketch_opt.remaining_args[0];
        sketch_opt.keep_source = sketch_opt.outdir;
      }
    }
    if (sketch_opt.dedup_comblco)
    {
      if (sketch->outdir_seen)
      {
        if (sketch_opt.fpath != NULL ||
            sketch_opt.num_remaining_args != 1 ||
            !file_exists_in_folder(sketch_opt.remaining_args[0], sketch_stat))
        {
          sketch_opt.dedup_raw_build_from_inputs = 1;
          if (sketch_opt.fpath == NULL && sketch_opt.num_remaining_args == 0)
            argp_error(state, "--dedup -o raw-build mode requires FASTA/FASTQ inputs or -l/--list.");
        }
        else
        {
          sketch_opt.dedup_copy_mode = 1;
          sketch_opt.dedup_source = sketch_opt.remaining_args[0];
        }
      }
      else
      {
        if (sketch_opt.num_remaining_args != 1)
          argp_error(state, "--dedup in-place mode requires exactly one target sketch directory.");
        sketch_opt.outdir = sketch_opt.remaining_args[0];
        sketch_opt.dedup_source = sketch_opt.outdir;
      }
    }
    if (sketch_opt.print_mode)
    {
      if (sketch->outdir_seen)
        argp_error(state, "--psmp/--psketch/--pindex/--ppos print to stdout and do not use -o/--outdir.");
      if (sketch_opt.fpath != NULL)
        argp_error(state, "--psmp/--psketch/--pindex/--ppos take one sketch directory argument, not -l/--list.");
      if (sketch_opt.num_remaining_args != 1)
        argp_error(state, "--psmp/--psketch/--pindex/--ppos require exactly one sketch directory.");
    }
    if (sketch->coden_pattern_seen && sketch->manual_pattern_seen)
    {
      argp_error(state, "-T/--use_coden_ctxobj selects the coden pattern; do not combine it with -C, -O, or -I. Use either -T or the manual -C/-O/-I pattern.");
    }
    int klen = sketch_opt.iolen + 2 * (sketch_opt.holen + sketch_opt.hclen);
    if (klen > 32)
    {
      printf("\nError: k-mer length %d should smaller than 32 \n\n", klen);
      exit(1);
    }
    if (mode_count == 0 && sketch_opt.index[0] == '\0' && sketch_opt.fpath == NULL && sketch_opt.num_remaining_args == 0)
    {
      printf("\nError: missing input sequences file \n\n");
      argp_state_help(state, stdout, ARGP_HELP_STD_HELP);
      argp_usage(state);
    }
    break;
  }
  default:
    return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

static struct argp argp_sketch =
    {
        opt_sketch,
        parse_sketch,
        0, //  "[arguments ...]",
        doc_sketch};

infile_tab_t *sketch_organize_infiles(sketch_opt_t *sketch_opt_val)
{
  int fmt_ck;
  if (sketch_opt_val->pipecmd == NULL)
    fmt_ck = 1; // need check format-- normal mode
  else
    fmt_ck = 0;
  // do it if fpath is not ""
  if (sketch_opt_val->fpath != NULL)
  {
    return organize_infile_list(sketch_opt_val->fpath, fmt_ck);
  }
  else if (sketch_opt_val->num_remaining_args > 0)
  {
    return organize_infile_frm_arg(sketch_opt_val->num_remaining_args, sketch_opt_val->remaining_args, fmt_ck);
  }
  else
  {
    perror("please specify the (meta)genome files");
    return NULL;
  }
};

static int sketch_stdin_input_count(infile_tab_t *infile_stat)
{
  int count = 0;
  if (!infile_stat)
    return 0;
  for (int i = 0; i < infile_stat->infile_num; ++i)
    if (strcmp(infile_stat->organized_infile_tab[i].fpath, "-") == 0)
      count++;
  return count;
}

extern uint32_t FILTER; // control dimensionality reducation level
extern uint32_t hash_id;
extern dim_sketch_stat_t comblco_stat_one;
extern void compute_sketch(sketch_opt_t *, infile_tab_t *);
extern void gen_inverted_index4comblco(const char *sketchdir);
extern int merge_comblco(sketch_opt_t *sketch_opt_val);
extern int append_comblco(sketch_opt_t *sketch_opt_val);
extern int remove_comblco_samples(sketch_opt_t *sketch_opt_val);
extern int keep_comblco_samples(sketch_opt_t *sketch_opt_val);
extern int dedup_comblco_samples(sketch_opt_t *sketch_opt_val);
extern int sketch_qc_comblco(sketch_opt_t *sketch_opt_val);
extern uint32_t get_sketching_id(uint32_t hclen, uint32_t holen, uint32_t iolen, uint32_t drfold, uint32_t FILTER);

static void sketch_free_infile_tab(infile_tab_t *infile_stat)
{
  if (!infile_stat)
    return;
  for (int i = 0; i < infile_stat->infile_num; ++i)
    free(infile_stat->organized_infile_tab[i].fpath);
  free(infile_stat->organized_infile_tab);
  free(infile_stat);
}

static bool sketch_valid_lco_dir(const char *path)
{
  return file_exists_in_folder(path, sketch_stat) &&
         file_exists_in_folder(path, idx_sketch_suffix) &&
         file_exists_in_folder(path, combined_sketch_suffix);
}

static bool sketch_dir_is_empty(const char *path)
{
  DIR *dir = opendir(path);
  if (!dir)
    return false;
  struct dirent *ent;
  while ((ent = readdir(dir)) != NULL)
  {
    if (strcmp(ent->d_name, ".") != 0 && strcmp(ent->d_name, "..") != 0)
    {
      closedir(dir);
      return false;
    }
  }
  closedir(dir);
  return true;
}

static void sketch_remove_tree(const char *path)
{
  struct stat st;
  if (lstat(path, &st) != 0)
    return;
  if (!S_ISDIR(st.st_mode))
  {
    if (unlink(path) != 0)
      err(errno, "%s(): cannot remove %s", __func__, path);
    return;
  }

  DIR *dir = opendir(path);
  if (!dir)
    err(errno, "%s(): cannot open %s", __func__, path);
  struct dirent *ent;
  while ((ent = readdir(dir)) != NULL)
  {
    if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0)
      continue;
    char *child = format_string("%s/%s", path, ent->d_name);
    if (!child)
      err(EXIT_FAILURE, "%s(): OOM remove path", __func__);
    sketch_remove_tree(child);
    free(child);
  }
  closedir(dir);
  if (rmdir(path) != 0)
    err(errno, "%s(): cannot remove %s", __func__, path);
}

static void sketch_fail_if_path_exists(const char *path)
{
  struct stat st;
  if (lstat(path, &st) == 0)
    errx(EINVAL, "%s(): temporary path already exists: %s", __func__, path);
  if (errno != ENOENT)
    err(errno, "%s(): cannot stat %s", __func__, path);
}

static char *sketch_raw_build_tmp_dir(const char *outdir, const char *label)
{
  char *path = format_string("%s.dedup_build.%ld.%s", outdir, (long)getpid(), label);
  if (!path)
    err(EXIT_FAILURE, "%s(): OOM temporary path", __func__);
  sketch_fail_if_path_exists(path);
  return path;
}

static void sketch_check_raw_dedup_output(const char *outdir)
{
  struct stat st;
  if (lstat(outdir, &st) != 0)
  {
    if (errno == ENOENT)
      return;
    err(errno, "%s(): cannot stat %s", __func__, outdir);
  }
  if (!S_ISDIR(st.st_mode))
    errx(EINVAL, "%s(): output path exists but is not a directory: %s", __func__, outdir);
  if (sketch_valid_lco_dir(outdir))
    errx(EINVAL,
         "%s(): --dedup -o cannot exactly update an existing sketch from raw inputs: %s. Exact dedup requires all candidate samples in one sketch; append into a temporary combined sketch and run --dedup on that sketch.",
         __func__, outdir);
  if (sketch_dir_is_empty(outdir))
    return;
  errx(EINVAL, "%s(): output directory exists but is neither an empty directory nor a valid sketch: %s",
       __func__, outdir);
}

static void sketch_promote_tmp_dir(const char *tmp_dir, const char *outdir)
{
  struct stat st;
  const bool out_exists = lstat(outdir, &st) == 0;
  if (!out_exists && errno != ENOENT)
    err(errno, "%s(): cannot stat %s", __func__, outdir);
  if (out_exists && !S_ISDIR(st.st_mode))
    errx(EINVAL, "%s(): output path exists but is not a directory: %s", __func__, outdir);

  char *backup = NULL;
  if (out_exists)
  {
    backup = format_string("%s.dedup_build.%ld.backup", outdir, (long)getpid());
    if (!backup)
      err(EXIT_FAILURE, "%s(): OOM backup path", __func__);
    sketch_fail_if_path_exists(backup);
    if (rename(outdir, backup) != 0)
      err(errno, "%s(): cannot move existing %s to %s", __func__, outdir, backup);
  }

  if (rename(tmp_dir, outdir) != 0)
  {
    const int saved_errno = errno;
    if (backup && rename(backup, outdir) != 0)
      warn("%s(): rollback failed moving %s back to %s", __func__, backup, outdir);
    errno = saved_errno;
    err(errno, "%s(): cannot promote %s to %s", __func__, tmp_dir, outdir);
  }

  if (backup)
  {
    sketch_remove_tree(backup);
    free(backup);
  }
}

static void sketch_prepare_output_stat(sketch_opt_t *opt, infile_tab_t *infile_stat)
{
  FILTER = UINT32_MAX >> opt->drfold;
  if (opt->coden_ctxobj_pattern)
  {
    hash_id = get_sketching_id(NUM_CODENS, 0, 0, opt->drfold, FILTER);
#if NUM_CODENS < 11
    klen = 3 * NUM_CODENS + 1;
#else
    klen = 32;
#endif
    comblco_stat_one.coden_len = NUM_CODENS;
    comblco_stat_one.hclen = 0;
    comblco_stat_one.holen = 0;
  }
  else
  {
    hash_id = get_sketching_id(opt->hclen, opt->holen, opt->iolen, opt->drfold, FILTER);
    klen = 2 * (opt->hclen + opt->holen) + opt->iolen;
    comblco_stat_one.coden_len = 0;
    comblco_stat_one.hclen = opt->hclen;
    comblco_stat_one.holen = opt->holen;
  }
  if (NUM_CODENS > 11 || klen > 32 || FILTER < 256)
    err(EINVAL, "%s(): NUM_CODENS(%d) or klen (%d) or FILTER (%u) is out of range (NUM_CODENS <=11 and klen <=32 and FILTER: 256..0xffffffff)",
        __func__, NUM_CODENS, klen, FILTER);

  printf("Sketching method hashid = %u\tctxobj_coden_len=%u\tklen=%u\tFILTER=%u\thclen=%d\n",
         hash_id, comblco_stat_one.coden_len, klen, FILTER, comblco_stat_one.hclen);
  comblco_stat_one.hash_id = hash_id;
  comblco_stat_one.koc = opt->abundance;
  comblco_stat_one.conflict = opt->conflict;
  comblco_stat_one.klen = klen;
  comblco_stat_one.drfold = opt->drfold;
  comblco_stat_one.infile_num = opt->asone
                                    ? (infile_stat->infile_num < 1 ? 0 : 1)
                                    : infile_stat->infile_num;
  const_comask_init(&comblco_stat_one);
  set_uint64kmer2generic_ctxobj(opt->coden_ctxobj_pattern);
}

static void sketch_compute_inputs_to_dir(sketch_opt_t *opt, infile_tab_t *infile_stat,
                                         const char *outdir)
{
  char *old_outdir = opt->outdir;
  opt->outdir = (char *)outdir;
  sketch_prepare_output_stat(opt, infile_stat);
  mkdir_p(opt->outdir);
  compute_sketch(opt, infile_stat);
  opt->outdir = old_outdir;
}

static int sketch_dedup_raw_build_from_inputs(sketch_opt_t *opt)
{
  infile_tab_t *infile_stat = sketch_organize_infiles(opt);
  if (!infile_stat || infile_stat->infile_num < 1)
    errx(EXIT_FAILURE, "not valid fas/fastq files!");

  int stdin_count = sketch_stdin_input_count(infile_stat);
  if (stdin_count > 1)
    errx(EXIT_FAILURE, "stdin input '-' can be used only once");
  if ((stdin_count > 0 || opt->pipecmd != NULL) && opt->split_mfa)
    errx(EXIT_FAILURE, "--splitmfa does not support '-' or --pipecmd streaming inputs");

  sketch_check_raw_dedup_output(opt->outdir);
  char *new_dir = sketch_raw_build_tmp_dir(opt->outdir, "new");
  char *dedup_dir = sketch_raw_build_tmp_dir(opt->outdir, "dedup");

  sketch_opt_t new_opt = *opt;
  new_opt.dedup_comblco = false;
  new_opt.dedup_raw_build_from_inputs = false;
  new_opt.dedup_copy_mode = false;
  sketch_compute_inputs_to_dir(&new_opt, infile_stat, new_dir);
  if (opt->dedup_index)
    gen_inverted_index4comblco(new_dir);

  sketch_opt_t dedup_opt = {
      .outdir = dedup_dir,
      .p = opt->p,
      .dedup_source = new_dir,
      .dedup_copy_mode = true,
      .drop_position = opt->drop_position,
      .dedup_cutoff = opt->dedup_cutoff,
      .dedup_max_afcut = opt->dedup_max_afcut,
      .dedup_ctxcut = opt->dedup_ctxcut,
      .dedup_metric = opt->dedup_metric,
      .dedup_index = opt->dedup_index,
      .dedup_index_max_ctx_freq = opt->dedup_index_max_ctx_freq,
      .dedup_index_min_votes = opt->dedup_index_min_votes,
      .dedup_index_sample_step = opt->dedup_index_sample_step,
  };
  const int kept_samples = dedup_comblco_samples(&dedup_opt);

  sketch_promote_tmp_dir(dedup_dir, opt->outdir);
  sketch_remove_tree(new_dir);

  printf("Built deduplicated sketch %s from %d input samples; total samples=%d\n",
         opt->outdir, infile_stat->infile_num, kept_samples);

  sketch_free_infile_tab(infile_stat);
  free(new_dir);
  free(dedup_dir);
  return kept_samples;
}

int cmd_sketch(struct argp_state *state)
{
  struct arg_sketch sketch = {
      0,
  };
  int argc = state->argc - state->next + 1;
  char **argv = &state->argv[state->next - 1];
  sketch.global = state->input;
  argp_parse(&argp_sketch, argc, argv, ARGP_IN_ORDER, &argc, &sketch);
  state->next += argc - 1;

  if (sketch_opt.merge_comblco)
  {
    int merge_count = merge_comblco(&sketch_opt);
  }
  else if (sketch_opt.print_mode)
  {
    const char *sketch_path = sketch_opt.remaining_args[0];
    if (!file_exists_in_folder(sketch_path, sketch_stat))
      errx(EXIT_FAILURE, "%s is not a valid sketch", sketch_path);
    if (sketch_opt.print_mode == 1)
      sketch_inspect_print_samples(sketch_path);
    else
      sketch_inspect_print_content(sketch_path, sketch_opt.print_mode - 1);
  }
  else if (sketch_opt.append_comblco)
  {
    int append_count = append_comblco(&sketch_opt);
  }
  else if (sketch_opt.remove_comblco)
  {
    int remove_count = remove_comblco_samples(&sketch_opt);
  }
  else if (sketch_opt.keep_comblco)
  {
    int keep_count = keep_comblco_samples(&sketch_opt);
  }
  else if (sketch_opt.dedup_comblco)
  {
    int dedup_count = sketch_opt.dedup_raw_build_from_inputs
                          ? sketch_dedup_raw_build_from_inputs(&sketch_opt)
                          : dedup_comblco_samples(&sketch_opt);
  }
  else if (sketch_opt.sketch_qc)
  {
    return sketch_qc_comblco(&sketch_opt);
  }
  else if (sketch_opt.index[0] != '\0')
  {
    gen_inverted_index4comblco(sketch_opt.index);
  }
  else
  {
    infile_tab_t *infile_stat = sketch_organize_infiles(&sketch_opt);
    if (infile_stat->infile_num)
    {
      int stdin_count = sketch_stdin_input_count(infile_stat);
      if (stdin_count > 1)
        errx(EXIT_FAILURE, "stdin input '-' can be used only once");
      if ((stdin_count > 0 || sketch_opt.pipecmd != NULL) && sketch_opt.split_mfa)
        errx(EXIT_FAILURE, "--splitmfa does not support '-' or --pipecmd streaming inputs");
      FILTER = UINT32_MAX >> sketch_opt.drfold;
      /* conditionally initilize some comblco_stat_one members*/
      if(sketch_opt.coden_ctxobj_pattern){
        hash_id = get_sketching_id(NUM_CODENS,0,0,sketch_opt.drfold, FILTER);
#if NUM_CODENS < 11
        klen = 3 * NUM_CODENS + 1 ; // klen is 3*NUM_CODENS+1
#else
        klen = 32 ;
#endif
//  klen = NUM_CODENS < 11 ? 3 * NUM_CODENS + 1 : 32  ; // klen is 3*NUM_CODENS+1 only when NUM_CODENS <=10
        comblco_stat_one.coden_len = NUM_CODENS; // set coden_len
        comblco_stat_one.hclen = 0;
        comblco_stat_one.holen = 0;
      }else{
        hash_id =  get_sketching_id(sketch_opt.hclen, sketch_opt.holen, sketch_opt.iolen, sketch_opt.drfold, FILTER);
        klen = 2 * (sketch_opt.hclen + sketch_opt.holen) + sketch_opt.iolen;
        comblco_stat_one.coden_len = 0; // no coden ctxobj pattern
        comblco_stat_one.hclen = sketch_opt.hclen;
        comblco_stat_one.holen = sketch_opt.holen;
      }
      if (NUM_CODENS > 11 || klen > 32 || FILTER < 256)
        err(EINVAL, "%s(): NUM_CODENS(%d) or klen (%d) or FILTER (%u) is out of range (NUM_CODENS <=11 and klen <=32 and FILTER: 256..0xffffffff)", __func__, NUM_CODENS, klen, FILTER);

      printf("Sketching method hashid = %u\tctxobj_coden_len=%u\tklen=%u\tFILTER=%u\thclen=%d\n", hash_id, comblco_stat_one.coden_len, klen, FILTER, comblco_stat_one.hclen);
      { /*initilize the rest comblco_stat_one member*/
        comblco_stat_one.hash_id = hash_id;
        comblco_stat_one.koc = sketch_opt.abundance;
        comblco_stat_one.conflict = sketch_opt.conflict;
        comblco_stat_one.klen = klen;     
        comblco_stat_one.drfold = sketch_opt.drfold;
        if(sketch_opt.asone)
          comblco_stat_one.infile_num =  infile_stat->infile_num < 1 ? 0: 1;
        else
          comblco_stat_one.infile_num = infile_stat->infile_num ;
      }
      // ensure initalize global masks
      const_comask_init(&comblco_stat_one);
      mkdir_p(sketch_opt.outdir);
      // initialize k-mer rearrange method: reorder_unituple_by_coden_pattern64() or uint64_kmer2ctxobj()
      set_uint64kmer2generic_ctxobj(sketch_opt.coden_ctxobj_pattern);
      compute_sketch(&sketch_opt, infile_stat);
    }
    else
    {
      printf("not valid fas/fastq files!\n");
    }
    free(infile_stat->organized_infile_tab);
  }
  return 1;
};
