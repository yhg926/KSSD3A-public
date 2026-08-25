#include "command_ani.h"
#include "command_sketch_wrapper.h"
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
#include <sys/stat.h>
#include <unistd.h>

/*** argp wrapper ***/
struct arg_ani
{
	struct arg_global *global;

	char *name;
	int manual_pattern_seen;
	int coden_pattern_seen;
};

enum
{
	ANI_SKETCH_USE_CODEN = 781,
	ANI_SKETCH_CTXLEN,
	ANI_SKETCH_OUTEROBJLEN,
	ANI_SKETCH_INNEROBJLEN,
	ANI_SKETCH_DRFOLD,
	ANI_SKETCH_LST_KMER_OCRS,
	ANI_SKETCH_NPERCENTILE,
	ANI_SKETCH_NCAP,
	ANI_SKETCH_READS_QC,
	ANI_SKETCH_ABUNDANCE,
	ANI_SKETCH_ANNO,
	ANI_SKETCH_CONFLICT,
	ANI_SKETCH_ASONE,
	ANI_SKETCH_SPLITMFA,
	ANI_SKETCH_PIPECMD,
	ANI_RAW_OUTPUT,
	ANI_UNIFIED_METRIC
};

enum
{
	ANI_GROUP_INPUT = 1,
	ANI_GROUP_MODE,
	ANI_GROUP_FILTER,
	ANI_GROUP_REPORT,
	ANI_GROUP_EXECUTION,
	ANI_GROUP_AUTO_MODE,
	ANI_GROUP_AUTOSKETCH
};

static struct argp_option opt_ani[] =
	{
		{0, 0, 0, 0, "Input sketches and auto-ANI lists:", ANI_GROUP_INPUT},
		{"ref", 'r', "<SKETCH>", 0, "Reference sketch directory.", ANI_GROUP_INPUT},
		{"query", 'q', "<SKETCH|FASTA|FASTQ[.gz]>", 0, "Query sketch directory or sequence file.", ANI_GROUP_INPUT},
		{"qraw", 'x', "<SKETCH|FASTA|FASTQ[.gz]>", 0, "Raw-read query sketch directory or sequence file; enables naive ANI mode.", ANI_GROUP_INPUT},
		{"reflist", 779, "<FILE>", 0, "Reference FASTA/FASTQ/sketch path list for auto ANI.", ANI_GROUP_INPUT},
		{"qrylist", 780, "<FILE>", 0, "Query FASTA/FASTQ/sketch path list for auto ANI.", ANI_GROUP_INPUT},
		{"index", 'i', "<FILE>", OPTION_HIDDEN, "Reserved; currently unused.", 2},
		//		{"model", 'M', "<FILE>", 0, "specify the trained ani model.", 2},

		{0, 0, 0, 0, "ANI mode and conflict handling:", ANI_GROUP_MODE},
		{"naive", 'v', 0, 0, "Use naive distance/ANI calculation.", ANI_GROUP_MODE},
		{"ignoreconflict", 778, 0, 0, "Ignore reference-side contexts that contain conflicting objects.", ANI_GROUP_MODE},

		{0, 0, 0, 0, "Filtering and metrics:", ANI_GROUP_FILTER},
		{"afcut", 'f', "<FLOAT>", 0, "Skip reports with max(Qry_align_fraction, Ref_align_fraction) below this value. [0.5; 0.2 for unassembled/qraw queries]", ANI_GROUP_FILTER},
		{"anicut", 'n', "<FLOAT>", 0, "Skip reports with ANI below this value. [0.95]", ANI_GROUP_FILTER},
		{"control", 'c', "<FLOAT>", 0, "Skip duplicated samples with distance below this value. [0]", ANI_GROUP_FILTER},
		{"ctxcut", 't', "<INT>", 0, "Skip reports with overlapped context count below this value. [3]", ANI_GROUP_FILTER},
		{"slmetrics", 's', "<+-1..8>", 0, "Metric: Best(1), Recalibrated(2), CtxMoE(3), Naive(4), MashD(5), AafD(6), MashD_if_far(7), AafD_if_far(8). For unassembled/qraw, 1..4 use Naive. In matrix mode positive reports distance and negative reports ANI; detail prints both. [1]", ANI_GROUP_FILTER},
		{"unified-metric", ANI_UNIFIED_METRIC, 0, 0, "Expert: in unassembled/qraw mode, honor -s instead of forcing 1..4 to Naive; Best/Recalibrated still fall back when unavailable.", ANI_GROUP_FILTER},

		{0, 0, 0, 0, "Reporting and output:", ANI_GROUP_REPORT},
		{"diagonal", 'd', 0, 0, "Set diagonal values.", ANI_GROUP_REPORT},
		{"exception", 'e', "<INT>", 0, "Distance value to use when skipped. [1]", ANI_GROUP_REPORT},
		{"glist", 'g', "<FILE>", 0, "Sample output file for KSSD set grouping.", ANI_GROUP_REPORT},
		{"outfmt", 'm', "<0/1/2>", 0, "Output format: detail(0), matrix(1), or triangle(2). [0]", ANI_GROUP_REPORT},
		{"outfile", 'o', "<FILE>", 0, "Output file. [STDOUT]", ANI_GROUP_REPORT},
		{"raw-output", ANI_RAW_OUTPUT, 0, 0, "Skip calibrated/best ANI computation; selected calibrated metrics fall back to raw distances when unavailable.", ANI_GROUP_REPORT},
		{"top", 'N', "<INT>", 0, "Report at most top N references per query. [all]", ANI_GROUP_REPORT},

		{0, 0, 0, 0, "Execution:", ANI_GROUP_EXECUTION},
		{"threads", 'p', "<INT>", 0, "Number of threads to use.", ANI_GROUP_EXECUTION},

		{0, 0, 0, 0, "Direct FASTA/FASTQ/sketch input mode:", ANI_GROUP_AUTO_MODE},
		{"pair", 777, 0, 0, "Treat positional inputs as reference then query inputs; optional when -r/-q/--qraw are omitted.", ANI_GROUP_AUTO_MODE},
		{"pipecmd", ANI_SKETCH_PIPECMD, "<cmd>", 0, "Stream each raw input through command; '{}' is replaced by the input path, otherwise the path is appended.", ANI_GROUP_AUTO_MODE},

		{0, 0, 0, 0, "Auto-sketch options for FASTA/FASTQ inputs:", ANI_GROUP_AUTOSKETCH},
		{"use_coden_ctxobj", ANI_SKETCH_USE_CODEN, 0, 0, "Use the coden context-object structure pattern. [default]", ANI_GROUP_AUTOSKETCH},
		{"DimRdcFold", ANI_SKETCH_DRFOLD, "<INT>", 0, "K-mer space downsampling rate 1/2^f. [8]", ANI_GROUP_AUTOSKETCH},
		{"ctxlen", ANI_SKETCH_CTXLEN, "<INT>", 0, "Manual half context length; disables default coden mode.", ANI_GROUP_AUTOSKETCH},
		{"outerobjlen", ANI_SKETCH_OUTEROBJLEN, "<INT>", 0, "Manual half outer object length; disables default coden mode.", ANI_GROUP_AUTOSKETCH},
		{"innerobjlen", ANI_SKETCH_INNEROBJLEN, "<INT>", 0, "Manual inner object length; disables default coden mode.", ANI_GROUP_AUTOSKETCH},
		{"LstKmerOcrs", ANI_SKETCH_LST_KMER_OCRS, "<INT>", 0, "Least k-mer occurrence required in FASTQ input. [1]", ANI_GROUP_AUTOSKETCH},
		{"npercentile", ANI_SKETCH_NPERCENTILE, "<0..1>", 0, "Weighted lower-tail count percentile among k-mers passing --LstKmerOcrs.", ANI_GROUP_AUTOSKETCH},
		{"ncap", ANI_SKETCH_NCAP, "<INT>", OPTION_HIDDEN, "Cap for effective lower k-mer occurrence cutoff; 0 disables. [0]", ANI_GROUP_AUTOSKETCH},
		{"readsQC", ANI_SKETCH_READS_QC, 0, 0, "Infer and apply read-count QC filtering.", ANI_GROUP_AUTOSKETCH},
		{"abundance", ANI_SKETCH_ABUNDANCE, 0, 0, "Use abundance-estimation mode.", ANI_GROUP_AUTOSKETCH},
		{"anno", ANI_SKETCH_ANNO, 0, 0, "Write FASTA/FASTQ header annotations to lcofiles.anno.", ANI_GROUP_AUTOSKETCH},
		{"conflict", ANI_SKETCH_CONFLICT, 0, 0, "Keep conflicting context-objects.", ANI_GROUP_AUTOSKETCH},
		{"asone", ANI_SKETCH_ASONE, 0, 0, "Treat each raw input path as one final genome.", ANI_GROUP_AUTOSKETCH},
		{"splitmfa", ANI_SKETCH_SPLITMFA, 0, 0, "Treat each multi-FASTA record as one genome.", ANI_GROUP_AUTOSKETCH},
		{0}};

static char doc_ani[] =
	"\n"
	"Estimate average nucleotide identity (ANI)."
	"\v"
	"Use -r/-q with sketches or FASTA/FASTQ(.gz) files, -r/--qraw for raw-read query sketches or FASTA/FASTQ(.gz) files, or positional FASTA/FASTQ/sketch inputs.\n"
	"For --qraw, query sequence files are sketched with conflicting objects kept; reference sketches may be conflict-free or conflict-containing.\n"
	"Use --ignoreconflict to skip reference-side conflict contexts instead of comparing them.\n"
	"When -r, -q, and --qraw are omitted, the first positional input is the reference and the rest are queries.\n"
	"Use --reflist/--qrylist to read reference and query FASTA/FASTQ/sketch paths from files.\n"
	"If any auto ANI input is a sketch, FASTA/FASTQ inputs are sketched with the first sketch's lcofiles.stat parameters.\n"
	"If no auto ANI input is a sketch, FASTA/FASTQ inputs use sketch defaults: coden mode and --DimRdcFold 8.\n"
	"Auto-sketching accepts sketch long options such as --DimRdcFold, --ctxlen, --readsQC, and --conflict.\n"
	"For sequence files, -q sketches then runs normal sketch ANI; --qraw sketches then runs raw-query ANI.\n"
	"Direct -q sequence inputs keep conflicts only when --conflict is set; --qraw sequence inputs keep query conflicts by default.\n"
	"Detail output uses the selected metric schema: Qry, Ref, ANI, Distance, Confidence, Selected_metric, diagnostics, Ref_annotation.\n"
	"With one sketch and -m1 or -m2, ANI computes a self full matrix or lower triangle.\n"
	"Use --raw-output to skip calibrated/best ANI computation; selected calibrated fields fall back to raw distances when unavailable.\n"
	"Reference inverted indexes are auto-created for detail ANI when ref genome count is > KSSD3A_AUTO_REF_INDEX_THRESHOLD [100]; set it to 0 to disable.\n"
	"Use '-' as one raw input to read FASTA/FASTQ from stdin, or --pipecmd CMD to stream each raw input through a command.\n"
	"\n"
	"Examples:\n"
	"  kssd3a ani -r ref_sketches -q qry_sketches -m 0 > ani.tsv\n"
	"  kssd3a ani -r ref_sketches -q query.fasta.gz -m 0 > ani.tsv\n"
	"  kssd3a ani -r ref_sketches --qraw raw_read_qry_sketches -m 0 > raw_ani.tsv\n"
	"  kssd3a ani -r ref_sketches --qraw reads.fastq.gz -m 0 > raw_ani.tsv\n"
	"  kssd3a ani -q sketches -m2 -s -1 -d > self_ani_triangle.tsv\n"
	"  kssd3a ani -f0 -n0 -o pair.tsv ref.fasta query.fasta\n"
	"  samtools fastq reads.bam | kssd3a ani --conflict \\\n"
	"    ref_read_sketch - -o raw_ani.tsv\n"
	"  kssd3a ani --conflict --pipecmd 'samtools fastq {}' \\\n"
	"    ref_read_sketch reads.bam -o raw_ani.tsv\n"
	"  kssd3a ani --DimRdcFold 8 --reflist refs.txt --qrylist qrys.txt -o ani.tsv\n"
	"  kssd3a ani ref_sketch query.fasta query_sketch";

ani_opt_t ani_opt = {
	.fmt = 0, // 0:detail, 1:matrix 2: triangle
	.c = 0.0, // control duplicated sample by skip distance < c;
	.p = 1,
	.d = 0,	   // diagonal
	.v = 0,	   // naive model
	.s = 1,	   // select metrics: 1:standard, 2:MashD, 3:AafD, 4:MashD_if_far, 5:AafD_if_far
	.pair = 0, // pairwise compute
	.unassembled = 0,
	.unified_metric = 0,
	.ignoreconflict = 0,
	.raw_output = 0,
	.ctxcut = 3,
	.afcut = 0.5,
	.afcut_set = false,
	.anicut = 0.95,
	.sketch_hclen = 11,
	.sketch_holen = 0,
	.sketch_iolen = 0,
	.sketch_drfold = 8,
	.sketch_kmerocrs = 1,
	.sketch_ncap = 0,
	.sketch_npercentile = 0.0,
	.sketch_reads_qc = false,
	.sketch_abundance = false,
	.sketch_asone = false,
	.sketch_conflict = false,
	.sketch_anno = false,
	.sketch_split_mfa = false,
	.sketch_coden_ctxobj_pattern = true,
	.ntop = -1, //all 
	.e = 1, // abort
	.index[0] = '\0',
	.refdir[0] = '\0',
	.qrydir[0] = '\0',
	.reflist[0] = '\0',
	.qrylist[0] = '\0',
	.sketch_pipecmd[0] = '\0',
	.outf[0] = '\0',
	.gl[0] = '\0',
	//	.model[0] = '\0',
	.num_remaining_args = 0, // int num_remaining_args; no option arguments num.
	.remaining_args = NULL	 // char **remaining_args; no option arguments array.
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
	if (errno != 0 || end == arg || *end != '\0' || !isfinite(value) || value < min_value || value > max_value)
		argp_error(state, "%s requires a number in range %.6g..%.6g", option_name, min_value, max_value);
	return value;
}

static void copy_path_arg(struct argp_state *state, const char *option_name,
						  char *dest, size_t dest_size, const char *arg)
{
	if (strlen(arg) >= dest_size)
		argp_error(state, "%s path is too long; maximum supported length is %zu bytes", option_name, dest_size - 1);
	snprintf(dest, dest_size, "%s", arg);
}

static error_t parse_ani(int key, char *arg, struct argp_state *state)
{
	struct arg_ani *ani = state->input;

	switch (key)
	{
	case 'm':
	{
		ani_opt.fmt = parse_int_range(state, "-m/--outfmt", arg, 0, 2);
		break;
	}
	case 'c':
	{
		ani_opt.c = parse_nonnegative_double(state, "-c/--control", arg);
		break;
	}
	case 'e':
	{
		ani_opt.e = parse_int_range(state, "-e/--exception", arg, 1, INT_MAX);
		break;
	}
	case 's':
	{
		ani_opt.s = parse_int_range(state, "-s/--slmetrics", arg, -8, 8);
		if (ani_opt.s == 0)
			argp_error(state, "-s/--slmetrics accepts -8..-1 or 1..8, but not 0");
		break;
	}
	case 'p':
	{
		ani_opt.p = parse_int_range(state, "-p/--threads", arg, 1, 65536);
		break;
	}
	case ANI_SKETCH_USE_CODEN:
	{
		ani->coden_pattern_seen = 1;
		ani_opt.sketch_coden_ctxobj_pattern = true;
		break;
	}
	case ANI_SKETCH_CTXLEN:
	{
		ani->manual_pattern_seen = 1;
		ani_opt.sketch_coden_ctxobj_pattern = false;
		ani_opt.sketch_hclen = parse_int_range(state, "--ctxlen", arg, 4, 16);
		break;
	}
	case ANI_SKETCH_OUTEROBJLEN:
	{
		ani->manual_pattern_seen = 1;
		ani_opt.sketch_coden_ctxobj_pattern = false;
		ani_opt.sketch_holen = parse_int_range(state, "--outerobjlen", arg, 0, 8);
		break;
	}
	case ANI_SKETCH_INNEROBJLEN:
	{
		ani->manual_pattern_seen = 1;
		ani_opt.sketch_coden_ctxobj_pattern = false;
		ani_opt.sketch_iolen = parse_int_range(state, "--innerobjlen", arg, 0, 8);
		break;
	}
	case ANI_SKETCH_DRFOLD:
	{
		ani_opt.sketch_drfold = parse_int_range(state, "--DimRdcFold", arg, 0, 24);
		break;
	}
	case ANI_SKETCH_LST_KMER_OCRS:
	{
		ani_opt.sketch_kmerocrs = parse_int_range(state, "--LstKmerOcrs", arg, 1, 65536);
		break;
	}
	case ANI_SKETCH_NPERCENTILE:
	{
		ani_opt.sketch_npercentile = parse_double_range(state, "--npercentile", arg, 0.0, 1.0);
		break;
	}
	case ANI_SKETCH_NCAP:
	{
		ani_opt.sketch_ncap = parse_int_range(state, "--ncap", arg, 0, 65536);
		break;
	}
	case ANI_SKETCH_READS_QC:
	{
		ani_opt.sketch_reads_qc = true;
		if (ani_opt.sketch_kmerocrs < 2)
			ani_opt.sketch_kmerocrs = 2;
		break;
	}
	case ANI_SKETCH_ABUNDANCE:
	{
		ani_opt.sketch_abundance = true;
		break;
	}
	case ANI_SKETCH_ANNO:
	{
		ani_opt.sketch_anno = true;
		break;
	}
	case ANI_SKETCH_CONFLICT:
	{
		ani_opt.sketch_conflict = true;
		break;
	}
	case ANI_SKETCH_ASONE:
	{
		ani_opt.sketch_asone = true;
		break;
	}
	case ANI_SKETCH_SPLITMFA:
	{
		ani_opt.sketch_split_mfa = true;
		break;
	}
	case ANI_SKETCH_PIPECMD:
	{
		copy_path_arg(state, "--pipecmd", ani_opt.sketch_pipecmd, sizeof(ani_opt.sketch_pipecmd), arg);
		break;
	}
	case ANI_RAW_OUTPUT:
	{
		ani_opt.raw_output = true;
		break;
	}
	case ANI_UNIFIED_METRIC:
	{
		ani_opt.unified_metric = true;
		break;
	}
	case 't':
	{
		ani_opt.ctxcut = parse_int_range(state, "-t/--ctxcut", arg, 0, INT_MAX);
		break;
	}
	case 'f':
	{
		ani_opt.afcut = (float)parse_nonnegative_double(state, "-f/--afcut", arg);
		ani_opt.afcut_set = true;
		break;
	}
	case 'n':
	{
		ani_opt.anicut = (float)parse_nonnegative_double(state, "-n/--anicut", arg);
		break;
	}
	case 'N':
	{
		ani_opt.ntop = parse_int_range(state, "-N/--top", arg, 1, INT_MAX);
		break;
	}
	case 'i':
	{
		copy_path_arg(state, "-i/--index", ani_opt.index, sizeof(ani_opt.index), arg);
		break;
	}
	case 'M':
	{
		copy_path_arg(state, "-M/--model", ani_opt.model, sizeof(ani_opt.model), arg);
		break;
	}
	case 'q':
	{
		copy_path_arg(state, "-q/--query", ani_opt.qrydir, sizeof(ani_opt.qrydir), arg);
		break;
	}
	case 780:
	{
		copy_path_arg(state, "--qrylist", ani_opt.qrylist, sizeof(ani_opt.qrylist), arg);
		break;
	}
	case 'r':
	{
		copy_path_arg(state, "-r/--ref", ani_opt.refdir, sizeof(ani_opt.refdir), arg);
		break;
	}
	case 779:
	{
		copy_path_arg(state, "--reflist", ani_opt.reflist, sizeof(ani_opt.reflist), arg);
		break;
	}
	case 'o':
	{
		copy_path_arg(state, "-o/--outfile", ani_opt.outf, sizeof(ani_opt.outf), arg);
		break;
	}
	case 'g':
	{
		copy_path_arg(state, "-g/--glist", ani_opt.gl, sizeof(ani_opt.gl), arg);
		break;
	}
	case 'd':
	{
		ani_opt.d = 1;
		break;
	}
	case 'v':
	{
		ani_opt.v = 1;
		break;
	}
	case 778:
	{
		ani_opt.ignoreconflict = 1;
		break;
	}
	case 'x':
	{
		copy_path_arg(state, "--qraw", ani_opt.qrydir, sizeof(ani_opt.qrydir), arg);
		ani_opt.unassembled = 1; 
		ani_opt.v = 1;
		break;
	}
	case 777:
	{
		ani_opt.pair = 1;
		break;
	}
	case ARGP_KEY_ARGS:
	{
		ani_opt.num_remaining_args = state->argc - state->next;
		ani_opt.remaining_args = state->argv + state->next;
		break;
	}
	case ARGP_KEY_END:
	{
		const bool has_reflist = ani_opt.reflist[0] != '\0';
		const bool has_qrylist = ani_opt.qrylist[0] != '\0';

		if (ani->coden_pattern_seen && ani->manual_pattern_seen)
			argp_error(state, "--use_coden_ctxobj selects the coden pattern; do not combine it with --ctxlen, --outerobjlen, or --innerobjlen.");
		if (!ani_opt.sketch_coden_ctxobj_pattern)
		{
			int sketch_klen = ani_opt.sketch_iolen + 2 * (ani_opt.sketch_holen + ani_opt.sketch_hclen);
			if (sketch_klen > 32)
				argp_error(state, "auto-sketch k-mer length %d should be at most 32", sketch_klen);
		}
		if (ani_opt.sketch_pipecmd[0] != '\0' && ani_opt.sketch_split_mfa)
			argp_error(state, "--splitmfa does not support --pipecmd streaming inputs");

		if (has_reflist || has_qrylist)
		{
			if (!has_reflist || !has_qrylist)
			{
				printf("\nError: --reflist and --qrylist must be used together.\n\n");
				argp_state_help(state, stdout, ARGP_HELP_STD_HELP);
				argp_usage(state);
			}
			if (ani_opt.refdir[0] != '\0' || ani_opt.qrydir[0] != '\0')
			{
				printf("\nError: --reflist/--qrylist use auto ANI inputs; do not combine them with -r, -q, or --qraw.\n\n");
				argp_state_help(state, stdout, ARGP_HELP_STD_HELP);
				argp_usage(state);
			}
			if (ani_opt.pair || ani_opt.num_remaining_args > 0)
			{
				printf("\nError: --reflist/--qrylist cannot be combined with --pair or positional inputs.\n\n");
				argp_state_help(state, stdout, ARGP_HELP_STD_HELP);
				argp_usage(state);
			}
		}
		else if (ani_opt.pair)
		{
			if (ani_opt.refdir[0] != '\0' || ani_opt.qrydir[0] != '\0')
			{
				printf("\nError: --pair uses positional FASTA/FASTQ/sketch inputs; do not combine it with -r, -q, or --qraw.\n\n");
				argp_state_help(state, stdout, ARGP_HELP_STD_HELP);
				argp_usage(state);
			}
			if (ani_opt.num_remaining_args < 2)
			{
				printf("\nError: --pair requires one reference input followed by at least one query input.\n\n");
				argp_state_help(state, stdout, ARGP_HELP_STD_HELP);
				argp_usage(state);
			}
		}
			else if (ani_opt.refdir[0] == '\0' && ani_opt.qrydir[0] == '\0' &&
					 ani_opt.num_remaining_args >= 2)
			{
				ani_opt.pair = 1;
			}
			else if (ani_opt.refdir[0] == '\0' && ani_opt.qrydir[0] == '\0' &&
					 ani_opt.num_remaining_args == 1 && ani_opt.fmt != 0)
			{
				copy_path_arg(state, "sketch", ani_opt.qrydir, sizeof(ani_opt.qrydir),
							  ani_opt.remaining_args[0]);
			}
			else if (ani_opt.qrydir[0] == '\0')
			{
				printf("\nError: Mandatory options: '-q' or '--qraw' are missing unless positional inputs are used.\n\n");
				argp_state_help(state, stdout, ARGP_HELP_STD_HELP);
				argp_usage(state);
			}
			else if (ani_opt.refdir[0] == '\0')
			{
				if (ani_opt.fmt == 0)
				{
					printf("\nError: -r/--ref is required with -q/--query or --qraw for detail output. Use -m1 or -m2 for one-sketch self matrices.\n\n");
					argp_state_help(state, stdout, ARGP_HELP_STD_HELP);
					argp_usage(state);
				}
			}
		if (ani_opt.s < -8 || ani_opt.s > 8 || ani_opt.s == 0)
			argp_error(state, "-s option should be within range 1..8 or -8..-1");

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

static struct argp argp_ani =
	{
		opt_ani,
		parse_ani,
			"[-r REF_SKETCH -q QRY | --qraw RAW_QRY | --reflist REF_LIST --qrylist QRY_LIST | [--pair] REF_INPUT QUERY_INPUT... | SELF_SKETCH with -m1/-m2]",
			doc_ani};
extern const char sorted_comb_ctxgid64obj32[];
extern void gen_inverted_index4comblco(const char *sketchdir);

get_generic_dist_from_features_fn get_generic_dist_from_features = NULL;
int ani_model_drfold = ANI_MODEL_REFERENCE_DRFOLD;

static int qraw_multi_query_index_threshold(void)
{
	const char *env = getenv("KSSD3A_QRAW_MULTI_THRESHOLD");
	if (!env || env[0] == '\0')
		return 20;

	char *end = NULL;
	errno = 0;
	long value = strtol(env, &end, 10);
	if (errno != 0 || end == env || *end != '\0' || value < 1 || value > INT_MAX)
		errx(EXIT_FAILURE, "invalid KSSD3A_QRAW_MULTI_THRESHOLD='%s' (use positive integer)", env);
	return (int)value;
}

static int auto_ref_index_threshold(void)
{
	const char *env = getenv("KSSD3A_AUTO_REF_INDEX_THRESHOLD");
	if (!env || env[0] == '\0')
		return 100;

	char *end = NULL;
	errno = 0;
	long value = strtol(env, &end, 10);
	if (errno != 0 || end == env || *end != '\0' || value < 0 || value > INT_MAX)
		errx(EXIT_FAILURE, "invalid KSSD3A_AUTO_REF_INDEX_THRESHOLD='%s' (use non-negative integer; 0 disables)", env);
	return (int)value;
}

static int ani_matrix_direct_threshold(void)
{
	const char *name = "KSSD3A_ANI_MATRIX_DIRECT_THRESHOLD";
	const char *env = getenv(name);
	if (!env || env[0] == '\0') {
		name = "KSSD3A_ANI_M1_DIRECT_THRESHOLD";
		env = getenv(name);
	}
	if (!env || env[0] == '\0')
		return 1000;

	char *end = NULL;
	errno = 0;
	long value = strtol(env, &end, 10);
	if (errno != 0 || end == env || *end != '\0' || value < 0 || value > INT_MAX)
		errx(EXIT_FAILURE, "invalid %s='%s' (use non-negative integer)", name, env);
	return (int)value;
}

static bool ani_matrix_should_use_index(int ref_infile_num, int qry_infile_num, int direct_threshold)
{
	if (ref_infile_num < 1 || qry_infile_num < 1)
		return false;
	if (direct_threshold == 0)
		return true;
	return ref_infile_num > direct_threshold || qry_infile_num > direct_threshold;
}

static void ensure_sorted_ref_index_for_ani_m1(const char *refdir)
{
	if (file_exists_in_folder(refdir, sorted_comb_ctxgid64obj32))
		return;
	gen_inverted_index4comblco(refdir);
	if (!file_exists_in_folder(refdir, sorted_comb_ctxgid64obj32))
		errx(EXIT_FAILURE, "failed to create reference index '%s/%s'",
			 refdir, sorted_comb_ctxgid64obj32);
}

static int sketch_infile_num(const char *sketch_dir)
{
	size_t stat_size = 0;
	char *stat_path = test_get_fullpath(sketch_dir, sketch_stat);
	dim_sketch_stat_t *stat = read_from_file(stat_path, &stat_size);
	free(stat_path);
	const int infile_num = stat->infile_num;
	free_read_from_file(stat, stat_size);
	return infile_num;
}

static bool force_ref_index_requested(void)
{
	const char *env = getenv("KSSD3A_FORCE_REF_INDEX");
	return env && env[0] != '\0' && strcmp(env, "0") != 0;
}

typedef struct
{
	char **items;
	int n;
	int cap;
} path_vec_t;

static void path_vec_push(path_vec_t *vec, const char *path)
{
	if (vec->n == vec->cap)
	{
		int new_cap = vec->cap ? vec->cap * 2 : 8;
		char **new_items = realloc(vec->items, (size_t)new_cap * sizeof(new_items[0]));
		if (!new_items)
			err(errno, "%s(): OOM path vector", __func__);
		vec->items = new_items;
		vec->cap = new_cap;
	}
	vec->items[vec->n] = strdup(path);
	if (!vec->items[vec->n])
		err(errno, "%s(): OOM path copy", __func__);
	vec->n++;
}

static void path_vec_free(path_vec_t *vec)
{
	for (int i = 0; i < vec->n; ++i)
		free(vec->items[i]);
	free(vec->items);
	vec->items = NULL;
	vec->n = 0;
	vec->cap = 0;
}

static char *trim_path_line(char *line)
{
	char *start = line;
	while (isspace((unsigned char)*start))
		start++;

	char *end = start + strlen(start);
	while (end > start && isspace((unsigned char)end[-1]))
		end--;
	*end = '\0';
	return start;
}

static void read_path_list_to_vec(const char *list_path, path_vec_t *out, const char *option_name)
{
	FILE *list = fopen(list_path, "r");
	if (!list)
		err(errno, "can't open %s file %s", option_name, list_path);

	char *line = malloc(LMAX);
	if (!line)
		err(errno, "%s(): OOM list line buffer", __func__);

	int line_no = 0;
	while (fgets(line, LMAX, list) != NULL)
	{
		line_no++;
		size_t len = strlen(line);
		if (len == LMAX - 1 && line[len - 1] != '\n' && line[len - 1] != '\r')
			errx(EXIT_FAILURE, "%s file %s line %d exceeds maximum line length %d",
				 option_name, list_path, line_no, LMAX - 1);

		char *path = trim_path_line(line);
		if (path[0] == '\0')
			continue;
		if (strlen(path) >= PATHLEN)
			errx(EXIT_FAILURE, "%s file %s line %d path exceeds maximum length %d",
				 option_name, list_path, line_no, PATHLEN - 1);
		path_vec_push(out, path);
	}
	if (ferror(list))
		err(errno, "failed while reading %s file %s", option_name, list_path);

	free(line);
	fclose(list);

	if (out->n < 1)
		errx(EXIT_FAILURE, "%s file %s has no input paths", option_name, list_path);
}

static bool is_lco_sketch_dir(const char *path)
{
	struct stat st;
	if (stat(path, &st) != 0 || !S_ISDIR(st.st_mode))
		return false;
	return file_exists_in_folder(path, sketch_stat) &&
		   file_exists_in_folder(path, combined_sketch_suffix) &&
		   file_exists_in_folder(path, idx_sketch_suffix);
}

static dim_sketch_stat_t read_lco_sketch_stat(const char *sketch_dir)
{
	size_t stat_size = 0;
	char *stat_path = test_get_fullpath(sketch_dir, sketch_stat);
	dim_sketch_stat_t *stat = read_from_file(stat_path, &stat_size);
	free(stat_path);
	if (stat_size < sizeof(*stat))
		err(EINVAL, "%s(): %s/%s has %zu bytes, expected at least %zu",
			__func__, sketch_dir, sketch_stat, stat_size, sizeof(*stat));
	dim_sketch_stat_t copy = *stat;
	free_read_from_file(stat, stat_size);
	return copy;
}

static bool same_ani_auto_sketch_stat(const dim_sketch_stat_t *a, const dim_sketch_stat_t *b)
{
	return a->hash_id == b->hash_id &&
		   a->koc == b->koc &&
		   a->conflict == b->conflict &&
		   a->coden_len == b->coden_len &&
		   a->klen == b->klen &&
		   a->hclen == b->hclen &&
		   a->holen == b->holen &&
		   a->drfold == b->drfold;
}

static void describe_sketch_stat_mismatch(
	const char *path,
	const dim_sketch_stat_t *expected,
	const dim_sketch_stat_t *actual)
{
	errx(EXIT_FAILURE,
		 "%s sketch parameters do not match first auto ANI sketch "
		 "(hash_id %u/%u, koc %d/%d, conflict %d/%d, coden_len %d/%d, "
		 "klen %d/%d, hclen %d/%d, holen %d/%d, drfold %d/%d)",
		 path,
		 actual->hash_id, expected->hash_id,
		 actual->koc, expected->koc,
		 actual->conflict, expected->conflict,
		 actual->coden_len, expected->coden_len,
		 actual->klen, expected->klen,
		 actual->hclen, expected->hclen,
		 actual->holen, expected->holen,
		 actual->drfold, expected->drfold);
}

extern uint32_t FILTER;
extern uint32_t hash_id;
extern dim_sketch_stat_t comblco_stat_one;
extern uint32_t get_sketching_id(uint32_t hclen, uint32_t holen, uint32_t iolen, uint32_t drfold, uint32_t FILTER);
extern void compute_sketch(sketch_opt_t *sketch_opt_val, infile_tab_t *infile_stat);
extern int merge_comblco(sketch_opt_t *sketch_opt_val);

static dim_sketch_stat_t default_auto_ani_sketch_stat(const ani_opt_t *opt)
{
	dim_sketch_stat_t stat = {0};
	stat.koc = opt->sketch_abundance;
	stat.conflict = opt->sketch_conflict;
	stat.drfold = opt->sketch_drfold;
	if (opt->sketch_coden_ctxobj_pattern)
	{
		stat.coden_len = NUM_CODENS;
#if NUM_CODENS < 11
		stat.klen = 3 * NUM_CODENS + 1;
#else
		stat.klen = 32;
#endif
		stat.hclen = 0;
		stat.holen = 0;
	}
	else
	{
		stat.coden_len = 0;
		stat.hclen = opt->sketch_hclen;
		stat.holen = opt->sketch_holen;
		stat.klen = opt->sketch_iolen + 2 * (opt->sketch_hclen + opt->sketch_holen);
	}
	const uint32_t filter = UINT32_MAX >> stat.drfold;
	if (stat.coden_len > 0)
		stat.hash_id = get_sketching_id((uint32_t)stat.coden_len, 0, 0, (uint32_t)stat.drfold, filter);
	else
		stat.hash_id = get_sketching_id((uint32_t)stat.hclen, (uint32_t)stat.holen,
										(uint32_t)opt->sketch_iolen, (uint32_t)stat.drfold, filter);
	stat.infile_num = 0;
	return stat;
}

static int stat_inner_object_len(const dim_sketch_stat_t *stat)
{
	return stat->klen - 2 * (stat->hclen + stat->holen);
}

static void validate_supported_auto_sketch_stat(const dim_sketch_stat_t *stat)
{
	if (stat->drfold < 0 || stat->drfold > 32)
		errx(EXIT_FAILURE, "unsupported sketch drfold %d in first auto ANI sketch", stat->drfold);
	if (stat->klen > 32)
		errx(EXIT_FAILURE, "first auto ANI sketch has unsupported klen=%d; expected at most 32", stat->klen);
	if (stat->coden_len > 0)
	{
		if (stat->coden_len != NUM_CODENS)
			errx(EXIT_FAILURE,
				 "first auto ANI sketch uses coden_len=%d, but this binary was built with NUM_CODENS=%d",
				 stat->coden_len, NUM_CODENS);
#if NUM_CODENS < 11
		const int expected_klen = 3 * NUM_CODENS + 1;
#else
		const int expected_klen = 32;
#endif
		if (stat->klen != expected_klen)
			errx(EXIT_FAILURE,
				 "first auto ANI sketch klen=%d is not compatible with coden_len=%d in this binary",
				 stat->klen, stat->coden_len);
	}
	else if (stat_inner_object_len(stat) < 0)
	{
		errx(EXIT_FAILURE,
			 "first auto ANI sketch has invalid manual context-object lengths: klen=%d hclen=%d holen=%d",
			 stat->klen, stat->hclen, stat->holen);
	}
}

static sketch_opt_t sketch_opt_from_stat(const dim_sketch_stat_t *stat, const char *outdir, const ani_opt_t *ani_opt_val)
{
	sketch_opt_t opt = {
		.hclen = stat->hclen,
		.holen = stat->holen,
		.iolen = stat_inner_object_len(stat),
		.drfold = stat->drfold,
		.kmerocrs = ani_opt_val->sketch_kmerocrs,
		.npercentile = ani_opt_val->sketch_npercentile,
		.ncap = ani_opt_val->sketch_ncap,
		.reads_qc = ani_opt_val->sketch_reads_qc,
		.p = ani_opt_val->p > 0 ? ani_opt_val->p : 1,
		.abundance = stat->koc,
		.asone = ani_opt_val->sketch_asone,
		.conflict = stat->conflict,
		.anno = ani_opt_val->sketch_anno,
		.compute_meta = true,
		.merge_comblco = false,
		.sketch_qc = false,
		.split_mfa = ani_opt_val->sketch_split_mfa,
		.coden_ctxobj_pattern = stat->coden_len > 0,
		.index[0] = '\0',
		.fpath = NULL,
		.outdir = (char *)outdir,
		.pipecmd = ani_opt_val->sketch_pipecmd[0] ? (char *)ani_opt_val->sketch_pipecmd : NULL,
		.num_remaining_args = 0,
		.remaining_args = NULL,
	};
	if (opt.coden_ctxobj_pattern)
		opt.iolen = 0;
	return opt;
}

static void apply_auto_sketch_stat(const dim_sketch_stat_t *stat, int infile_num)
{
	FILTER = UINT32_MAX >> stat->drfold;
	hash_id = stat->hash_id;
	comblco_stat_one = *stat;
	comblco_stat_one.infile_num = infile_num;
	const_comask_init(&comblco_stat_one);
	set_uint64kmer2generic_ctxobj(stat->coden_len > 0);
}

static void free_infile_tab(infile_tab_t *tab)
{
	if (!tab)
		return;
	for (int i = 0; i < tab->infile_num; ++i)
		free(tab->organized_infile_tab[i].fpath);
	free(tab->organized_infile_tab);
	free(tab);
}

static void remove_tree(const char *path)
{
	DIR *dir = opendir(path);
	if (!dir)
		return;
	struct dirent *ent;
	while ((ent = readdir(dir)) != NULL)
	{
		if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0)
			continue;
		char child[PATHLEN * 2];
		snprintf(child, sizeof(child), "%s/%s", path, ent->d_name);
		struct stat st;
		if (lstat(child, &st) != 0)
			continue;
		if (S_ISDIR(st.st_mode))
			remove_tree(child);
		else
			unlink(child);
	}
	closedir(dir);
	rmdir(path);
}

static void sketch_sequence_arg_to_dir(
	const char *input,
	const dim_sketch_stat_t *template_stat,
	const char *outdir,
	const ani_opt_t *ani_opt_val)
{
	char *arg = (char *)input;
	int fmt_ck = (ani_opt_val->sketch_pipecmd[0] == '\0' && strcmp(input, "-") != 0) ? 1 : 0;
	infile_tab_t *tab = organize_infile_frm_arg(1, &arg, fmt_ck);
	if (!tab || tab->infile_num < 1)
		errx(EXIT_FAILURE, "%s is not a valid FASTA/FASTQ input for auto ANI", input);

	int infile_num = ani_opt_val->sketch_asone ? (tab->infile_num < 1 ? 0 : 1) : tab->infile_num;
	apply_auto_sketch_stat(template_stat, infile_num);
	sketch_opt_t opt = sketch_opt_from_stat(template_stat, outdir, ani_opt_val);
	mkdir_p(outdir);
	compute_sketch(&opt, tab);
	free_infile_tab(tab);
}

static char *merge_auto_sketches(path_vec_t *sketch_dirs, const char *outdir)
{
	if (sketch_dirs->n == 1)
	{
		char *copy = strdup(sketch_dirs->items[0]);
		if (!copy)
			err(errno, "%s(): OOM sketch path copy", __func__);
		return copy;
	}

	sketch_opt_t merge_opt = {
		.outdir = (char *)outdir,
		.num_remaining_args = sketch_dirs->n,
		.remaining_args = sketch_dirs->items,
	};
	mkdir_p(outdir);
	merge_comblco(&merge_opt);
	char *copy = strdup(outdir);
	if (!copy)
		err(errno, "%s(): OOM merged sketch path copy", __func__);
	return copy;
}

static bool paths_refer_to_same_file(const char *a, const char *b)
{
	struct stat sa;
	struct stat sb;
	return a && b &&
		   stat(a, &sa) == 0 &&
		   stat(b, &sb) == 0 &&
		   sa.st_dev == sb.st_dev &&
		   sa.st_ino == sb.st_ino;
}

static int run_sketch_ani(ani_opt_t *opt)
{
	if (opt->unassembled && !opt->afcut_set)
		opt->afcut = 0.2f;

	size_t file_size;
	char *qry_stat_path = test_get_fullpath(opt->qrydir, sketch_stat);
	dim_sketch_stat_t *qry_dim_sketch_stat = read_from_file(qry_stat_path, &file_size);
	free(qry_stat_path);
	const int qry_infile_num = qry_dim_sketch_stat->infile_num;

	ani_model_drfold = qry_dim_sketch_stat->drfold;
	const_comask_init(qry_dim_sketch_stat);
	const bool force_ref_index = force_ref_index_requested();
	const int matrix_direct_n = opt->fmt ? ani_matrix_direct_threshold() : 0;
	if (opt->refdir[0] == '\0')
	{
		if (opt->fmt == 0)
			errx(EXIT_FAILURE, "one-sketch ANI detail output is not supported; use -m1 full matrix or -m2 triangle");
		if (opt->fmt == 1 &&
			(force_ref_index || ani_matrix_should_use_index(qry_infile_num, qry_infile_num, matrix_direct_n)))
		{
			if (force_ref_index && !file_exists_in_folder(opt->qrydir, sorted_comb_ctxgid64obj32))
				errx(EXIT_FAILURE, "KSSD3A_FORCE_REF_INDEX requires reference index '%s/%s'",
					 opt->qrydir, sorted_comb_ctxgid64obj32);
			else
				ensure_sorted_ref_index_for_ani_m1(opt->qrydir);
			free_read_from_file(qry_dim_sketch_stat, file_size);
			if (comb_sortedsketch64_indexed_self_full(opt))
				return 1;
			comb_sortedsketch64_self_matrix(opt);
			return 1;
		}
		if (opt->fmt == 2 &&
			(force_ref_index || ani_matrix_should_use_index(qry_infile_num, qry_infile_num, matrix_direct_n)))
		{
			free_read_from_file(qry_dim_sketch_stat, file_size);
			if (comb_sortedsketch64_indexed_self_triangle(opt))
				return 1;
			comb_sortedsketch64_self_matrix(opt);
			return 1;
		}
		free_read_from_file(qry_dim_sketch_stat, file_size);
		comb_sortedsketch64_self_matrix(opt);
		return 1;
	}
	if (opt->fmt && paths_refer_to_same_file(opt->refdir, opt->qrydir))
	{
		if (opt->fmt == 1 &&
			(force_ref_index || ani_matrix_should_use_index(qry_infile_num, qry_infile_num, matrix_direct_n)))
		{
			if (force_ref_index && !file_exists_in_folder(opt->qrydir, sorted_comb_ctxgid64obj32))
				errx(EXIT_FAILURE, "KSSD3A_FORCE_REF_INDEX requires reference index '%s/%s'",
					 opt->qrydir, sorted_comb_ctxgid64obj32);
			else
				ensure_sorted_ref_index_for_ani_m1(opt->qrydir);
			free_read_from_file(qry_dim_sketch_stat, file_size);
			if (comb_sortedsketch64_indexed_self_full(opt))
				return 1;
			comb_sortedsketch64_self_matrix(opt);
			return 1;
		}
		if (opt->fmt == 2 &&
			(force_ref_index || ani_matrix_should_use_index(qry_infile_num, qry_infile_num, matrix_direct_n)))
		{
			free_read_from_file(qry_dim_sketch_stat, file_size);
			if (comb_sortedsketch64_indexed_self_triangle(opt))
				return 1;
			comb_sortedsketch64_self_matrix(opt);
			return 1;
		}
		free_read_from_file(qry_dim_sketch_stat, file_size);
		comb_sortedsketch64_self_matrix(opt);
		return 1;
	}
	if (opt->fmt == 2)
		errx(EXIT_FAILURE, "ANI triangle output requires one sketch or identical -r/-q sketches");
	const int ref_infile_num = sketch_infile_num(opt->refdir);
	const bool use_small_query_stream =
		!force_ref_index && opt->fmt == 0 && ref_infile_num > 0 &&
		qry_infile_num < ref_infile_num;
	if (use_small_query_stream && qry_infile_num == 1)
	{
		free_read_from_file(qry_dim_sketch_stat, file_size);
		return stream_ref_sketches_one_qraw_lookup(opt);
	}
	const char *multi_qraw_mode = getenv("KSSD3A_QRAW_MULTI");
	if (opt->fmt == 0 && qry_infile_num > 1)
	{
		bool use_query_index = use_small_query_stream &&
							   qry_infile_num <= qraw_multi_query_index_threshold();
		if (opt->unassembled && !force_ref_index && multi_qraw_mode && multi_qraw_mode[0] != '\0')
		{
			if (strcmp(multi_qraw_mode, "sortedindex") == 0)
				use_query_index = true;
			else if (strcmp(multi_qraw_mode, "legacy") == 0 ||
					 strcmp(multi_qraw_mode, "refindex") == 0)
				use_query_index = false;
			else
				errx(EXIT_FAILURE, "invalid KSSD3A_QRAW_MULTI='%s' (use sortedindex, refindex, or legacy)",
					 multi_qraw_mode);
		}
		if (use_query_index)
		{
			free_read_from_file(qry_dim_sketch_stat, file_size);
			return stream_ref_sketches_multi_qraw_sortedindex(opt);
		}
	}
	free_read_from_file(qry_dim_sketch_stat, file_size);
	const bool ref_index_exists = file_exists_in_folder(opt->refdir, sorted_comb_ctxgid64obj32);
	if (opt->fmt == 1)
	{
		if (force_ref_index)
		{
			if (!ref_index_exists)
				errx(EXIT_FAILURE, "KSSD3A_FORCE_REF_INDEX requires reference index '%s/%s'",
					 opt->refdir, sorted_comb_ctxgid64obj32);
			return mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(opt);
		}
		if (ani_matrix_should_use_index(ref_infile_num, qry_infile_num, matrix_direct_n))
		{
			ensure_sorted_ref_index_for_ani_m1(opt->refdir);
			return mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(opt);
		}
	}
	else if (ref_index_exists)
		return mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(opt);
	const int auto_index_n = auto_ref_index_threshold();
	if (!opt->v && !opt->unassembled && opt->fmt == 0 &&
		auto_index_n > 0 && ref_infile_num > auto_index_n)
	{
		gen_inverted_index4comblco(opt->refdir);
		return mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(opt);
	}
	if (force_ref_index)
		errx(EXIT_FAILURE, "KSSD3A_FORCE_REF_INDEX requires reference index '%s/%s'",
			 opt->refdir, sorted_comb_ctxgid64obj32);

	if (opt->fmt)
	{
		comb_sortedsketch64Xcomb_sortedsketch64(opt);
		return 1;
	}
	comb_manysmall_sortedsketch64Xcomb_fewlarge_sortedsketch64_filter_and_sort_survivors(opt);
	return 1;
}

static bool select_auto_template_from_inputs(const path_vec_t *inputs, dim_sketch_stat_t *template_stat)
{
	for (int i = 0; i < inputs->n; ++i)
	{
		const char *path = inputs->items[i];
		if (!is_lco_sketch_dir(path))
			continue;
		*template_stat = read_lco_sketch_stat(path);
		validate_supported_auto_sketch_stat(template_stat);
		return true;
	}
	return false;
}

static dim_sketch_stat_t select_auto_template_stat(const ani_opt_t *opt, const path_vec_t *ref_inputs, const path_vec_t *qry_inputs)
{
	dim_sketch_stat_t template_stat = default_auto_ani_sketch_stat(opt);
	bool have_template = select_auto_template_from_inputs(ref_inputs, &template_stat);
	if (!have_template)
		have_template = select_auto_template_from_inputs(qry_inputs, &template_stat);
	if (!have_template)
		validate_supported_auto_sketch_stat(&template_stat);
	return template_stat;
}

static void validate_auto_sketch_inputs(const path_vec_t *inputs, const dim_sketch_stat_t *template_stat)
{
	for (int i = 0; i < inputs->n; ++i)
	{
		const char *path = inputs->items[i];
		if (!is_lco_sketch_dir(path))
			continue;
		dim_sketch_stat_t stat = read_lco_sketch_stat(path);
		if (!same_ani_auto_sketch_stat(template_stat, &stat))
			describe_sketch_stat_mismatch(path, template_stat, &stat);
	}
}

static int path_vec_stdin_count(const path_vec_t *inputs)
{
	int count = 0;
	for (int i = 0; i < inputs->n; ++i)
		if (strcmp(inputs->items[i], "-") == 0)
			count++;
	return count;
}

static bool auto_input_is_streamed_raw(const ani_opt_t *opt, const char *path)
{
	return strcmp(path, "-") == 0 || opt->sketch_pipecmd[0] != '\0';
}

static bool auto_input_is_fastq_path(const char *path)
{
	if (strcmp(path, "-") == 0 || is_lco_sketch_dir(path))
		return false;
	return isOK_fmt_infile((char *)path, fastq_fmt, FQ_FMT_SZ);
}

static bool auto_query_side_is_unassembled(const ani_opt_t *opt, const path_vec_t *qry_inputs)
{
	for (int i = 0; i < qry_inputs->n; ++i)
	{
		const char *path = qry_inputs->items[i];
		if (is_lco_sketch_dir(path))
			continue;
		if (auto_input_is_fastq_path(path))
			return true;
		if (opt->sketch_conflict && auto_input_is_streamed_raw(opt, path))
			return true;
	}
	return false;
}

static void prepare_auto_input_side(
	const path_vec_t *inputs,
	path_vec_t *sketches,
	const char *side_prefix,
	const dim_sketch_stat_t *template_stat,
	const char *tmp_root,
	const ani_opt_t *ani_opt_val)
{
	for (int i = 0; i < inputs->n; ++i)
	{
		const char *input = inputs->items[i];
		if (is_lco_sketch_dir(input))
		{
			path_vec_push(sketches, input);
			continue;
		}

		char outdir[PATHLEN * 2];
		int written = snprintf(outdir, sizeof(outdir), "%s/%sseq_%d", tmp_root, side_prefix, i);
		if (written < 0 || (size_t)written >= sizeof(outdir))
			errx(EXIT_FAILURE, "auto ANI temporary output path is too long for %s", input);
		sketch_sequence_arg_to_dir(input, template_stat, outdir, ani_opt_val);
		path_vec_push(sketches, outdir);
	}
}

static void copy_runtime_path(char *dest, size_t dest_size, const char *path, const char *label)
{
	if (strlen(path) >= dest_size)
		errx(EXIT_FAILURE, "%s path is too long; maximum supported length is %zu bytes", label, dest_size - 1);
	snprintf(dest, dest_size, "%s", path);
}

static int run_auto_ani_from_sides(ani_opt_t *opt, const path_vec_t *ref_inputs, const path_vec_t *qry_inputs)
{
	if (ref_inputs->n < 1)
		errx(EXIT_FAILURE, "auto ANI requires at least one reference input");
	if (qry_inputs->n < 1)
		errx(EXIT_FAILURE, "auto ANI requires at least one query input");
	if (path_vec_stdin_count(ref_inputs) + path_vec_stdin_count(qry_inputs) > 1)
		errx(EXIT_FAILURE, "stdin input '-' can be used only once in auto ANI");

	dim_sketch_stat_t template_stat = select_auto_template_stat(opt, ref_inputs, qry_inputs);
	validate_auto_sketch_inputs(ref_inputs, &template_stat);
	validate_auto_sketch_inputs(qry_inputs, &template_stat);

	char tmp_template[] = "/tmp/kssd3a_ani_auto_XXXXXX";
	char *tmp_root = mkdtemp(tmp_template);
	if (!tmp_root)
		err(errno, "%s(): mkdtemp", __func__);

	path_vec_t ref_sketches = {0};
	path_vec_t qry_sketches = {0};
	prepare_auto_input_side(ref_inputs, &ref_sketches, "ref", &template_stat, tmp_root, opt);
	prepare_auto_input_side(qry_inputs, &qry_sketches, "qry", &template_stat, tmp_root, opt);
	if (ref_sketches.n < 1 || qry_sketches.n < 1)
		errx(EXIT_FAILURE, "auto ANI failed to prepare reference/query sketches");

	char ref_merge[PATHLEN * 2];
	char qry_merge[PATHLEN * 2];
	int written = snprintf(ref_merge, sizeof(ref_merge), "%s/ref_merged", tmp_root);
	if (written < 0 || (size_t)written >= sizeof(ref_merge))
		errx(EXIT_FAILURE, "auto ANI reference merge path is too long");
	written = snprintf(qry_merge, sizeof(qry_merge), "%s/qry_merged", tmp_root);
	if (written < 0 || (size_t)written >= sizeof(qry_merge))
		errx(EXIT_FAILURE, "auto ANI query merge path is too long");
	char *refdir = merge_auto_sketches(&ref_sketches, ref_merge);
	char *qrydir = merge_auto_sketches(&qry_sketches, qry_merge);

	char old_refdir[PATHLEN];
	char old_qrydir[PATHLEN];
	const bool old_unassembled = opt->unassembled;
	const bool old_v = opt->v;
	get_generic_dist_from_features_fn old_dist_fn = get_generic_dist_from_features;
	bool auto_unassembled = auto_query_side_is_unassembled(opt, qry_inputs);
	snprintf(old_refdir, sizeof(old_refdir), "%s", opt->refdir);
	snprintf(old_qrydir, sizeof(old_qrydir), "%s", opt->qrydir);
	copy_runtime_path(opt->refdir, sizeof(opt->refdir), refdir, "auto ANI reference sketch");
	copy_runtime_path(opt->qrydir, sizeof(opt->qrydir), qrydir, "auto ANI query sketch");
	opt->unassembled = auto_unassembled;
	if (auto_unassembled)
	{
		opt->v = true;
		get_generic_dist_from_features = get_naive_dist;
	}

	int rc = run_sketch_ani(opt);

	snprintf(opt->refdir, sizeof(opt->refdir), "%s", old_refdir);
	snprintf(opt->qrydir, sizeof(opt->qrydir), "%s", old_qrydir);
	opt->unassembled = old_unassembled;
	opt->v = old_v;
	get_generic_dist_from_features = old_dist_fn;
	free(refdir);
	free(qrydir);
	path_vec_free(&ref_sketches);
	path_vec_free(&qry_sketches);
	remove_tree(tmp_root);
	return rc;
}

static int run_auto_positional_ani(ani_opt_t *opt)
{
	if (opt->num_remaining_args < 2)
		errx(EXIT_FAILURE, "auto ANI requires one reference input followed by at least one query input");

	path_vec_t ref_inputs = {0};
	path_vec_t qry_inputs = {0};
	path_vec_push(&ref_inputs, opt->remaining_args[0]);
	for (int i = 1; i < opt->num_remaining_args; ++i)
		path_vec_push(&qry_inputs, opt->remaining_args[i]);

	int rc = run_auto_ani_from_sides(opt, &ref_inputs, &qry_inputs);
	path_vec_free(&ref_inputs);
	path_vec_free(&qry_inputs);
	return rc;
}

static int run_auto_list_ani(ani_opt_t *opt)
{
	path_vec_t ref_inputs = {0};
	path_vec_t qry_inputs = {0};
	read_path_list_to_vec(opt->reflist, &ref_inputs, "--reflist");
	read_path_list_to_vec(opt->qrylist, &qry_inputs, "--qrylist");

	int rc = run_auto_ani_from_sides(opt, &ref_inputs, &qry_inputs);
	path_vec_free(&ref_inputs);
	path_vec_free(&qry_inputs);
	return rc;
}

static bool path_is_directory(const char *path)
{
	struct stat st;
	return path && stat(path, &st) == 0 && S_ISDIR(st.st_mode);
}

static int run_direct_query_sequence_ani(ani_opt_t *opt, bool force_raw_query, const char *option_name)
{
	if (!is_lco_sketch_dir(opt->refdir))
		errx(EXIT_FAILURE,
			 "%s sequence input requires -r/--ref to be a sketch directory; "
			 "for direct reference sequence input use positional auto ANI",
			 option_name);
	if (path_is_directory(opt->qrydir))
		errx(EXIT_FAILURE, "%s input '%s' is a directory but not a valid sketch", option_name, opt->qrydir);

	dim_sketch_stat_t qry_template = read_lco_sketch_stat(opt->refdir);
	validate_supported_auto_sketch_stat(&qry_template);

	char input[PATHLEN];
	char old_qrydir[PATHLEN];
	snprintf(input, sizeof(input), "%s", opt->qrydir);
	snprintf(old_qrydir, sizeof(old_qrydir), "%s", opt->qrydir);
	const bool raw_query = force_raw_query;
	const bool keep_query_conflicts = force_raw_query || opt->sketch_conflict;
	if (keep_query_conflicts)
		qry_template.conflict = true;

	char tmp_template[] = "/tmp/kssd3a_ani_query_XXXXXX";
	char *tmp_root = mkdtemp(tmp_template);
	if (!tmp_root)
		err(errno, "%s(): mkdtemp", __func__);

	char qry_out[PATHLEN * 2];
	int written = snprintf(qry_out, sizeof(qry_out), "%s/qry_seq", tmp_root);
	if (written < 0 || (size_t)written >= sizeof(qry_out))
		errx(EXIT_FAILURE, "temporary %s sketch output path is too long", option_name);

	sketch_sequence_arg_to_dir(input, &qry_template, qry_out, opt);

	const bool old_unassembled = opt->unassembled;
	const bool old_v = opt->v;
	char old_refdir[PATHLEN];
	char old_outf[PATHLEN];
	snprintf(old_refdir, sizeof(old_refdir), "%s", opt->refdir);
	snprintf(old_outf, sizeof(old_outf), "%s", opt->outf);
	get_generic_dist_from_features_fn old_dist_fn = get_generic_dist_from_features;
	copy_runtime_path(opt->qrydir, sizeof(opt->qrydir), qry_out, "temporary query sketch");
	opt->unassembled = raw_query;
	if (raw_query)
	{
		opt->v = true;
		get_generic_dist_from_features = get_naive_dist;
	}

	int rc = run_sketch_ani(opt);

	snprintf(opt->qrydir, sizeof(opt->qrydir), "%s", old_qrydir);
	snprintf(opt->refdir, sizeof(opt->refdir), "%s", old_refdir);
	snprintf(opt->outf, sizeof(opt->outf), "%s", old_outf);
	opt->unassembled = old_unassembled;
	opt->v = old_v;
	get_generic_dist_from_features = old_dist_fn;
	remove_tree(tmp_root);
	return rc;
}

int cmd_ani(struct argp_state *state)
{

	struct arg_ani ani = {
		0,
	};
	int argc = state->argc - state->next + 1;
	char **argv = &state->argv[state->next - 1];
	ani.global = state->input;
	argp_parse(&argp_ani, argc, argv, ARGP_IN_ORDER, &argc, &ani);
	state->next += argc - 1;

	// instantilize distance fn according selection of naive model or not
	ani_model_drfold = ani_opt.sketch_drfold;
	get_generic_dist_from_features = ani_opt.v ? get_naive_dist : lm3ways_dist_from_features;

	if (ani_opt.reflist[0] != '\0' || ani_opt.qrylist[0] != '\0')
		return run_auto_list_ani(&ani_opt);
	else if (ani_opt.qrydir[0] != '\0')
	{
		if (!is_lco_sketch_dir(ani_opt.qrydir))
			return run_direct_query_sequence_ani(&ani_opt,
												 ani_opt.unassembled,
												 ani_opt.unassembled ? "--qraw" : "-q/--query");
		return run_sketch_ani(&ani_opt);
	}
	else if (ani_opt.num_remaining_args >= 2)
		return run_auto_positional_ani(&ani_opt);
	return 1;
}
