#include "command_composite.h"
#include "global_basic.h"
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
#include <limits.h>

/*** argp wrapper ***/
struct arg_composite
{
  struct arg_global* global;

  char* name;
};

static struct argp_option opt_composite[] =
{
	{"ref",'r',"<DIR>", 0, "Species-specific pan unique k-mer database directory.",1},
	{"query",'q',"<DIR>", 0, "Query sketches with abundances.",2},
	{"outfile",'o',"<PATH>",0,"Profile output file; with --binVec, output directory.",3},
	{"threads",'p',"<INT>", 0, "Number of threads to use.",4},
	{"binVec",'b',0,0,"Output species abundances in binary vector format (.abv).",5},
	{"idxbv",'i',0,0,"Build index of abundance binary vectors.",6},
	{"search",'s',"<0-2>",0,"Search similar abundance binary vectors using cosine(0), L1 norm(1), or L2 norm(2).",7},
	{"readabv",'d',0,0,"Read an .abv file.",8},
  { 0 }
};

static char doc_composite[] =
  "\n"
  "Estimate metagenomic composition and work with abundance binary vectors."
  "\v"
  "Use -q for abundance estimation, -i for indexing, -s for search, or -d to read .abv files.\n"
  "\n"
  "Examples:\n"
  "  kssd3a composite -r markerdb -q query_abundance_sketches -o profile_out\n"
  "  kssd3a composite -r profile_db -i\n"
  "  kssd3a composite -r profile_db -s 0 query.abv"
  ;

composite_opt_t composite_opt ={
	.b = 0, //write out abundance binary vector (1) or not (0)
	.i = 0, // index abundance binary vectors (1) or not (0)
	.s = -1, 	
	.d = 0, // read .abv file
	.p = 1,
	.refdir[0] = '\0',
	.qrydir[0] = '\0',
	.outdir = "./",
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

static void copy_path_arg(struct argp_state *state, const char *option_name,
                          char *dest, size_t dest_size, const char *arg)
{
  if (strlen(arg) >= dest_size)
    argp_error(state, "%s path is too long; maximum supported length is %zu bytes", option_name, dest_size - 1);
  snprintf(dest, dest_size, "%s", arg);
}

static error_t parse_composite(int key, char* arg, struct argp_state* state) {
  struct arg_composite* composite = state->input;
  assert( composite );
  assert( composite->global );
	
  switch(key)
  {
		case 'b':
		{
			composite_opt.b = 1;
			break;
		}
		case 'i':
		{
			 composite_opt.i = 1;
				break;
		}
		case 's':
		{
				composite_opt.s = parse_int_range(state, "-s/--search", arg, 0, 2);
				break;
		}
		case 'd':
		{
			composite_opt.d = 1;
			break;
		}
		case 'p':
		{
			composite_opt.p = parse_int_range(state, "-p/--threads", arg, 1, 65536);
			break;
		}
		case 'r':
		{
			copy_path_arg(state, "-r/--ref", composite_opt.refdir, sizeof(composite_opt.refdir), arg);
			break;
		}
		case 'q':
		{
			copy_path_arg(state, "-q/--query", composite_opt.qrydir, sizeof(composite_opt.qrydir), arg);
			break;
		}
		case 'o':
		{
			copy_path_arg(state, "-o/--outfile", composite_opt.outdir, sizeof(composite_opt.outdir), arg);
			break;
		}
		case ARGP_KEY_ARGS:
		{
			composite_opt.num_remaining_args = state->argc - state->next;
			composite_opt.remaining_args  = state->argv + state->next;
			break;
		}
    case ARGP_KEY_NO_ARGS:
    {	
			if(state->argc<2)
			{
      	printf("\v");
				argp_state_help(state,stdout,ARGP_HELP_SHORT_USAGE);
				printf("\v");
      	argp_state_help(state,stdout,ARGP_HELP_LONG);
      	printf("\v");
      	return EINVAL;
			}
    }
		break;
    default:
      return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

static struct argp argp_composite =
{
  opt_composite,
  parse_composite,
	0,//  "[arguments ...]",
  doc_composite
};

int cmd_composite(struct argp_state* state)
{
  struct arg_composite composite = { 0, };
  int    argc = state->argc - state->next + 1;
  char** argv = &state->argv[state->next - 1];
  composite.global = state->input;
	
  argp_parse(&argp_composite, argc, argv, ARGP_IN_ORDER, &argc, &composite);
	
  state->next += argc - 1;
	if( composite_opt.refdir[0] != '\0' ){
		if (composite_opt.qrydir[0] != '\0') //if(argc >1)
			return get_species_abundance (&composite_opt);
		else if (composite_opt.i) 
			return index_abv (&composite_opt);
		else if (composite_opt.s != -1){
			if( composite_opt.s >=0 && composite_opt.s <3 && composite_opt.num_remaining_args >0 )	return abv_search(&composite_opt);
			else printf("\vUsage: %s composite -r <ref> -s <0|1|2> <query.abv>\n", state->name);
		}
		else printf("\vUsage: %s composite -r <ref> < mode: -q | -i | -s >\n",state->name);
	}
	else if(composite_opt.d){
		if(composite_opt.num_remaining_args < 1)	printf("\vUsage: %s composite -d <query.abv>\n", state->name);
		else return read_abv (&composite_opt);	
	}
	else
		printf("\vUsage: %s composite -r <ref> < mode: -q | -i | -s >\n", state->name);	

	return 1;
}//cmd_composite();
