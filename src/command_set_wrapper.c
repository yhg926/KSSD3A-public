#include "command_set_wrapper.h"
#include "command_set.h"
#include "command_operate.h"
/*** argp wraper ***/
#define SET_ARG_KEY_MODE 890

struct arg_set
{
  struct arg_global* global;

  char* name;
};

static struct argp_option opt_set[] =
{
  {"union",'u', 0,  0, "Compute the union set of sketches.",1 },
	{"subtract",'s',"<pan>", 0,"Subtract the pan-sketch from each input sketch.",2 },
	{"intsect",'i',"<pan>", 0, "Intersect each input sketch with the pan-sketch.",2},
	{"intersect",334,"<pan>", 0, "Alias for --intsect.",2},
	{"key",SET_ARG_KEY_MODE,"<full|ctx>", 0, "Set operation key. full compares context+object; ctx compares context only for --intersect/--subtract. [full]",2},
	{"uniq_union",'q',0,  0, "Compute the unique union set of sketches.",3 },
    {"markerdb",333,0,  0, "Generate marker database instead of unique union set. Requires -q.",4 },
//	{"combin_pan",'c',0,  0, "combine pan files to combco file.",4 },
	{"threads",'p',"<INT>",  0, "Number of threads.",5 },
	{"print",'P',0,  0, "Print genome names.",5 },
	{"psketch",777,0,  0, "Print sketch content.",5 },
	{"pindex",888,0,  0, "Print context/genome/object index content.",5 },
	{"ppos",889,0,  0, "Print sketch positions as sample, sketch entry, and zero-based position.",5 },
	{"grouping",'g',"<file.tsv>",0,"Group genomes by an input category file.",5},
	{"outdir",'o',"<path>",0,"Output directory.",6},
  { 0 }
};

static char doc_set[] =
  "\n"
  "Run set operations on combined sketches."
  "\v"
  "Choose one operation such as --union, --uniq_union, --intersect, --intsect, or --subtract.\n"
  "Use --key ctx with --intersect/--subtract on -T long sketches to match by context only while preserving original context-object records in output.\n"
  "\n"
  "Examples:\n"
  "  kssd3a set --union -o union_sketch input_sketches\n"
  "  kssd3a set --uniq_union --markerdb -o markerdb input_sketches\n"
  "  kssd3a set --intersect pan_sketch --key ctx -o intersected input_sketches\n"
  "  kssd3a set --subtract pan_sketch --key ctx -o subtracted input_sketches"
  ;


set_opt_t set_opt = {
.operation = -1,//0:subtract,1:intersect,2 union, 3 uniq_union, 4 combin_pan
.key_mode = SET_KEY_FULL,
.q2markerdb = 0, // when -q set, generate markerdb instead of uniq union set, only for lco sketch 
.p = 1,
.P = 0,
.show = 0,
.num_remaining_args = 0,
.remaining_args = NULL,			
.insketchpath[0] = '\0',
.pansketchpath[0]='\0',
.subsetf[0] = '\0',
.outdir = "./"
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

static error_t parse_set(int key, char* arg, struct argp_state* state) {

  struct arg_set* set = state->input;
  assert( set );
  assert( set->global );
	
  switch(key)
  {

    case 'u':
		{ 
			if (set_opt.operation != -1 ) printf("set operation is already set, -u is ignored.\n");
			else set_opt.operation = 2 ;
			break;
		}
		case 'q':
    {
			
      if (set_opt.operation != -1 ) printf("set operation is already set, -q is ignored.\n");
      else set_opt.operation = 3 ;
      break;
    }
		case 's':
		{
			if (set_opt.operation != -1) printf("set operation is already set, -s is ignored.\n");
			else {
				set_opt.operation = 0;
				copy_path_arg(state, "-s/--subtract", set_opt.pansketchpath, sizeof(set_opt.pansketchpath), arg);
			}
			break;
		}
		case 'i':
		case 334:
		{
		
			if (set_opt.operation != -1 ) printf("set operation is already set, --intersect is ignored.\n");
			else {
				set_opt.operation = 1 ;            
				copy_path_arg(state, "-i/--intersect", set_opt.pansketchpath, sizeof(set_opt.pansketchpath), arg);
			}
			break;
		}
		case SET_ARG_KEY_MODE:
		{
			if (strcmp(arg, "full") == 0)
				set_opt.key_mode = SET_KEY_FULL;
			else if (strcmp(arg, "ctx") == 0)
				set_opt.key_mode = SET_KEY_CTX;
			else
				argp_error(state, "--key must be one of: full, ctx");
			break;
		}
		case 'c':
		{
			if (set_opt.operation != -1 ) printf("set operation is already set, -c is ignored.\n");
			else set_opt.operation = 4 ;
			break;
		}
		case 'o':
		{
		
			copy_path_arg(state, "-o/--outdir", set_opt.outdir, sizeof(set_opt.outdir), arg);
			 
			break;
		}
		case 'p':
		{
			set_opt.p = parse_int_range(state, "-p/--threads", arg, 1, 65536);
			break;
		}
		case 'P':
		{
			set_opt.P = 1;	
			break;
		}
		case 'g':
		{
			copy_path_arg(state, "-g/--grouping", set_opt.subsetf, sizeof(set_opt.subsetf), arg);
			break;
		}
		case 333:
		{
			set_opt.q2markerdb = 1; 
			break;
		}
		case 777:
		{
			set_opt.show = 1;
			break;	
		}
		case 888:
        {
            set_opt.show = 2;
            break;
        }
		case 889:
        {
            set_opt.show = 3;
            break;
        }
		case ARGP_KEY_ARGS:
			copy_path_arg(state, "<combined sketch>", set_opt.insketchpath, sizeof(set_opt.insketchpath), state->argv[state->next]);
			set_opt.num_remaining_args = state->argc - state->next;
			set_opt.remaining_args  = state->argv + state->next;
			break;
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

static struct argp argp_set =
{
  opt_set,
  parse_set,
	"<combined sketch>", //0  "[arguments ...]",
  doc_set
};

int cmd_set(struct argp_state* state)
{
  struct arg_set set = { 0, };
  int    argc = state->argc - state->next + 1;
  char** argv = &state->argv[state->next - 1];
  char*  argv0 =  argv[0];

  set.global = state->input;
	argv[0] = malloc(strlen(state->name) + strlen(" set") + 1);

  if(!argv[0])
    argp_failure(state, 1, ENOMEM, 0);
	sprintf(argv[0], "%s set", state->name);
  argp_parse(&argp_set, argc, argv, ARGP_IN_ORDER, &argc, &set);

	free(argv[0]);
  argv[0] = argv0;
  state->next += argc - 1;
	// operation and arg control
	if(argc >1){	
		if(set_opt.key_mode == SET_KEY_CTX && set_opt.operation != 0 && set_opt.operation != 1)
			errx(EXIT_FAILURE, "--key ctx currently supports only --intersect/--intsect and --subtract");
		if(set_opt.operation == 2){
			if(file_exists_in_folder(set_opt.insketchpath,co_dstat) )
				return sketch_union(&set_opt); 
			else if(file_exists_in_folder(set_opt.insketchpath,sketch_stat))
				return lsketch_union(&set_opt);			
		}
		else if(set_opt.operation == 3){
	  	if(file_exists_in_folder(set_opt.insketchpath,co_dstat))
				return uniq_sketch_union(&set_opt) ;
      else if(file_exists_in_folder(set_opt.insketchpath,sketch_stat))
        return lsketch_union(&set_opt);
		}
		else if(set_opt.operation == 4){
			return combin_pans(&set_opt);
		}
		else if(set_opt.operation == 0 || set_opt.operation == 1 ){
			if(file_exists_in_folder(set_opt.pansketchpath,co_dstat)) {
				if(set_opt.key_mode == SET_KEY_CTX)
					errx(EXIT_FAILURE, "--key ctx requires long -T sketches and does not support legacy cofiles sketches");
				return sketch_operate(&set_opt) ;
			}
			else if (file_exists_in_folder(set_opt.pansketchpath,sketch_stat))
				 return lsketch_operate(&set_opt) ;
		}
		else {
			if(set_opt.P) {
				if(file_exists_in_folder(set_opt.insketchpath,co_dstat)) print_gnames(&set_opt);
				else if(file_exists_in_folder(set_opt.insketchpath,sketch_stat)) print_lco_gnames(&set_opt);
				else printf("%s is not a valid sketch\n",set_opt.insketchpath );

			}
			else if(set_opt.show > 0){
					if(file_exists_in_folder(set_opt.insketchpath,sketch_stat)) show_content(&set_opt);
			}
			else if (set_opt.subsetf[0]!='\0') {
				if(file_exists_in_folder(set_opt.insketchpath,co_dstat))
					return grouping_genomes(&set_opt); // combin_subset_pans(set_opt.subsetf);
				else if(file_exists_in_folder(set_opt.insketchpath,sketch_stat))
					return lgrouping_genomes(&set_opt);

			}
			else printf("set operation use : -u, -q, -i or -s\n");
			return -1 ;
		}
	}
	else
		return -1;
}
