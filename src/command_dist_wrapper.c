#include "command_dist_wrapper.h"
#include "command_dist.h"
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

#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
#endif

/*** argp wraper ***/
struct arg_dist
{
  struct arg_global* global;

  char* name;
};

static struct argp_option opt_dist[] =
{
	{"halfKmerlength",'k',"INT",0, "Half k-mer length, from 2 to 15. [8]" },
	{"threadN",'p',"INT",0,"Number of threads to use. [all threads]"},
	{"list",'l',"FILE",0,"File containing input sequence paths."},
	{"DimRdcLevel",'L',"INT|FILE",0,"Dimensionality-reduction level, or path to a .shuf file. [2]"},
	{"maxMemory",'m',"NUM",0,"Maximum memory to use, in GB."},
	{"LstKmerOcrs",'n',"INT",0,"Least k-mer occurrence required in FASTQ input."},
	{"quality",'Q',"INT",0,"Filter k-mers with lowest base quality below this Phred score. [0]"},
	{"reference_dir",'r',"<path>",0,"Reference genome or database directory to search against."},
//{"distance",'d',0,0,"caculate pairwise distance."},
	{"outdir",'o',"<path>",0,"Output directory for result files." },
	{"neighborN_max",'N',"INT",0,"Maximum number of nearest reference genomes. [1]"},
	{"mutDist_max",'D',"FLT",0,"Maximum mutation distance to output. [1]"},
	{"metric",'M',"0/1",0,"Output metric: Jaccard(0) or containment(1). [0]"},
	{"containment",'C',"0/1/2",0,"Containment denominator: Min(0), Ref(1), or Qry(2). [0]"},
	{"outfields",'O',"0/1/2",0,"Output fields: Distance(0), Q-values(1), confidence intervals(2). Higher values include earlier fields. [2]"},
	{"correction",333,"0/1",0,"Correct shared k-mer counts. [0]" },
  {"abundance",'A',0,0,"Enable abundance-estimation mode."},
	{"dedup",'u',0,0,"Ignore repeated k-mers in the reference."},
	{"keepcofile",888,0,0,"Keep intermediate .co files."},
	{"pipecmd",'P',"<cmd>",0,"Pipe command for reading input."},
	{"keepskf",777,0,0,"Keep shared k-mer count file. [false]"},
	{"skf",'f',"<skfpath>",0,"Shared k-mer count file path."},
	{"byread",555,0,0,"Sketch the input by read. [false]"},
//{"onlyMashD",222,0,0,"only print mash distance."},
//{"stage2",999,0,0,"input is intermedia .co files."},
  { 0 }
};

static char doc_dist[] =
  "\n"
  "Legacy cofile-style distance command. For current KSSD3A sketch directories, use `matrix` for distances or `ani` for ANI estimates."
  "\v"
  "This compatibility command expects older cofiles.stat/mcofiles.stat-style inputs and may not accept sketch directories produced by the current `sketch` command.\n"
  "If -r is omitted, KSSD3 computes pairwise distances among query inputs.\n"
  "\n"
  "Examples:\n"
  "  kssd3a dist -r ref_sketches qry_sketches\n"
  "  kssd3a dist -M 1 -C 0 -r ref_sketches qry_sketches\n"
  "  kssd3a dist -k 8 -L 2 -p 8 genomes/*.fasta"
  ;
//options collector for command dist
dist_opt_val_t dist_opt_val = 
{ 
.k = 8,	//half kmer len: k
.p = 0,  // threads num: p
.dr_level = 2,  // dr_level;//dimension reduction level
.dr_file = "", //char dr_file[PATHLEN];//dimension reduction file
.mmry = 0, // double mmry; //maxMemory;
.fmt = "mfa", 
.refpath = "", //char refpath[PATHLEN]; //reference sequences path
.fpath = "",	//char fpath[PATHLEN]; //query files path
.outdir = ".", //char outdir[PATHLEN]; // results dir
.kmerocrs = 1, //int kmerocrs;fastq file Kmer least occurence  
.kmerqlty = 0, //int kmerqlty; fastq file Kmer quality threshold
.keepco = false, // bool keepco; wether keep intermidia .co file or not 
.stage2 = false, // input is intermeida .co folder
.num_neigb = 0, //neighborN_max, 0 means all
.mut_dist_max = 1, //mutDist_max
.metric = Jcd,
.ctm = 0, // containment option: min / query / reference
.outfields = CI,
.correction = false,
.abundance = false, // no abundance
.u = false, // k-mer dedup off
.pipecmd = "", // no pipe command 
.shared_kmerpath="", // share_kmer_ct file path
.keep_shared_kmer=false, //not keep share_kmer_ct file 
.byread=false,	
.num_remaining_args = 0, //int num_remaining_args; no option arguments num. 
.remaining_args = NULL //char **remaining_args; no option arguments array.
} ;

const char outdir_name[] = "kssd_rslt" ;

static error_t parse_dist(int key, char* arg, struct argp_state* state) {
  struct arg_dist* dist = state->input;
  assert( dist );
  assert( dist->global );
  switch(key)
  {
		case 'k':
			dist_opt_val.k = atoi(arg);	
			break;	
		case 'p':
		{	
#ifdef _OPENMP
				dist_opt_val.p = atoi(arg) ;
#else
   		warnx("This version of kssd was built without OpenMP and "
          "thus does not support multi threading. Ignoring -p %d",atoi(arg));
      break;
#endif
		}
			break;

		case 'm':
		{	
			double sys_mm = get_sys_mmry();
			double rqst_mm = atoi(arg);
			if ( rqst_mm > sys_mm ){
				warnx("Memory request is larger than system available %f. Ignoring -m %f",sys_mm,rqst_mm);
				dist_opt_val.mmry = sys_mm;
			}
			else 
				dist_opt_val.mmry = rqst_mm ;
		}
		break;
		
		case 'r':
		{
			if(strlen(arg) > PATHLEN) {
        err(errno,"the list path should not longer than %d", PATHLEN);
        exit(EXIT_FAILURE);
      };
			strcpy(dist_opt_val.refpath,arg);			
			break;
		}
		case 'l':
		{
			if(strlen(arg) > PATHLEN) {
				err(errno,"the list path should not longer than %d", PATHLEN); 
				exit(EXIT_FAILURE);
			};		
			strcpy(dist_opt_val.fpath,arg);
			break;
		}
		case 'L':
		{
			struct stat path_stat;
			if( stat(arg,&path_stat) >=0 && S_ISREG(path_stat.st_mode)){
				if(strlen(arg) < PATHLEN )
					strcpy(dist_opt_val.dr_file,arg);
				else
					 err(errno,"-L argument path should not longer than %d",PATHLEN);			
			}
			else{
				if ( atoi(arg) >= dist_opt_val.k - 2 || atoi(arg) < 0 ) 
					err(errno,"-L: dimension reduction level should never larger than Kmer length - 2,"
								" which is %d here",dist_opt_val.k - 2 );  			
				dist_opt_val.dr_level = atoi(arg);
			} 
			break;
		}
		case 'n':
		{
			if( atoi(arg) > 65536 ) {
				dist_opt_val.kmerocrs = 65536;
				warnx("-n argument is larger than Max, it has been set to 7, ignorned -n %d ",atoi(arg));
			}else if( atoi(arg) < 1 ){
				dist_opt_val.kmerocrs = 1;
				 warnx("-n argument is smaller than Min, it has been set to 1, ignorned -n %d ",atoi(arg));
			} 

			else dist_opt_val.kmerocrs = atoi(arg);
		}
			break;
		case 'Q':
		{
			dist_opt_val.kmerqlty = atoi(arg);
		}	
			break;
    case 'o':
		{
      if(strlen(arg) > PATHLEN) {
        err(errno,"the outdir path should not longer than %d", PATHLEN);
        exit(EXIT_FAILURE);
      };
      strcpy(dist_opt_val.outdir,arg);
			break;
    }
		case 'N':
		{
			dist_opt_val.num_neigb = atoi(arg);
			break;
		}
		case 'D':
		{
			dist_opt_val.mut_dist_max = atof(arg);
			break;
		}
		case 'A':
		{
			dist_opt_val.abundance = true;
			break;
		}
		case 'u':
    {
      dist_opt_val.u = true;
      break;
    }
		case 555:
		{
			dist_opt_val.byread = true;
			break;
		}
		case 'P':
		{
			if(strlen(arg) > PATHLEN) {
        err(errno,"the piped command should not longer than %d", PATHLEN);
        exit(EXIT_FAILURE);
      };
      strcpy(dist_opt_val.pipecmd,arg);
      break;
		}
		case 'M':
		{
			dist_opt_val.metric = atoi(arg) ;	
			break;
		}
		case 'C':
		{
			dist_opt_val.ctm = atoi(arg);
			if(dist_opt_val.ctm >2 || dist_opt_val.ctm < 0 ){
				printf("containment option %d should be 0,1 or2\n",dist_opt_val.ctm);
				exit(1);
			}
			break;
		}
		case 'O':
    {
      dist_opt_val.outfields = atoi(arg) ;
      break;
    }
		case 'f':
		{
			strcpy(dist_opt_val.shared_kmerpath,arg);
			break;
		}
		case 888:
		{
			dist_opt_val.keepco = true;
			break;	
		}
		case 999:
		{
			dist_opt_val.stage2 = true;
      break;
		}
		case 333:
    {
      dist_opt_val.correction = atoi(arg) ;
      break;
    }
		case 777:
		{
			dist_opt_val.keep_shared_kmer = true;
			break;
		}
		case ARGP_KEY_ARGS:
				dist_opt_val.num_remaining_args = state->argc - state->next;
        dist_opt_val.remaining_args  = state->argv + state->next;	
			break;
    case ARGP_KEY_NO_ARGS:
    {
			if(state->argc < 2)
      	{
				printf("\v");
				argp_state_help(state,stdout,ARGP_HELP_SHORT_USAGE);
				printf("\v");
      	argp_state_help(state,stdout,ARGP_HELP_LONG);
      	printf("\v");
				exit(0);
				};
      	return EINVAL;
    }
   break;
   default:
		{
#ifdef _OPENMP
		if(dist_opt_val.p == 0)	
    	dist_opt_val.p =  omp_get_num_procs();
#else
		if(dist_opt_val.p == 0)
    	dist_opt_val.p = 1; 
#endif
    return ARGP_ERR_UNKNOWN;
		}
 }
	//prohibit cases
//	if( (dist_opt_val.metric == Bth) && ( dist_opt_val.num_neigb != 0 || dist_opt_val.mut_dist_max < 1 ) ) err(errno,"when set -M 2, -D and -N should not be set.");

  return 0;
}

static struct argp argp_dist =
{
  opt_dist,
  parse_dist,
  "[-r <reference>] [<query>]",
  doc_dist
};

int cmd_dist(struct argp_state* state)
{
  struct arg_dist dist = { 0, };
  int    argc = state->argc - state->next + 1;
  char** argv = &state->argv[state->next - 1];
  char*  argv0 =  argv[0];

  dist.global = state->input;
  argv[0] = malloc(strlen(state->name) + strlen(" dist") + 1);

  if(!argv[0])
    argp_failure(state, 1, ENOMEM, 0);

  sprintf(argv[0], "%s dist", state->name);
  argp_parse(&argp_dist, argc, argv, ARGP_IN_ORDER, &argc, &dist);
  free(argv[0]);
  argv[0] = argv0;
  state->next += argc - 1;

	return dist_dispatch(&dist_opt_val);
}

const char *mk_dist_rslt_dir (const char *parentdirpath, const char * outdirpath )
{
  struct stat dstat;
  const char *outfullpath = malloc(PATHLEN *sizeof(char));
  sprintf((char *)outfullpath,"%s/%s",parentdirpath,outdirpath);

  if(stat(parentdirpath, &dstat) == 0 && S_ISDIR(dstat.st_mode)){
    if( stat(outfullpath, &dstat) == 0 ){
			errno = EEXIST; 
      err(errno,"%s",outfullpath);
		}
    else{
     // printf("Making results dir: %s\n",outfullpath);
      mkdir(outfullpath,0777);
    }
  }
  else {
    mkdir(parentdirpath,0777);
    mkdir(outfullpath,0777);
  }
  return outfullpath;
};
