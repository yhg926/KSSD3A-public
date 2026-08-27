#ifndef COMMAND_ANI
#define COMMAND_ANI
#include <argp.h>
#include <math.h>
#include <tgmath.h>
#include <assert.h>
#include <err.h>
#include <stdint.h>

#include "global_basic.h"
#include "kssdlib_sort.h"
#include "sketch_rearrange.h"
#include "model_ani.h"
#include "model_refaf_hgb.h"
#include "dna_popcount.h"
#include "kssdlib_hash.h"

typedef struct ani_opt
{
	int fmt;   // print out format: 0:detail; 1, aafD; 2. 1-ani
	double c;  // minimal distance to enrolled sketches
	int p;	   // threads
	bool d;	   // diagnal
	bool v;	   // naive model ?
	bool pair; // pairwise compute
	bool unassembled; // query sketch is unassembled;
	bool unified_metric; // honor -s in unassembled mode instead of forcing naive
	bool ignoreconflict; // ignore reference contexts with conflicting objects
	bool raw_output; // skip calibrated/best ANI; print NULLs in unified detail fields
	bool estimate_coverage; // append abundance-derived coverage columns in detail output
	int e;
	int s; // select metrics;
	int ntop; // report at most top N references for each query
	int ctxcut;
	float afcut;
	bool afcut_set;
	float anicut;
	int sketch_hclen;
	int sketch_holen;
	int sketch_iolen;
	int sketch_drfold;
	int sketch_kmerocrs;
	int sketch_ncap;
	double sketch_npercentile;
	bool sketch_reads_qc;
	bool sketch_abundance;
	bool sketch_asone;
	bool sketch_conflict;
	bool sketch_anno;
	bool sketch_split_mfa;
	bool sketch_coden_ctxobj_pattern;
	char index[PATHLEN];
	char qrydir[PATHLEN];
	char refdir[PATHLEN];
	char qrylist[PATHLEN];
	char reflist[PATHLEN];
	char sketch_pipecmd[PATHLEN];
	char outf[PATHLEN];
	char gl[PATHLEN]; // genome list with selection code
	char model[PATHLEN];
	int num_remaining_args;
	char **remaining_args;
} ani_opt_t;

static inline double ani_report_af_value(double af_qry, double af_ref)
{
	return af_qry > af_ref ? af_qry : af_ref;
}

static inline bool ani_report_af_pass(const ani_opt_t *ani_opt, double af_qry, double af_ref)
{
	return ani_report_af_value(af_qry, af_ref) >= (double)ani_opt->afcut;
}

static inline uint32_t ani_report_af_needed_ctx(const ani_opt_t *ani_opt,
												uint32_t qry_ctx, uint32_t ref_ctx)
{
	const uint32_t denom = qry_ctx < ref_ctx ? qry_ctx : ref_ctx;
	const double need_d = (double)ani_opt->afcut * (double)denom;
	if (need_d <= 0.0)
		return 0;
	if (need_d >= (double)UINT32_MAX)
		return UINT32_MAX;
	return (uint32_t)ceil(need_d);
}

typedef struct
{
	uint32_t num_ctx; // including with confilict obj
	uint32_t num_conflictobj;
} num_ctx_cfltobj_t;

typedef struct
{
	uint64_t arrlen;
	uint32_t num_ctx;
	uint32_t num_conflictobj;
	double numgids_perctx;
	int infile_num;
	num_ctx_cfltobj_t *num_ctx_cfltobj_arr;
} sort_sketch_summary_t;

typedef struct
{
	uint32_t id;
	float ani;
} idani_t;

typedef struct
{
	uint32_t diff_obj;
	uint32_t diff_obj_section;
} obj_section_t;

typedef struct
{
	uint32_t num_ctx;
	uint32_t num_mut2_ctx;
} ctx_mut2_t;

int cmd_ani(struct argp_state *state);
int compute_ani(ani_opt_t *ani_opt);
int mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(ani_opt_t *ani_opt);
int sparse_mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(ani_opt_t *ani_opt);
int stream_ref_sketches_one_qraw_lookup(ani_opt_t *ani_opt);
int stream_ref_sketches_multi_qraw_sortedindex(ani_opt_t *ani_opt);
int compare_u32(const void *a, const void *b);
void sort_arrays(uint64_t *a, uint32_t *b, uint32_t n);
ctxgidobj_t *comb_sortedsketch64_2sortedcomb_ctxgid64obj32(unify_sketch_t *result);
sort_sketch_summary_t *summarize_ctxgidobj_arr(ctxgidobj_t *ctxgidobj_arr, uint64_t *sketch_index, uint32_t arrlen, int infile_num);
void free_sort_sketch_summary(sort_sketch_summary_t *sort_sketch_summary);
// void sorted_ctxgidobj_arr2triangle (ctxgidobj_t* ctxgidobj_arr, sort_sketch_summary_t *sort_sketch_summary);
void sorted_ctxgidobj_arrXcomb_sortedsketch64(unify_sketch_t *qry_result, ctxgidobj_t *ctxgidobj_arr, sort_sketch_summary_t *sort_sketch_summary);
//void comb_sortedsketch64Xcomb_sortedsketch64_sorted_per_q(ani_opt_t *ani_opt);
	void comb_sortedsketch64Xcomb_sortedsketch64_filter_and_sort_survivors(ani_opt_t *ani_opt);
	void comb_manysmall_sortedsketch64Xcomb_fewlarge_sortedsketch64_filter_and_sort_survivors(ani_opt_t *ani_opt);

	void comb_sortedsketch64Xcomb_sortedsketch64(ani_opt_t *ani_opt);
	bool comb_sortedsketch64_indexed_self_full(ani_opt_t *ani_opt);
	bool comb_sortedsketch64_indexed_self_triangle(ani_opt_t *ani_opt);
	void comb_sortedsketch64_self_matrix(ani_opt_t *ani_opt);
	void simple_sortedsketch64Xcomb_sortedsketch64(simple_sketch_t *simple_sketch, infile_tab_t *genomes_infiletab, ani_opt_t *ani_opt);

size_t *find_first_occurrences_AT_ctxgidobj_arr(const uint64_t *a, size_t a_size, const ctxgidobj_t *b, size_t b_size);
// void get_ani_features_from_two_sorted_ctxobj64 (const uint64_t *a, size_t n,  const uint64_t *b, size_t m, ani_features_t* ani_features);
// print functions
void ani_block_print(
	int ref_infile_num, int qry_gid_offset, int this_block_size,
	uint64_t *ref_sketch_index, uint64_t *qry_sketch_index,
	const uint32_t *ref_ctx_count, const uint32_t *qry_ctx_count,
	ctx_mut2_t *ctx, obj_section_t *obj,
	char (*refname)[PATHLEN], char (*qryfname)[PATHLEN],
	char (*refanno)[PATHLEN],
	const infile_meta_t *qry_infile_meta,
	const infile_meta_t *ref_infile_meta,
	uint32_t *num_passid_block, idani_t **sort_idani_block,
	FILE *outfp, ani_opt_t *ani_opt, int matrix_mode);


//fenceposts search
#ifndef KSSD_FENCEPOSTS_H
#define KSSD_FENCEPOSTS_H

#include <stdint.h>
#include <stddef.h>

#ifndef GID_NBITS
#  define GID_NBITS 20  /* you can override in your project config */
#endif

/* Heuristic chooser for k (buckets = 1<<k). Clamp keeps memory/cache sane. */
int kssd_choose_k_fenceposts(size_t b_size, size_t a_size);

/* Build fenceposts F of size (1<<k)+1 for array b (sorted by key(b)=ctxgid>>GID_NBITS).
 * F[t] = first index i such that top-k bits of key(b[i]) >= t. F[2^k] = b_size.
 * Returns 0 on success, -1 on error (bad args).
 */
int kssd_build_fenceposts_ctxgid(const ctxgidobj_t *b, size_t b_size, int k, size_t *F);

/* Batch query: find leftmost index in b for each a[i] with key(a[i])=(a[i]>>nobjbits).
 * Returns malloc'd array of length a_size (caller free), or NULL on OOM.
 * indices[i] = SIZE_MAX if not found.
 */
size_t *kssd_find_first_occurrences_fenceposts(const uint64_t *a, size_t a_size,
                                               const ctxgidobj_t *b, size_t b_size,
                                               const size_t *F, int k, unsigned nobjbits);

#endif /* KSSD_FENCEPOSTS_H */


// inline functions
// 1. 3-way linear model distance (see model_ani.h): static inline double lm3ways_dist_from_features(ani_features_t *features)
// 2. naive distance: for where 3-way linear model is not applicable e.g. unassembled genomes , or Eukaryotic genomes?

#define DIVIDOR (7)
#define NAIVE_DIST_INTERCEPT_RAMP_TAU (0.01)
static inline double get_naive_dist(ani_features_t *features)
{
	if (features->XnY_ctx == 0)
		return 1;
	double ratio = (double)(features->N_diff_obj_section + EPSILON) / (features->N_diff_obj + EPSILON);
	double dist0 = (double)features->N_diff_obj / (features->XnY_ctx + features->N_diff_obj);
	double final_dist = 1 - pow((1 - dist0), ratio);
	if(	final_dist <= 0)
		return 0;

#if NUM_CODENS < 11   //9 or 10
	double scale = 0.1557143;  // 0.1557143 ~ 1/7
	double intercept = 9.961e-5;
#else 
	double scale = 0.1544286;
	double intercept = 7.133e-4;
#endif
	double intercept_weight = final_dist / NAIVE_DIST_INTERCEPT_RAMP_TAU;
	if (intercept_weight > 1.0)
		intercept_weight = 1.0;
	double predict_dist = final_dist * scale + intercept * intercept_weight;
	return predict_dist;
}

static inline double get_p_dist(ani_features_t *features)
{
	if (features->XnY_ctx == 0 || Bitslen.obj == 0)
		return 1;
	return (double)features->N_diff_obj_section /
		   ((double)features->XnY_ctx * ((double)Bitslen.obj / 2.0));
}

// 3. generic distance function from 1. or 2.
typedef double (*get_generic_dist_from_features_fn)(ani_features_t *features);
extern get_generic_dist_from_features_fn get_generic_dist_from_features;
// called by comb_sortedsketch64Xcomb_sortedsketch64()
// cautions: a and b must be obj conflition removed arrays.
static inline void get_ani_features_from_two_sorted_ctxobj64(const uint64_t *a, size_t n,
															 const uint64_t *b, size_t m, ani_features_t *ani_features)
{
	uint8_t nobjbits = Bitslen.obj;
	uint64_t objmask = (1UL << nobjbits) - 1;
	size_t i = 0, j = 0;
	memset(ani_features, 0, sizeof(ani_features_t));

	while (i < n && j < m)
	{

		if (a[i] >> nobjbits == b[j] >> nobjbits)
		{
			ani_features->XnY_ctx++;
			uint32_t has_diff_obj = (uint32_t)(a[i] & objmask) ^ (b[j] & objmask);
			if (has_diff_obj)
			{
				ani_features->N_diff_obj++;
				/* old count method
				int num_diff_obj_section = 0;
				for (int k = 0; k < nobjbits / 2; k++)
				{
					if (has_diff_obj & (3U << (2 * k)))
						num_diff_obj_section++;
				}
				*/
				int num_diff_obj_section = dna_popcount(has_diff_obj);
				ani_features->N_diff_obj_section += num_diff_obj_section;
				if (num_diff_obj_section > 1)
					ani_features->N_mut2_ctx++;
			}
			i++;
			j++;
		}
		else if (a[i] < b[j])
		{ // else if (a[i] >> nobjbits < b[j] >> nobjbits)
			i++;
		}
		else
		{
			j++;
		}
	}
}
static inline int min_diff_sections_ctxobj64_runs(const uint64_t *a, size_t a_begin, size_t a_end,
                                                  const uint64_t *b, size_t b_begin, size_t b_end,
                                                  uint64_t objmask)
{
    int min_diff_sections = NUM_CODENS + 1;

    for (size_t ai = a_begin; ai < a_end; ++ai) {
        const uint32_t objA = (uint32_t)(a[ai] & objmask);
        for (size_t bi = b_begin; bi < b_end; ++bi) {
            const uint32_t diff = objA ^ (uint32_t)(b[bi] & objmask);
            if (diff == 0) return 0;
            const int d = dna_popcount(diff);
            if (d < min_diff_sections) min_diff_sections = d;
        }
    }
    return min_diff_sections;
}

// a and b may both contain multiple objects for the same context.
// Both arrays must be sorted ascending by ctx-object key.
static inline void get_ani_features_ctx_min_over_conflicts_both_filtered(const uint64_t *a, size_t n,const uint64_t *b, size_t m, bool ignore_b_conflict, ani_features_t *ani_features)
{
    const uint8_t  nobjbits = Bitslen.obj;
    const uint64_t objmask  = ((1ULL << nobjbits) - 1ULL);
    memset(ani_features, 0, sizeof *ani_features);

    size_t i = 0, j = 0;

    while (i < n && j < m) {
        const uint64_t ctxA = a[i] >> nobjbits;
        const uint64_t ctxB = b[j] >> nobjbits;

        if (ctxA < ctxB) {
            do { ++i; } while (i < n && (a[i] >> nobjbits) == ctxA);
            continue;
        }
        if (ctxB < ctxA) {
            do { ++j; } while (j < m && (b[j] >> nobjbits) == ctxB);
            continue;
        }

        const size_t a_begin = i;
        do { ++i; } while (i < n && (a[i] >> nobjbits) == ctxA);
        const size_t a_end = i;

        const size_t b_begin = j;
        do { ++j; } while (j < m && (b[j] >> nobjbits) == ctxB);
        const size_t b_end = j;

        if (ignore_b_conflict && b_end - b_begin > 1)
            continue;

        ani_features->XnY_ctx++;
        const int min_diff_sections = min_diff_sections_ctxobj64_runs(a, a_begin, a_end, b, b_begin, b_end, objmask);

        if (min_diff_sections > 0) {
            ani_features->N_diff_obj++;
            ani_features->N_diff_obj_section += min_diff_sections;
            if (min_diff_sections > 1) {
                ani_features->N_mut2_ctx++;
            }
        }
    }
}

static inline void get_ani_features_ctx_min_over_conflicts_both(const uint64_t *a, size_t n,const uint64_t *b, size_t m, ani_features_t *ani_features)
{
    get_ani_features_ctx_min_over_conflicts_both_filtered(a, n, b, m, false, ani_features);
}

// Compatibility wrapper for older call sites; the both-side implementation is
// identical for conflict-free b, and also works when b keeps conflicts.
static inline void get_ani_features_ctx_min_over_conflicts_a_only(const uint64_t *a, size_t n,const uint64_t *b, size_t m, ani_features_t *ani_features)
{
    get_ani_features_ctx_min_over_conflicts_both(a, n, b, m, ani_features);
}

// Descending comparator by ANI (highest first)
static int compare_idani_desc(const void *a, const void *b)
{
	const idani_t *itemA = (const idani_t *)a;
	const idani_t *itemB = (const idani_t *)b;
	return (itemA->ani < itemB->ani) - (itemA->ani > itemB->ani);
}

// ---------------- in-place Top-N selection (min-heap on first N) ----------------
static inline void heap_sift_down_min_idani(idani_t *h, int n, int i)
{
    for (;;) {
        int l = i * 2 + 1;
        int r = l + 1;
        int s = i;
        if (l < n && h[l].ani < h[s].ani) s = l;
        if (r < n && h[r].ani < h[s].ani) s = r;
        if (s == i) break;
        idani_t tmp = h[i]; h[i] = h[s]; h[s] = tmp;
        i = s;
    }
}

static inline void heapify_min_idani(idani_t *h, int n)
{
    for (int i = (n >> 1) - 1; i >= 0; --i)
        heap_sift_down_min_idani(h, n, i);
}

// Keep top N elements in arr[0..N) (unordered tail). If sort_top!=0, sort arr[0..N) descending.
static inline void partial_topN_inplace_idani(idani_t *arr, int m, int N, int sort_top)
{
    if (N <= 0 || m <= 1) return;

    if (N >= m) {
        qsort(arr, m, sizeof(*arr), compare_idani_desc);
        return;
    }

    // Build min-heap in the first N slots
    heapify_min_idani(arr, N);

    // Stream remaining items, keep only top N
    for (int i = N; i < m; ++i) {
        if (arr[i].ani <= arr[0].ani) continue;

        // swap candidate into heap root; move old root to tail
        idani_t tmp = arr[i];
        arr[i] = arr[0];
        arr[0] = tmp;

        heap_sift_down_min_idani(arr, N, 0);
    }

    if (sort_top)
        qsort(arr, N, sizeof(*arr), compare_idani_desc);
}

// ---------------- recommended wrapper: choose qsort vs top-N ----------------
// After this call, the best min(topN,m) items are in arr[0..out_n), sorted descending.
// Returns out_n.
static inline int idani_select_topN_sorted(idani_t *arr, int m, int topN)
{
    if (m <= 0 || topN <= 0 ) return 0;
    if (topN > m) topN = m;
    // Heuristic: small m or topN not much smaller than m => full sort is fine
    if (m <= 128 || topN >= (m / 2)) {
        qsort(arr, m, sizeof(*arr), compare_idani_desc);
        return topN;
    }

    partial_topN_inplace_idani(arr, m, topN, /*sort_top=*/1);
    return topN;
}

#define MCTX(L, X, Y) (ctx[(size_t)((L) * (X) + (Y))])
#define MOBJ(L, X, Y) (obj[(size_t)((L) * (X) + (Y))])
static inline void count_ctx_obj_frm_comb_sketch_section(ctx_mut2_t *ctx, obj_section_t *obj, ctxgidobj_t *ctxgidobj_arr, size_t ref_sksize, int ref_gnum, const uint32_t *ref_ctx_count, const uint32_t *qry_ctx_count_block, int section_gnum, uint64_t *section_sk, uint64_t *section_skidx, uint32_t *num_passid_block, idani_t **sort_idani_block, ani_opt_t *ani_opt)
{
	uint32_t nobjbits = Bitslen.obj;
	uint64_t gidmask = UINT64_MAX >> (64 - GID_NBITS), objmask = (1UL << nobjbits) - 1;
	uint32_t ctxcut = (uint32_t)ani_opt->ctxcut;
	float anicut = ani_opt->anicut;
	int ntop = ani_opt->ntop > 0 ? ani_opt->ntop : ref_gnum;
	int fmt = ani_opt->fmt; 
	// for fencepost methods
	int k = 17;
	size_t buckets = (size_t)1u << k;
	size_t *F= (size_t*) malloc((buckets + 1) * sizeof(*F) );
	if(kssd_build_fenceposts_ctxgid(ctxgidobj_arr, ref_sksize,k,F) != 0) {
		free(F);
		return;
	}

#pragma omp parallel for num_threads(ani_opt->p) schedule(guided)
	for (int i = 0; i < section_gnum; i++)
	{
		size_t a_begin = (size_t)(section_skidx[i] - section_skidx[0]);
		size_t a_size = (size_t)(section_skidx[i + 1] - section_skidx[i]);
		if (a_size == 0) {
			if (fmt == 0)
				num_passid_block[i] = 0;
			continue;
		}
		uint64_t *a = section_sk + a_begin;
		//fencepost
		
		size_t *idx= kssd_find_first_occurrences_fenceposts(a, a_size, ctxgidobj_arr, ref_sksize,F,k,nobjbits);
		if(!idx) {
			if (fmt == 0)
				num_passid_block[i] = 0;
			continue;
		}

		//size_t *idx = find_first_occurrences_AT_ctxgidobj_arr(a, a_size, ctxgidobj_arr, ref_sksize);
		
		// record the minimum number of different objects when a has same context with different objects
        //uint8_t *min_obj_sec_confict_a = malloc(ref_gnum);
		int s = 1;
		for (int j = 0; j < a_size; j+=s)
		{
			// consecuted s same context with different objects from array a are anlysised in one batch   
			for(s = 1;j+s < a_size && idx[j] == idx[j+s];s++);
			
			if (idx[j] == SIZE_MAX) continue;
			// Skip when no findings, or conflict objects (adjacent elements with the same context)
			// if ((j > 0 && (a[j] >> Bitslen.obj) == (a[j - 1] >> Bitslen.obj)) ||
			//	(j < a_size - 1 && (a[j] >> Bitslen.obj) == (a[j + 1] >> Bitslen.obj)))
			//	continue;
			size_t d = idx[j];
			while (d < ref_sksize && (ctxgidobj_arr[d].ctxgid >> GID_NBITS) == (a[j] >> nobjbits))
			{
				const uint64_t ctxgid = ctxgidobj_arr[d].ctxgid;
				uint32_t gid = ctxgid & gidmask;
				const size_t ref_run_begin = d;
				do { ++d; } while (d < ref_sksize && ctxgidobj_arr[d].ctxgid == ctxgid);
				const size_t ref_run_end = d;
				if (ani_opt->ignoreconflict && ref_run_end - ref_run_begin > 1)
					continue;
				MCTX(ref_gnum, i, gid).num_ctx++;
				int min_diff_obj_section = NUM_CODENS + 1; 

				for(int n = 0; n < s; n++){
					uint32_t obj_q = (uint32_t)(a[j+n] & objmask);
					for (size_t r = ref_run_begin; r < ref_run_end; r++) {
						uint32_t has_diff_obj_n = obj_q ^ ctxgidobj_arr[r].obj;
						if (has_diff_obj_n == 0){
							min_diff_obj_section = 0;
							break;
						}
						int diff_obj_section_n = dna_popcount(has_diff_obj_n);
						if (diff_obj_section_n < min_diff_obj_section)
							min_diff_obj_section = diff_obj_section_n;
					}
					if (min_diff_obj_section == 0) break;
				}

				if(min_diff_obj_section){
					MOBJ(ref_gnum, i, gid).diff_obj++;
					MOBJ(ref_gnum, i, gid).diff_obj_section += min_diff_obj_section;
					if (min_diff_obj_section > 1)
						MCTX(ref_gnum, i, gid).num_mut2_ctx++;

				}
			}
		}
		free(idx);
		if (fmt == 0)
		{ // print details  >0:matrix no need to sort
			// sorting ref gid by ani descendingly
			num_passid_block[i] = 0;
			for (int j = 0; j < ref_gnum; j++)
			{
				uint32_t shared_ctx = MCTX(ref_gnum, i, j).num_ctx;
				if (shared_ctx < ctxcut)
					continue;
				if (ref_ctx_count[j] == 0 || qry_ctx_count_block[i] == 0 ||
					!ani_report_af_pass(ani_opt,
										(double)shared_ctx / (double)qry_ctx_count_block[i],
										(double)shared_ctx / (double)ref_ctx_count[j]))
					continue;
				float dist = (float)MOBJ(ref_gnum, i, j).diff_obj / MCTX(ref_gnum, i, j).num_ctx / DIVIDOR;
				float ani = 1 - dist;
				if (ani < anicut)
					continue;
				sort_idani_block[i][num_passid_block[i]].id = j;
				sort_idani_block[i][num_passid_block[i]].ani = ani;
				num_passid_block[i]++;
			}
			int out_n = idani_select_topN_sorted(sort_idani_block[i], (int)num_passid_block[i], ntop);
			num_passid_block[i] = (uint32_t)out_n;  // downstream prints only topN now
			//qsort(sort_idani_block[i], num_passid_block[i], sizeof(idani_t), compare_idani_desc);
			//num_passid_block[i] = num_passid_block[i] < ntop ? num_passid_block[i] : ntop;
		}
	}
	free(F);
}

static inline void count_ctx_obj_frm_comb_sketch_section_lower(
	ctx_mut2_t *ctx, obj_section_t *obj,
	ctxgidobj_t *ctxgidobj_arr, size_t ref_sksize, int ref_gnum,
	int qry_gid_offset, int section_gnum,
	uint64_t *section_sk, uint64_t *section_skidx,
	ani_opt_t *ani_opt)
{
	const uint32_t nobjbits = Bitslen.obj;
	const uint64_t gidmask = UINT64_MAX >> (64 - GID_NBITS);
	const uint64_t objmask = (1UL << nobjbits) - 1;
	const int k = 17;
	const size_t buckets = (size_t)1u << k;
	size_t *F = malloc((buckets + 1) * sizeof(*F));
	if (!F)
		err(EXIT_FAILURE, "%s(): malloc fenceposts", __func__);
	if (kssd_build_fenceposts_ctxgid(ctxgidobj_arr, ref_sksize, k, F) != 0) {
		free(F);
		errx(EXIT_FAILURE, "%s(): invalid sorted context index", __func__);
	}

#pragma omp parallel for num_threads(ani_opt->p) schedule(guided)
	for (int i = 0; i < section_gnum; ++i) {
		const int qry_gid = qry_gid_offset + i;
		const size_t a_begin = (size_t)(section_skidx[i] - section_skidx[0]);
		const size_t a_size = (size_t)(section_skidx[i + 1] - section_skidx[i]);
		if (a_size == 0)
			continue;
		uint64_t *a = section_sk + a_begin;
		size_t *idx = kssd_find_first_occurrences_fenceposts(
			a, a_size, ctxgidobj_arr, ref_sksize, F, k, nobjbits);
		if (!idx)
			continue;

		int s = 1;
		for (size_t j = 0; j < a_size; j += (size_t)s) {
			for (s = 1; j + (size_t)s < a_size && idx[j] == idx[j + (size_t)s]; ++s)
				;
			if (idx[j] == SIZE_MAX)
				continue;
			size_t d = idx[j];
			const uint64_t qry_ctx = a[j] >> nobjbits;
			while (d < ref_sksize &&
				   (ctxgidobj_arr[d].ctxgid >> GID_NBITS) == qry_ctx) {
				const uint64_t ctxgid = ctxgidobj_arr[d].ctxgid;
				const uint32_t gid = (uint32_t)(ctxgid & gidmask);
				const size_t ref_run_begin = d;
				do {
					++d;
				} while (d < ref_sksize && ctxgidobj_arr[d].ctxgid == ctxgid);
				const size_t ref_run_end = d;
				if ((int)gid >= qry_gid) {
					if ((int)gid > qry_gid)
						break;
					continue;
				}
				if (ani_opt->ignoreconflict && ref_run_end - ref_run_begin > 1)
					continue;

				MCTX(ref_gnum, i, gid).num_ctx++;
				int min_diff_obj_section = NUM_CODENS + 1;
				for (int n = 0; n < s; ++n) {
					const uint32_t obj_q = (uint32_t)(a[j + (size_t)n] & objmask);
					for (size_t r = ref_run_begin; r < ref_run_end; ++r) {
						const uint32_t has_diff_obj_n = obj_q ^ ctxgidobj_arr[r].obj;
						if (has_diff_obj_n == 0) {
							min_diff_obj_section = 0;
							break;
						}
						const int diff_obj_section_n = dna_popcount(has_diff_obj_n);
						if (diff_obj_section_n < min_diff_obj_section)
							min_diff_obj_section = diff_obj_section_n;
					}
					if (min_diff_obj_section == 0)
						break;
				}
				if (min_diff_obj_section) {
					MOBJ(ref_gnum, i, gid).diff_obj++;
					MOBJ(ref_gnum, i, gid).diff_obj_section += min_diff_obj_section;
					if (min_diff_obj_section > 1)
						MCTX(ref_gnum, i, gid).num_mut2_ctx++;
				}
			}
		}
		free(idx);
	}
	free(F);
}




#endif
