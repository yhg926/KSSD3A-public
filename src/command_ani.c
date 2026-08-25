#include "command_ani.h"
#include "command_matrix.h"
#include "global_basic.h"
#include "kssdlib_sort.h"
#include "sketch_rearrange.h"
// #include "command_sketch.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <err.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <libgen.h>
#include <dirent.h>
#include <omp.h>
#include <stdatomic.h>
#include <ctype.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include "../klib/kstring.h" // from klib
#include "../klib/khash.h"
#include "../klib/kvec.h"      // klib: dynamic arrays
// #include "../klib/khashl.h"
#define GID_NBITS 20 // 2^20, 1M
#define CONFLICT_OBJ UINT32_MAX
// pulic vars
const char gid_obj_prefix[] = "gidobj", ctx_idx_prefix[] = "ctx.index";
extern const char sorted_comb_ctxgid64obj32[];
extern double C9O7_98[6], C9O7_96[6];
size_t file_size;

const char unified_detail_header[] = "Qry\tRef\tANI\tDistance\tConfidence\tSelected_metric\tXnY_ctx\tQry_align_fraction\tblastn_Qry_align_fraction\tRef_align_fraction\tblastn_Ref_align_fraction\tN_diff_obj\tN_diff_obj_section\tN_mut2_ctx\tRef_annotation";
#define ANI_SELECTED_METRIC_COUNT 9
const char select_metrics_header[ANI_SELECTED_METRIC_COUNT][20] = {
	"BestDist", "RecalDist", "CtxMoE", "Naive",
	"MashD", "AafD", "MashD_if_far", "AafD_if_far", "p_dist"
};

typedef enum {
	ANI_CONF_LOW = 0,
	ANI_CONF_MEDIUM = 1,
	ANI_CONF_HIGH = 2
} ani_confidence_t;

static inline const char *ani_confidence_label(unsigned char confidence)
{
	switch ((ani_confidence_t)confidence) {
	case ANI_CONF_HIGH:
		return "high";
	case ANI_CONF_MEDIUM:
		return "medium";
	case ANI_CONF_LOW:
	default:
		return "low";
	}
}

static inline unsigned char ani_confidence_from_values(double raw_ani, double calibrated_ani, double ref_af)
{
	const double drop = raw_ani - calibrated_ani;
	if (ref_af < 0.5 || (raw_ani >= 0.95 && drop > 0.02))
		return ANI_CONF_LOW;
	if (fabs(drop) > 0.01)
		return ANI_CONF_MEDIUM;
	return ANI_CONF_HIGH;
}

static inline const infile_meta_t *infile_meta_at(const unify_sketch_t *sketch,
															uint32_t idx)
{
	if (!sketch || !sketch->infile_meta || idx >= (uint32_t)sketch->infile_num)
		return NULL;
	return &sketch->infile_meta[idx];
}

static inline bool infile_meta_complete_like_assembly(const infile_meta_t *s)
{
	if (!s || s->meta_fmt_version != KSSD_INFILE_META_VERSION ||
		s->infile_fmt != KSSD_INFILE_FMT_FASTA ||
		s->total_length_bp == 0 || s->median_length_bp == 0)
		return false;
	return s->total_length_bp >= 3500000ULL &&
		   s->record_count < 5 &&
		   s->median_length_bp > 100000U &&
		   (double)s->total_length_bp / (double)s->median_length_bp < 5.0;
}

static inline bool has_assembly_meta_record(const infile_meta_t *s)
{
	return s && s->meta_fmt_version == KSSD_INFILE_META_VERSION &&
		   s->infile_fmt == KSSD_INFILE_FMT_FASTA &&
		   s->total_length_bp > 0 && s->record_count > 0;
}

static inline double bounded_ani(double value)
{
	if (!isfinite(value))
		return NAN;
	if (value < 0.0)
		return 0.0;
	if (value > 1.0)
		return 1.0;
	return value;
}

static inline double bounded_align_fraction(double value)
{
	if (!isfinite(value))
		return value;
	if (value < 0.0)
		return 0.0;
	if (value > 1.0)
		return 1.0;
	return value;
}

static inline double ctx_k_for_ani(void)
{
	return Bitslen.ctx > 0 ? (double)Bitslen.ctx / 2.0 : 22.0;
}

static inline double mash_ani_from_counts(double overlap, double qry_ctx, double ref_ctx)
{
	if (overlap <= 0.0 || qry_ctx <= 0.0 || ref_ctx <= 0.0)
		return 0.0;
	const double denom = qry_ctx + ref_ctx - overlap;
	if (denom <= 0.0)
		return 0.0;
	const double jaccard = overlap / denom;
	if (jaccard <= 0.0)
		return 0.0;
	return bounded_ani(1.0 + log(2.0 * jaccard / (1.0 + jaccard)) / ctx_k_for_ani());
}

static inline double aaf_ani_from_counts(double overlap, double qry_ctx, double ref_ctx)
{
	if (overlap <= 0.0 || qry_ctx <= 0.0 || ref_ctx <= 0.0)
		return 0.0;
	const double min_ctx = qry_ctx < ref_ctx ? qry_ctx : ref_ctx;
	if (min_ctx <= 0.0)
		return 0.0;
	const double containment = overlap / min_ctx;
	if (containment <= 0.0)
		return 0.0;
	return bounded_ani(1.0 + log(containment) / ctx_k_for_ani());
}

static inline void print_ani_detail_header(FILE *outfp, const ani_opt_t *ani_opt, bool include_selected_metric)
{
	(void)ani_opt;
	(void)include_selected_metric;
	fprintf(outfp, "%s\n", unified_detail_header);
}

static inline const char *annotation_or_na(const char *annotation)
{
	return annotation && annotation[0] ? annotation : "NA";
}

static inline const char *annotation_at(char (*annotations)[PATHLEN], uint32_t idx)
{
	return annotation_or_na(annotations ? annotations[idx] : NULL);
}

static inline const char *unify_annotation_at(const unify_sketch_t *sketch, uint32_t idx)
{
	return annotation_at(sketch ? sketch->annotation : NULL, idx);
}

static char (*read_optional_sketch_annotations(const char *sketch_dir, int infile_num))[PATHLEN]
{
	if (infile_num <= 0 || !file_exists_in_folder(sketch_dir, sketch_anno_stat))
		return NULL;
	size_t anno_file_size = 0;
	char *anno_path = test_get_fullpath(sketch_dir, sketch_anno_stat);
	char (*annotations)[PATHLEN] =
		read_from_file(anno_path, &anno_file_size);
	free(anno_path);
	const size_t expected_size = (size_t)infile_num * PATHLEN;
	if (anno_file_size != expected_size)
		err(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
			__func__, sketch_dir, sketch_anno_stat, anno_file_size, expected_size);
	return annotations;
}

static infile_meta_t *read_optional_sketch_infile_meta_stats(const char *sketch_dir, int infile_num)
{
	if (infile_num <= 0 || !file_exists_in_folder(sketch_dir, sketch_infile_meta_stat))
		return NULL;
	size_t meta_file_size = 0;
	char *meta_path = test_get_fullpath(sketch_dir, sketch_infile_meta_stat);
	infile_meta_t *stats = read_from_file(meta_path, &meta_file_size);
	free(meta_path);
	const size_t expected_size = (size_t)infile_num * sizeof(stats[0]);
	if (meta_file_size != expected_size)
		err(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
			__func__, sketch_dir, sketch_infile_meta_stat, meta_file_size, expected_size);
	return stats;
}

static inline bool ani_best_guard_enabled(const ani_opt_t *ani_opt)
{
	return ani_opt && ani_opt->fmt == 0 && !ani_opt->unassembled && !ani_opt->v;
}

static void load_infile_meta_for_best_guard(unify_sketch_t *sketch, const char *sketch_dir,
											   const ani_opt_t *ani_opt)
{
	if (!ani_best_guard_enabled(ani_opt) || !sketch || sketch->infile_meta)
		return;
	sketch->infile_meta = read_optional_sketch_infile_meta_stats(sketch_dir, sketch->infile_num);
}

static unsigned ani_query_parse_flags(const ani_opt_t *ani_opt)
{
	return ani_best_guard_enabled(ani_opt) ? SKETCH_PARSE_INFILE_META : SKETCH_PARSE_NONE;
}

static unsigned ani_ref_parse_flags(const ani_opt_t *ani_opt)
{
	unsigned flags = SKETCH_PARSE_NONE;
	if (ani_opt && ani_opt->fmt == 0)
		flags |= SKETCH_PARSE_ANNOTATION;
	if (ani_best_guard_enabled(ani_opt))
		flags |= SKETCH_PARSE_INFILE_META;
	return flags;
}

static bool force_ref_index_requested(void)
{
	const char *env = getenv("KSSD3A_FORCE_REF_INDEX");
	return env && env[0] != '\0' && strcmp(env, "0") != 0;
}

static void *read_mmap_ro_file(const char *file_path, size_t *file_size)
{
	int fd = open(file_path, O_RDONLY);
	if (fd == -1)
		err(EXIT_FAILURE, "%s(): open %s failed", __func__, file_path);

	struct stat st;
	if (fstat(fd, &st) != 0)
	{
		close(fd);
		err(EXIT_FAILURE, "%s(): fstat %s failed", __func__, file_path);
	}
	if (st.st_size <= 0)
	{
		close(fd);
		errx(EXIT_FAILURE, "%s(): file is empty: %s", __func__, file_path);
	}

	void *buffer = mmap(NULL, (size_t)st.st_size, PROT_READ, MAP_SHARED, fd, 0);
	if (buffer == MAP_FAILED)
	{
		close(fd);
		err(EXIT_FAILURE, "%s(): mmap %s failed", __func__, file_path);
	}
	close(fd);
#ifdef MADV_WILLNEED
	(void)madvise(buffer, (size_t)st.st_size, MADV_WILLNEED);
#endif
	if (file_size)
		*file_size = (size_t)st.st_size;
	return buffer;
}

static ctxgidobj_t *read_reference_sorted_index(const char *file_path, size_t *file_size,
												bool *is_mmap)
{
	*is_mmap = force_ref_index_requested();
	if (*is_mmap)
		return read_mmap_ro_file(file_path, file_size);
	return read_from_file(file_path, file_size);
}

static void free_reference_sorted_index(ctxgidobj_t *index, size_t file_size, bool is_mmap)
{
	if (is_mmap)
	{
		if (munmap(index, file_size) == -1)
			err(EXIT_FAILURE, "%s(): munmap failed", __func__);
		return;
	}
	free_read_from_file(index, file_size);
}

static uint32_t count_ctx_runs_sorted_ctxobj64_local(const uint64_t *a, size_t n, bool ignoreconflict)
{
	const uint8_t nobjbits = Bitslen.obj;
	uint32_t runs = 0;
	for (size_t i = 0; i < n; ) {
		const uint64_t ctx = a[i] >> nobjbits;
		const size_t begin = i;
		do { ++i; } while (i < n && (a[i] >> nobjbits) == ctx);
		if (!ignoreconflict || i - begin == 1)
			runs++;
	}
	return runs;
}

static uint32_t *ref_ctx_counts_from_ctxgidobj(const ctxgidobj_t *arr, size_t arrlen, int ref_infile_num,
											   const uint64_t *ref_sketch_index, bool ignoreconflict)
{
	uint32_t *counts = calloc((size_t)ref_infile_num, sizeof(*counts));
	if (!counts)
		err(EXIT_FAILURE, "%s(): calloc counts", __func__);

	for (int gid = 0; gid < ref_infile_num; gid++)
		counts[gid] = (uint32_t)(ref_sketch_index[gid + 1] - ref_sketch_index[gid]);

	for (size_t i = 0; i < arrlen; ) {
		const uint64_t ctxgid = arr[i].ctxgid;
		const uint32_t gid = (uint32_t)(ctxgid & ((1ULL << GID_NBITS) - 1ULL));
		size_t run_len = 0;
		do { ++i; ++run_len; } while (i < arrlen && arr[i].ctxgid == ctxgid);
		if (run_len > 1)
			counts[gid] -= (uint32_t)(ignoreconflict ? run_len : run_len - 1);
	}
	return counts;
}

static void fill_ctx_counts_for_query_block(uint32_t *qry_ctx_count, int offset_gid, int this_block_size,
											const uint64_t *qry_sketch_index, const uint64_t *tmp_ctxobj,
											bool qry_conflict)
{
	const uint64_t block_base = qry_sketch_index[offset_gid];
	for (int i = 0; i < this_block_size; i++) {
		const int gid = offset_gid + i;
		const uint64_t begin = qry_sketch_index[gid] - block_base;
		const uint64_t end = qry_sketch_index[gid + 1] - block_base;
		const size_t len = (size_t)(end - begin);
		if (len == 0) {
			qry_ctx_count[gid] = 0;
			continue;
		}
		qry_ctx_count[gid] = qry_conflict
			? count_ctx_runs_sorted_ctxobj64_local(tmp_ctxobj + begin, len, false)
			: (uint32_t)len;
	}
}

#define BLOCK_SIZE (4096) // #of qry genomes per batch, for mem_eff handling
int mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(ani_opt_t *ani_opt)
{
	
	dim_sketch_stat_t *ref_dim_sketch_stat = read_from_file(test_get_fullpath(ani_opt->refdir, sketch_stat), &file_size);
	int ref_infile_num = ref_dim_sketch_stat->infile_num;
	// read index
	size_t ctxgidobj_arr_fsize;
	uint64_t *ref_sketch_index = read_from_file(test_get_fullpath(ani_opt->refdir, idx_sketch_suffix), &file_size);
	assert(file_size == (ref_infile_num + 1) * sizeof(ref_sketch_index[0]));
	size_t ref_sketch_size = ref_sketch_index[ref_infile_num];
	char *sorted_index_path = test_get_fullpath(ani_opt->refdir, sorted_comb_ctxgid64obj32);
	bool sorted_index_is_mmap = false;
	ctxgidobj_t *sortedcomb_ctxgid64obj32 = read_reference_sorted_index(sorted_index_path, &ctxgidobj_arr_fsize, &sorted_index_is_mmap);
	free(sorted_index_path);
	assert(ctxgidobj_arr_fsize == ref_sketch_size * sizeof(sortedcomb_ctxgid64obj32[0]));
	uint32_t *ref_ctx_count = ref_ctx_counts_from_ctxgidobj(sortedcomb_ctxgid64obj32, ref_sketch_size, ref_infile_num, ref_sketch_index, ani_opt->ignoreconflict);

	dim_sketch_stat_t *qry_dim_sketch_stat = read_from_file(test_get_fullpath(ani_opt->qrydir, sketch_stat), &file_size);
	int qry_infile_num = qry_dim_sketch_stat->infile_num;
	assert(qry_dim_sketch_stat->hash_id == ref_dim_sketch_stat->hash_id);
	uint64_t *qry_sketch_index = read_from_file(test_get_fullpath(ani_opt->qrydir, idx_sketch_suffix), &file_size);
	size_t qry_sketch_size = qry_sketch_index[qry_infile_num];
	uint32_t *qry_ctx_count = calloc((size_t)qry_infile_num, sizeof(*qry_ctx_count));
	if (!qry_ctx_count)
		err(EXIT_FAILURE, "%s(): calloc qry_ctx_count", __func__);

	int block_size = qry_infile_num < BLOCK_SIZE ? qry_infile_num : BLOCK_SIZE;
	int offset_gid = 0;
	const char *qry_comb_sketch_path = test_get_fullpath(ani_opt->qrydir, combined_sketch_suffix);
	FILE *fp = fopen(qry_comb_sketch_path, "rb");
	if (fp == NULL)
		err(errno, "%s", qry_comb_sketch_path);
	uint64_t *tmp_ctxobj = NULL;
	size_t tmp_ctxobj_capacity = 0;
	// uint32_t *ctx = malloc(ref_infile_num * block_size * sizeof(uint32_t));
	// uint32_t *obj = malloc( ref_infile_num * block_size * sizeof(uint32_t))  ;
	ctx_mut2_t *ctx = malloc(ref_infile_num * block_size * sizeof(ctx_mut2_t));
	obj_section_t *obj = malloc(ref_infile_num * block_size * sizeof(obj_section_t));
	// for order id by descending ani
	uint32_t *num_passid_block = malloc(block_size * sizeof(uint32_t));
	idani_t **sort_idani_block = malloc(block_size * sizeof(idani_t *));
	for (int i = 0; i < block_size; i++)
		sort_idani_block[i] = malloc(ref_infile_num * sizeof(idani_t));

	char (*refname)[PATHLEN] = (char (*)[PATHLEN])(ref_dim_sketch_stat + 1);
	char (*qryname)[PATHLEN] = (char (*)[PATHLEN])(qry_dim_sketch_stat + 1);
	char (*refanno)[PATHLEN] = read_optional_sketch_annotations(ani_opt->refdir, ref_infile_num);
	const bool enable_best_guard = ani_best_guard_enabled(ani_opt);
	infile_meta_t *qry_infile_meta =
		enable_best_guard ? read_optional_sketch_infile_meta_stats(ani_opt->qrydir, qry_infile_num) : NULL;
	infile_meta_t *ref_infile_meta =
		enable_best_guard ? read_optional_sketch_infile_meta_stats(ani_opt->refdir, ref_infile_num) : NULL;

	FILE *outfp = ani_opt->outf[0] == '\0' ? stdout : fopen(ani_opt->outf, "w");
	if (outfp == NULL)
		err(errno, "%s", ani_opt->outf);

	/* load model
	if (ani_opt->model[0] == '\0')
		err(EXIT_FAILURE, "%s(): need specify model file using -M ", __func__);
	init_model(ani_opt->model); // f8C9O7_model.xgb
	*/

	// printf header
	if (ani_opt->fmt)
	{ // matrix format
		for (int i = 0; i < ref_infile_num; i++)
			fprintf(outfp, "\t%s", refname[i]);
		fprintf(outfp, "\n");
	}
	else
		print_ani_detail_header(outfp, ani_opt, true);

	while (offset_gid < qry_infile_num)
	{

		int this_block_size = qry_infile_num - offset_gid;
		if (this_block_size > block_size)
			this_block_size = block_size;
		size_t this_sketch_size = (size_t)(qry_sketch_index[offset_gid + this_block_size] - qry_sketch_index[offset_gid]);
		if (this_sketch_size > tmp_ctxobj_capacity)
		{
			uint64_t *new_tmp_ctxobj = realloc(tmp_ctxobj, this_sketch_size * sizeof(uint64_t));
			if (!new_tmp_ctxobj)
				err(EXIT_FAILURE, "%s(): realloc tmp_ctxobj", __func__);
			tmp_ctxobj = new_tmp_ctxobj;
			tmp_ctxobj_capacity = this_sketch_size;
		}
		size_t read_sketch_size = fread(tmp_ctxobj, sizeof(uint64_t), this_sketch_size, fp);
		uint64_t *this_sketch_index = qry_sketch_index + offset_gid;
		assert(this_sketch_size == read_sketch_size);
		fill_ctx_counts_for_query_block(qry_ctx_count, offset_gid, this_block_size, qry_sketch_index, tmp_ctxobj, qry_dim_sketch_stat->conflict);

		memset(ctx, 0, ref_infile_num * block_size * sizeof(ctx_mut2_t));
		memset(obj, 0, ref_infile_num * block_size * sizeof(obj_section_t)); // memset(obj,0,ref_infile_num * block_size * sizeof(uint32_t));
		count_ctx_obj_frm_comb_sketch_section(ctx, obj, sortedcomb_ctxgid64obj32, ref_sketch_size, ref_infile_num, ref_ctx_count, qry_ctx_count + offset_gid, this_block_size, tmp_ctxobj, this_sketch_index, num_passid_block, sort_idani_block, ani_opt);
		ani_block_print(ref_infile_num, offset_gid, this_block_size, ref_sketch_index, qry_sketch_index, ref_ctx_count, qry_ctx_count, ctx, obj, refname, qryname, refanno, qry_infile_meta, ref_infile_meta, num_passid_block, sort_idani_block, outfp, ani_opt, ani_opt->fmt);

		offset_gid += this_block_size;
	}

	for (int i = 0; i < block_size; i++)
		free(sort_idani_block[i]);
	if (refanno)
		free_read_from_file(refanno, (size_t)ref_infile_num * PATHLEN);
	if (qry_infile_meta)
		free_read_from_file(qry_infile_meta, (size_t)qry_infile_num * sizeof(qry_infile_meta[0]));
	if (ref_infile_meta)
		free_read_from_file(ref_infile_meta, (size_t)ref_infile_num * sizeof(ref_infile_meta[0]));
	free_all(ref_dim_sketch_stat, ref_sketch_index, ref_ctx_count, qry_dim_sketch_stat, qry_sketch_index, qry_ctx_count, tmp_ctxobj, ctx, obj, num_passid_block, sort_idani_block, NULL);
	free_reference_sorted_index(sortedcomb_ctxgid64obj32, ctxgidobj_arr_fsize, sorted_index_is_mmap);
	fclose(fp);
	if (outfp != stdout)
		fclose(outfp);
	// clean xgb model
	// cleanup_model();
	return ctxgidobj_arr_fsize;
}

//
#define DIFF_OBJ_BITS 1
size_t dedup_with_ctxobj_counts(uint32_t *arr, size_t n, co_distance_t **ctxobj_cnt)
{
	*ctxobj_cnt = NULL;
	if (n == 0)
		return 0;

	co_distance_t *tmp_ctxobj_cnt = malloc(n * sizeof(co_distance_t));
	if (!tmp_ctxobj_cnt)
		err(EXIT_FAILURE, "%s(): tmp_ctxobj_cnt malloc failure", __func__);

	size_t j = 0;
	tmp_ctxobj_cnt[0].ctx_ct = 1;
	tmp_ctxobj_cnt[0].diff_obj = arr[0] % 2;
	arr[0] >>= DIFF_OBJ_BITS;

	for (size_t i = 1; i < n; i++)
	{
		if ((arr[i] >> DIFF_OBJ_BITS) == arr[j])
		{
			tmp_ctxobj_cnt[j].ctx_ct++;
			tmp_ctxobj_cnt[j].diff_obj += (arr[i] % 2);
		}
		else
		{
			j++;
			tmp_ctxobj_cnt[j].ctx_ct = 1;
			tmp_ctxobj_cnt[j].diff_obj = arr[i] % 2;
			arr[j] = arr[i] >> DIFF_OBJ_BITS;
		}
	}
	// Trim arrays
	co_distance_t *ctxobj_tmp = realloc(tmp_ctxobj_cnt, (j + 1) * sizeof(co_distance_t));
	*ctxobj_cnt = ctxobj_tmp ? ctxobj_tmp : tmp_ctxobj_cnt;

	return j + 1;
}

ctxgidobj_t *comb_sortedsketch64_2sortedcomb_ctxgid64obj32(unify_sketch_t *ref_result)
{
	// const_comask_init(&ref_result->stats.lco_stat_val);
	uint64_t sketch_size = ref_result->sketch_index[ref_result->infile_num];
	if (sketch_size > (float)UINT32_MAX * LD_FCTR)
		err(EXIT_FAILURE, "%s():sketch_index maximun %lu exceed UINT32_MAX*LF;%f", __func__, sketch_size, (float)UINT32_MAX * LD_FCTR);
	if (ref_result->infile_num >= (1 << GID_NBITS))
		err(EXIT_FAILURE, "%s(): genome numer %d exceed maximum:%u", __func__, ref_result->infile_num, 1 << GID_NBITS);
	if (GID_NBITS + 4 * hclen > 64)
		err(EXIT_FAILURE, "%s(): context_bits_len(%d)+gid_bits_len(%d) exceed 64", __func__, 4 * hclen, GID_NBITS);
	ctxgidobj_t *ctxgidobj_arr = ctxobj64_2ctxgidobj(ref_result->sketch_index, ref_result->comb_sketch, ref_result->infile_num, sketch_size);
	ctxgidobj_sort_array(ctxgidobj_arr, sketch_size);
	return ctxgidobj_arr;
}

sort_sketch_summary_t *summarize_ctxgidobj_arr(ctxgidobj_t *ctxgidobj_arr, uint64_t *sketch_index, uint32_t arrlen, int infile_num)
{

	uint64_t gidmask = UINT64_MAX >> (64 - GID_NBITS);
	num_ctx_cfltobj_t *num_ctx_cfltobj_arr = malloc(infile_num * sizeof(num_ctx_cfltobj_t));
	for (int i = 0; i < infile_num; i++)
	{
		num_ctx_cfltobj_arr[i].num_ctx = sketch_index[i + 1] - sketch_index[i];
		num_ctx_cfltobj_arr[i].num_conflictobj = 0;
	}

	uint32_t num_ctx = 1;
	for (uint32_t i = 1; i < arrlen; i++)
	{
		if (ctxgidobj_arr[i].ctxgid >> GID_NBITS != ctxgidobj_arr[i - 1].ctxgid >> GID_NBITS)
			num_ctx++;
		if (ctxgidobj_arr[i].ctxgid == ctxgidobj_arr[i - 1].ctxgid)
		{
			uint32_t gid = (uint32_t)ctxgidobj_arr[i].ctxgid & gidmask;
			if (num_ctx_cfltobj_arr[gid].num_ctx == (sketch_index[gid + 1] - sketch_index[gid]))
				num_ctx_cfltobj_arr[gid].num_conflictobj++;

			num_ctx_cfltobj_arr[gid].num_ctx--;
		}
	}
	//	for(int i = 0 ; i < infile_num; i++) num_conflictobj +=  num_ctx_cfltobj_arr[i].num_conflictobj;
	sort_sketch_summary_t *sort_sketch_summary = malloc(sizeof(sort_sketch_summary_t));
	sort_sketch_summary->arrlen = arrlen;
	sort_sketch_summary->numgids_perctx = (double)arrlen / num_ctx;
	sort_sketch_summary->num_ctx = num_ctx;
	sort_sketch_summary->infile_num = infile_num;
	sort_sketch_summary->num_ctx_cfltobj_arr = num_ctx_cfltobj_arr;
	return sort_sketch_summary;
}

void free_sort_sketch_summary(sort_sketch_summary_t *sort_sketch_summary)
{
	free(sort_sketch_summary->num_ctx_cfltobj_arr);
	free(sort_sketch_summary);
}

#define CTX(X, Y) (ctx[(size_t)(((X) * ((X) + 1)) / 2 + (Y))])
#define OBJ(X, Y) (obj[(size_t)(((X) * ((X) + 1)) / 2 + (Y))])

/* use sketch variants to caculate distance */
// 1.global sorted comb_sketch64 (i.e. sorted_ctxgidobj_arr or inverted index):
// slow when dist matrix is large, low memory cache efficient, but may be very fast when matrix is small?
void sorted_ctxgidobj_arr2triangle(ctxgidobj_t *ctxgidobj_arr, sort_sketch_summary_t *sort_sketch_summary)
{
	uint64_t gidmask = UINT64_MAX >> (64 - GID_NBITS);
	uint32_t arrlen = sort_sketch_summary->arrlen;
	int infile_num = sort_sketch_summary->infile_num;
	uint16_t *ctx = calloc((size_t)infile_num * (infile_num + 1) / 2, sizeof(uint16_t));
	uint16_t *obj = calloc((size_t)infile_num * (infile_num + 1) / 2, sizeof(uint16_t));

	for (uint32_t i = 0, j; i < arrlen - 1; i = j)
	{
		j = i + 1;
		// find range i..j;
		for (; j < arrlen && (ctxgidobj_arr[i].ctxgid >> GID_NBITS == ctxgidobj_arr[j].ctxgid >> GID_NBITS); j++)
			;

#pragma omp parallel for num_threads(32) schedule(guided)
		for (uint32_t a = i + 1; a < j; a++)
		{
			if (ctxgidobj_arr[a].ctxgid == ctxgidobj_arr[a - 1].ctxgid || ctxgidobj_arr[a].ctxgid == ctxgidobj_arr[a + 1].ctxgid)
				continue; // with confclit object
			uint32_t x = ctxgidobj_arr[a].ctxgid & gidmask;
#pragma omp parallel for num_threads(32) schedule(guided)
			for (uint32_t b = i; b < a; b++)
			{

				if ((b > 0 && ctxgidobj_arr[b].ctxgid == ctxgidobj_arr[b - 1].ctxgid) || ctxgidobj_arr[b].ctxgid == ctxgidobj_arr[b + 1].ctxgid)
					continue;
				uint32_t y = ctxgidobj_arr[b].ctxgid & gidmask;
				CTX(x, y)
				++;
				if (ctxgidobj_arr[a].obj != ctxgidobj_arr[b].obj)
					OBJ(x, y)
				++;
			}
		}
		printf("\ri=%d", i);
	}
	for (int x = 1; x < infile_num; x++)
	{

		for (int y = 0; y < x; y++)
		{
			if (CTX(x, y) > 0)
			{
				printf("%d\t%d\t%d\t%d\t%f\n", x, y, CTX(x, y), OBJ(x, y), (float)OBJ(x, y) / CTX(x, y));
			}
		}
	}
}

// sparse_mem_eff.. seems slower than mem_eff,
int sparse_mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(ani_opt_t *ani_opt)
{
	// initialize
	dim_sketch_stat_t *ref_dim_sketch_stat = read_from_file(test_get_fullpath(ani_opt->refdir, sketch_stat), &file_size);
	const_comask_init(ref_dim_sketch_stat);

	uint64_t gidmask = UINT64_MAX >> (64 - GID_NBITS);
	uint64_t objmask = (1UL << Bitslen.obj) - 1;

	int ref_infile_num = ref_dim_sketch_stat->infile_num;
	// read index
	size_t ctxgidobj_arr_fsize;
	uint64_t *ref_sketch_index = read_from_file(test_get_fullpath(ani_opt->refdir, idx_sketch_suffix), &file_size);
	assert(file_size == (ref_infile_num + 1) * sizeof(ref_sketch_index[0]));
	size_t ref_sketch_size = ref_sketch_index[ref_infile_num];
	ctxgidobj_t *sortedcomb_ctxgid64obj32 = read_from_file(test_get_fullpath(ani_opt->refdir, sorted_comb_ctxgid64obj32), &ctxgidobj_arr_fsize);
	assert(ctxgidobj_arr_fsize == ref_sketch_size * sizeof(sortedcomb_ctxgid64obj32[0]));

	dim_sketch_stat_t *qry_dim_sketch_stat = read_from_file(test_get_fullpath(ani_opt->qrydir, sketch_stat), &file_size);
	int qry_infile_num = qry_dim_sketch_stat->infile_num;
	assert(qry_dim_sketch_stat->hash_id == ref_dim_sketch_stat->hash_id);
	uint64_t *qry_sketch_index = read_from_file(test_get_fullpath(ani_opt->qrydir, idx_sketch_suffix), &file_size);
	size_t qry_sketch_size = qry_sketch_index[qry_infile_num];

	int block_size = BLOCK_SIZE;
	int offset_gid = 0;
	const char *qry_comb_sketch_path = test_get_fullpath(ani_opt->qrydir, combined_sketch_suffix);
	FILE *fp = fopen(qry_comb_sketch_path, "rb");
	if (fp == NULL)
		err(errno, "%s", qry_comb_sketch_path);
	uint64_t *tmp_ctxobj = malloc(qry_sketch_size * sizeof(uint64_t));
	// arry of dynamic gids array per query, note gid <<= DIFF_OBJ_BITS + (ref_obj==qry_obj?0:1);  >>sparse only
	//    Vector *ref_gids_perqry_arr  = malloc( block_size * sizeof(Vector));
	// ref_gids_perqry_arr[i][0] is capacity, ref_gids_perqry_arr[i][1] is length of the inner array
	uint32_t **ref_gids_perqry_arr = malloc(block_size * sizeof(uint32_t *));
	co_distance_t **ctxobj_cnt_perqry_arr = malloc(block_size * sizeof(co_distance_t *));
	size_t *lens = malloc(block_size * sizeof(size_t));
	for (int i = 0; i < block_size; i++)
	{
		ref_gids_perqry_arr[i] = malloc(6 * sizeof(uint32_t)); // 4 for innitialized capacity
		ref_gids_perqry_arr[i][0] = 4;
		ref_gids_perqry_arr[i][1] = 0;
	}

	for (int b = 0; b <= qry_infile_num / block_size; b++)
	{

		int this_block_size = (b == qry_infile_num / block_size) ? (qry_infile_num % block_size) : block_size;
		int this_sketch_size = qry_sketch_index[offset_gid + this_block_size] - qry_sketch_index[offset_gid];
		int read_sketch_size = fread(tmp_ctxobj, sizeof(uint64_t), this_sketch_size, fp);
		uint64_t *this_sketch_index = qry_sketch_index + offset_gid;
		assert(this_sketch_size == read_sketch_size);

#pragma omp parallel for num_threads(ani_opt->p) schedule(guided)
		for (int i = 0; i < this_block_size; i++)
		{
			uint64_t *a = tmp_ctxobj + (this_sketch_index[i] - this_sketch_index[0]);
			size_t a_size = this_sketch_index[i + 1] - this_sketch_index[i];
			size_t *idx = find_first_occurrences_AT_ctxgidobj_arr(a, a_size, sortedcomb_ctxgid64obj32, ref_sketch_size);

			for (int j = 0; j < a_size; j++)
			{
				if (idx[j] == SIZE_MAX)
					continue;
				// skip conlict object;
				if ((j > 0) && ((a[j] >> Bitslen.obj) == (a[j - 1] >> Bitslen.obj)))
					continue;
				if ((j < a_size - 1) && ((a[j] >> Bitslen.obj) == (a[j + 1] >> Bitslen.obj)))
					continue;

				for (int d = idx[j];; d++)
				{
					if ((sortedcomb_ctxgid64obj32[d].ctxgid >> Bitslen.gid) != (a[j] >> Bitslen.obj))
						break;
					uint32_t gid01 = sortedcomb_ctxgid64obj32[d].ctxgid & gidmask << DIFF_OBJ_BITS;
					if ((a[j] & objmask) != sortedcomb_ctxgid64obj32[d].obj)
						gid01 |= 1;
					if (ref_gids_perqry_arr[i][1] + 2 == ref_gids_perqry_arr[i][0])
					{
						ref_gids_perqry_arr[i][0] += 100;
						ref_gids_perqry_arr[i] = realloc(ref_gids_perqry_arr[i], ref_gids_perqry_arr[i][0] * sizeof(ref_gids_perqry_arr[i][0]));
					}
					ref_gids_perqry_arr[i][2 + ref_gids_perqry_arr[i][1]] = gid01;
					ref_gids_perqry_arr[i][1]++;
					// vector_push(&ref_gids_perqry_arr[i],&gid01);
				}
			}
			free(idx);
			qsort(ref_gids_perqry_arr[i] + 2, ref_gids_perqry_arr[i][1], sizeof(uint32_t), qsort_comparator_uint32);
			lens[i] = dedup_with_ctxobj_counts(ref_gids_perqry_arr[i] + 2, ref_gids_perqry_arr[i][1], &ctxobj_cnt_perqry_arr[i]);
		}

		for (int i = 0; i < this_block_size; i++)
		{
			int qry_gid = b * block_size + i;
			int qry_sketch_size = ref_sketch_index[qry_gid + 1] - ref_sketch_index[qry_gid];
			// sparse only code >>
			for (int l = 0; l < lens[i]; l++)
			{
				double dist = ctxobj_cnt_perqry_arr[i][l].diff_obj / ctxobj_cnt_perqry_arr[i][l].ctx_ct;
				double ani = (1 - dist);
#define VEC_GET_AS(type, vec, idx) (((type *)(vec).data)[idx])
				//			if(ctxobj_cnt_perqry_arr[i][l].ctx_ct > 20)
				//          	printf ("%d|%d|%d|%d|%d|%lf|%lf\t",qry_gid,VEC_GET_AS(uint32_t, ref_gids_perqry_arr[i], l),qry_sketch_size,ctxobj_cnt_perqry_arr[i][l].ctx_ct,ctxobj_cnt_perqry_arr[i][l].diff_obj,dist,ani);
			}
			//        printf("\n");
			//<<
		}
		//	for(int i =0 ; i< this_block_size; i++) vector_free(&ref_gids_perqry_arr[i]);//ref_gids_perqry_arr[i].size = 0;
#pragma omp parallel for num_threads(ani_opt->p) schedule(guided)
		for (int i = 0; i < this_block_size; i++)
		{
			ref_gids_perqry_arr[i][1] = 0;
			//			ref_gids_perqry_arr[i][0] = 4;
			free(ctxobj_cnt_perqry_arr[i]);
		}
		offset_gid += this_block_size;
	} // all blocks loop end

	for (int i = 0; i < block_size; i++)
		free(ref_gids_perqry_arr[i]);
	free_all(ref_gids_perqry_arr, ctxobj_cnt_perqry_arr, NULL); ////sparse only code >>
	free_all(ref_dim_sketch_stat, ref_sketch_index, qry_dim_sketch_stat, qry_sketch_index, tmp_ctxobj, NULL);
	free_read_from_file(sortedcomb_ctxgid64obj32, ctxgidobj_arr_fsize);
	fclose(fp);

	return ctxgidobj_arr_fsize;
}

// 2. inverted index(i.e. global sorted_ctxgidobj_arr) X common index(i.e. genome-wise sorted comb_sortedsketch64 ) :
//  ** the fatest method when dist matrix is sparse
void sorted_ctxgidobj_arrXcomb_sortedsketch64(unify_sketch_t *qry_result, ctxgidobj_t *ctxgidobj_arr, sort_sketch_summary_t *sort_sketch_summary)
{

	uint64_t *qry_sketch_index = qry_result->sketch_index;
	uint64_t *qry_comb_sketch = qry_result->comb_sketch;
	int qry_infile_num = qry_result->infile_num;

	uint64_t gidmask = (1UL << Bitslen.gid) - 1;
	uint64_t objmask = (1UL << Bitslen.obj) - 1;
	uint32_t ref_arrlen = sort_sketch_summary->arrlen;
	int ref_infile_num = sort_sketch_summary->infile_num;
	// for self comparision only
	assert(qry_infile_num == ref_infile_num);

	uint16_t *ctx = calloc((size_t)ref_infile_num * (ref_infile_num + 1) / 2, sizeof(uint16_t));
	uint16_t *obj = calloc((size_t)ref_infile_num * (ref_infile_num + 1) / 2, sizeof(uint16_t));
#pragma omp parallel for num_threads(32) schedule(guided)
	for (uint32_t rn = 0; rn < qry_infile_num; rn++)
	{
		//		if(rn > 0) break;
		uint64_t *a = qry_comb_sketch + qry_sketch_index[rn];
		size_t a_size = qry_sketch_index[rn + 1] - qry_sketch_index[rn];
		size_t *idx = find_first_occurrences_AT_ctxgidobj_arr(a, a_size, ctxgidobj_arr, ref_arrlen);
		for (int i = 0; i < a_size; i++)
		{
			if (idx[i] == SIZE_MAX)
				continue;
			// skip conlict object;
			if ((i > 0) && ((a[i] >> Bitslen.obj) == (a[i - 1] >> Bitslen.obj)))
				continue;
			if ((i < a_size - 1) && ((a[i] >> Bitslen.obj) == (a[i + 1] >> Bitslen.obj)))
				continue;

			for (int d = idx[i];; d++)
			{
				uint32_t gid = ctxgidobj_arr[d].ctxgid & gidmask;
				if (gid > rn || (ctxgidobj_arr[d].ctxgid >> Bitslen.gid) != (a[i] >> Bitslen.obj))
					break;
				CTX(rn, gid)
				++;
				if ((a[i] & objmask) != ctxgidobj_arr[d].obj) // wrong: if(a[i] & objmask != ctxgidobj_arr[d].obj) ...
					OBJ(rn, gid)
				++;
			}
		}
		free(idx);
		/*
		#pragma omp critical
				{
					for(int i = 0 ; i <= rn;i++){
						if( CTX(rn,i) > 0 ) printf("\t%d|%d|%d|%d|%f",i,rn, CTX(rn,i),OBJ(rn,i),(float) OBJ(rn,i)/CTX(rn, i));
					}
					printf("\n");
				}
		*/
	}

	for (uint32_t rn = 0; rn < qry_infile_num; rn++)
	{
		printf("%u", rn);
		for (int i = 0; i < rn; i++)
		{
			if (CTX(rn, i) > 0)
				printf("\t%d|%d|%d|%f", i, CTX(rn, i), OBJ(rn, i), (float)OBJ(rn, i) / CTX(rn, i));
		}
		printf("\n");
	}

} // end

// 3. common index X common index (i.e. genome-wise sorted comb_sortedsketch64 ):
//     code hits: use paris wise small sortted arrays overlapping
/*  Advatages:
		** immediate output (no need precompute dist matrix )
		** most convience and memory efficience, no invert indxeing needed
	** when genome are highly similar (dense dist matrix) speed is even faster than inverted index(i.e. global sorted_ctxgidobj_arr) X common index
		** small sortted arrays overlapping is ~ 2X-3X faster than hashtable lookup based overlapping.
*/

typedef struct {
    uint32_t rn;
    double   ani;          /* raw model ANI used internally */
    double   metric;       /* selected distance */
    double   selected_ani;
    double   moe_dist;
    double   naive_dist;
    double   p_dist;
    double   mash_dist;
    double   aaf_dist;
    double   calibrated_ani;
    double   best_ani;
    unsigned char confidence;
    unsigned char best_guarded;
    unsigned char af_pass;
    double   af_qry, blastn_af_qry, af_ref, blastn_af_ref;
    int      XnY_ctx, N_diff_obj, N_diff_obj_section, N_mut2_ctx;
} ani_row_t;

static int cmp_ani_desc(const void *pa, const void *pb) {
    const ani_row_t *a = (const ani_row_t*)pa;
    const ani_row_t *b = (const ani_row_t*)pb;
    if (a->selected_ani < b->selected_ani) return  1;
    if (a->selected_ani > b->selected_ani) return -1;
    if (a->af_qry < b->af_qry) return  1;
    if (a->af_qry > b->af_qry) return -1;
    return (a->rn > b->rn) - (a->rn < b->rn);
}

/* typedef to avoid anonymous-struct warnings from kvec_t params */
typedef kvec_t(ani_row_t) kv_ani_row_t;

/* --- tiny helpers --- */
static inline const char *ani_best_confidence_label(const ani_row_t *r)
{
	return r->best_guarded ? "guarded_low_confidence" : ani_confidence_label(r->confidence);
}

static inline void fill_row_calibration(ani_row_t *r, bool unassembled,
                                        bool enable_best_guard,
                                        const infile_meta_t *qry_asm,
                                        const infile_meta_t *ref_asm)
{
    if (unassembled) {
        r->calibrated_ani = r->ani;
        r->best_ani = r->ani;
        r->confidence = ANI_CONF_HIGH;
        r->best_guarded = 0;
        return;
    }
    r->calibrated_ani = refaf_hgb_predict_ani(r->ani, r->af_ref, (unsigned int)r->XnY_ctx,
                                              (unsigned int)r->N_diff_obj,
                                              (unsigned int)r->N_diff_obj_section,
                                              (unsigned int)r->N_mut2_ctx);
    r->confidence = ani_confidence_from_values(r->ani, r->calibrated_ani, r->af_ref);
    r->best_ani = r->calibrated_ani;
    r->best_guarded = 0;

    const bool query_assembly = has_assembly_meta_record(qry_asm);
    const double min_af = r->af_qry < r->af_ref ? r->af_qry : r->af_ref;
    if (!enable_best_guard ||
        !query_assembly ||
        !infile_meta_complete_like_assembly(qry_asm) ||
        r->ani < 0.958 ||
        r->calibrated_ani < 0.962 ||
        min_af < 0.25 || min_af >= 0.50 ||
        r->XnY_ctx <= 0 || r->af_qry <= 0.0 || r->af_ref <= 0.0)
        return;

    const double qry_ctx = (double)r->XnY_ctx / r->af_qry;
    const double ref_ctx = (double)r->XnY_ctx / r->af_ref;
    double exact = (double)r->XnY_ctx - (double)r->N_diff_obj;
    if (exact < 0.0)
        exact = 0.0;
    const double ctx_exact_mean =
        (mash_ani_from_counts((double)r->XnY_ctx, qry_ctx, ref_ctx) +
         aaf_ani_from_counts((double)r->XnY_ctx, qry_ctx, ref_ctx) +
         mash_ani_from_counts(exact, qry_ctx, ref_ctx) +
         aaf_ani_from_counts(exact, qry_ctx, ref_ctx)) / 4.0;
    if (r->ani - ctx_exact_mean < 0.015)
        return;

    r->best_ani = aaf_ani_from_counts(exact, qry_ctx, ref_ctx);
    r->best_guarded = 1;
}

static inline int selected_metric_request_id(const ani_opt_t *ani_opt)
{
	int s = ani_opt ? abs(ani_opt->s) : 1;
	if (s < 1 || s > ANI_SELECTED_METRIC_COUNT)
		s = 1;
	return s;
}

static inline bool ani_value_available(double value)
{
	return value > 0.0 && value <= 1.0;
}

static inline bool row_best_available(const ani_row_t *row, const ani_opt_t *ani_opt)
{
	return row && ani_opt && !ani_opt->raw_output && !ani_opt->unassembled &&
		   ani_value_available(row->best_ani);
}

static inline bool row_recal_available(const ani_row_t *row, const ani_opt_t *ani_opt)
{
	return row && ani_opt && !ani_opt->raw_output && !ani_opt->unassembled &&
		   ani_value_available(row->calibrated_ani);
}

static inline int selected_metric_effective_id(const ani_row_t *row, const ani_opt_t *ani_opt)
{
	int s = selected_metric_request_id(ani_opt);
	if (ani_opt && ani_opt->unassembled && !ani_opt->unified_metric && s <= 4)
		return 4;
	if (s == 1 && !row_best_available(row, ani_opt))
		s = 2;
	if (s == 2 && !row_recal_available(row, ani_opt))
		s = 3;
	return s;
}

static inline double selected_metric_distance_from_row(const ani_row_t *row, const ani_opt_t *ani_opt)
{
	switch (selected_metric_effective_id(row, ani_opt)) {
	case 1:
		return 1.0 - row->best_ani;
	case 2:
		return 1.0 - row->calibrated_ani;
	case 3:
		return row->moe_dist;
	case 4:
		return row->naive_dist;
	case 5:
		return row->mash_dist;
	case 6:
		return row->aaf_dist;
	case 7:
		return row->af_pass ? row->moe_dist : row->mash_dist;
	case 8:
		return row->af_pass ? row->moe_dist : row->aaf_dist;
	case 9:
		return row->p_dist;
	default:
		return row->metric;
	}
}

static inline const char *selected_metric_name_for_row(const ani_row_t *row, const ani_opt_t *ani_opt)
{
	return select_metrics_header[selected_metric_effective_id(row, ani_opt) - 1];
}

static inline const char *selected_metric_confidence_label(const ani_row_t *row, const ani_opt_t *ani_opt)
{
	if (!ani_opt || ani_opt->raw_output)
		return "raw";
	if (ani_opt->unassembled && !ani_opt->unified_metric &&
		selected_metric_request_id(ani_opt) <= 4)
		return "unassembled";
	switch (selected_metric_effective_id(row, ani_opt)) {
	case 1:
		return ani_best_confidence_label(row);
	case 2:
		return ani_confidence_label(row->confidence);
	default:
		return "raw";
	}
}

static inline void fill_row_base_distances(ani_row_t *row, const ani_features_t *features,
										   const ani_opt_t *ani_opt,
										   uint32_t qry_ctx, uint32_t ref_ctx,
										   double af_q, double af_r)
{
	row->af_pass = ani_report_af_pass(ani_opt, af_q, af_r);
	if (row->af_pass) {
		ani_features_t tmp = *features;
		row->moe_dist = lm3ways_dist_from_features(&tmp);
		tmp = *features;
		row->naive_dist = get_naive_dist(&tmp);
		tmp = *features;
		row->p_dist = get_p_dist(&tmp);
	} else {
		row->moe_dist = ani_opt->e;
		row->naive_dist = ani_opt->e;
		row->p_dist = ani_opt->e;
	}
	row->mash_dist = get_mashD(Bitslen.ctx / 2, ref_ctx, qry_ctx, features->XnY_ctx);
	row->aaf_dist = get_aafD(Bitslen.ctx / 2, ref_ctx, qry_ctx, features->XnY_ctx);
}

static inline void finalize_row_selected_metric(ani_row_t *row, const ani_opt_t *ani_opt)
{
	row->metric = selected_metric_distance_from_row(row, ani_opt);
	row->selected_ani = 1.0 - row->metric;
}

static inline ani_row_t make_selected_output_row(uint32_t rn,
												 const ani_features_t *features,
												 const ani_opt_t *ani_opt,
												 uint32_t qry_ctx,
												 uint32_t ref_ctx,
												 double af_qry,
												 double blastn_af_qry,
												 double af_ref,
												 double blastn_af_ref,
												 const infile_meta_t *qry_asm,
												 const infile_meta_t *ref_asm)
{
	ani_row_t row;
	memset(&row, 0, sizeof(row));
	row.rn = rn;
	row.af_qry = af_qry;
	row.blastn_af_qry = blastn_af_qry;
	row.af_ref = af_ref;
	row.blastn_af_ref = blastn_af_ref;
	row.XnY_ctx = features->XnY_ctx;
	row.N_diff_obj = features->N_diff_obj;
	row.N_diff_obj_section = features->N_diff_obj_section;
	row.N_mut2_ctx = features->N_mut2_ctx;
	fill_row_base_distances(&row, features, ani_opt, qry_ctx, ref_ctx, af_qry, af_ref);
	row.ani = 1.0 - (ani_opt->v ? row.naive_dist : row.moe_dist);
	if (!ani_opt->raw_output && row.af_pass)
		fill_row_calibration(&row, ani_opt->unassembled, ani_best_guard_enabled(ani_opt),
							 qry_asm, ref_asm);
	finalize_row_selected_metric(&row, ani_opt);
	return row;
}

static inline void append_unified_detail_row(kstring_t *ks_out,
                                             const ani_opt_t *ani_opt,
                                             const char *qry_name,
                                             const char *ref_name,
                                             const ani_row_t *r,
                                             const char *annotation)
{
    const double selected_distance = r->metric;
    const double selected_similarity = r->selected_ani;
    const char *metric_name = selected_metric_name_for_row(r, ani_opt);
    const char *confidence = selected_metric_confidence_label(r, ani_opt);
    const char *ref_annotation = annotation_or_na(annotation);
    const double blastn_af_qry = bounded_align_fraction(r->blastn_af_qry);
    const double blastn_af_ref = bounded_align_fraction(r->blastn_af_ref);

    ksprintf(ks_out, "%s\t%s\t%lf\t%lf\t%s\t%s\t%d\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%s\n",
             qry_name, ref_name, selected_similarity, selected_distance, confidence, metric_name,
             r->XnY_ctx, r->af_qry, blastn_af_qry, r->af_ref, blastn_af_ref,
             r->N_diff_obj, r->N_diff_obj_section, r->N_mut2_ctx, ref_annotation);
}

static inline void print_unified_detail_row(FILE *outfp,
                                            const ani_opt_t *ani_opt,
                                            const char *qry_name,
                                            const char *ref_name,
                                            const ani_row_t *r,
                                            const char *annotation)
{
    kstring_t ks = (kstring_t){0, 0, 0};
    append_unified_detail_row(&ks, ani_opt, qry_name, ref_name, r, annotation);
    if (ks.l)
        fwrite(ks.s, 1, ks.l, outfp);
	free(ks.s);
}

static inline double ani_matrix_exception_value(const ani_opt_t *ani_opt)
{
	return ani_opt->s < 0 ? 1.0 - ani_opt->e : ani_opt->e;
}

static inline double ani_matrix_diagonal_value(const ani_opt_t *ani_opt)
{
	return ani_opt->s < 0 ? 1.0 : 0.0;
}

static double ani_matrix_pair_value(const unify_sketch_t *qry_result, uint32_t qn,
									const unify_sketch_t *ref_result, uint32_t rn,
									const ani_opt_t *ani_opt)
{
	uint64_t *arr_qry = qry_result->comb_sketch + qry_result->sketch_index[qn];
	size_t len_qry = qry_result->sketch_index[qn + 1] - qry_result->sketch_index[qn];
	const uint32_t qry_ctx = qry_result->conflict
		? count_ctx_runs_sorted_ctxobj64_local(arr_qry, len_qry, false)
		: (uint32_t)len_qry;
	uint64_t *arr_ref = ref_result->comb_sketch + ref_result->sketch_index[rn];
	size_t len_ref = ref_result->sketch_index[rn + 1] - ref_result->sketch_index[rn];
	const uint32_t ref_ctx = ref_result->conflict
		? count_ctx_runs_sorted_ctxobj64_local(arr_ref, len_ref, ani_opt->ignoreconflict)
		: (uint32_t)len_ref;

	if (qry_ctx == 0 || ref_ctx == 0)
		return ani_matrix_exception_value(ani_opt);

	ani_features_t ani_features;
	get_ani_features_ctx_min_over_conflicts_both_filtered(
		arr_qry, len_qry, arr_ref, len_ref,
		ani_opt->ignoreconflict && ref_result->conflict,
		&ani_features);
	const double af_qry = (double)ani_features.XnY_ctx / (double)qry_ctx;
	const double af_ref = (double)ani_features.XnY_ctx / (double)ref_ctx;
	ani_features_t tmp = ani_features;
	tmp.X_ctx = qry_ctx;
	const double blastn_af_qry = lm3ways_af_ANIb_from_features(&tmp);
	tmp = ani_features;
	tmp.X_ctx = ref_ctx;
	const double blastn_af_ref = lm3ways_af_ANIb_from_features(&tmp);
	ani_row_t outrow = make_selected_output_row(rn, &ani_features, ani_opt,
												qry_ctx, ref_ctx,
												af_qry, blastn_af_qry,
												af_ref, blastn_af_ref,
												infile_meta_at(qry_result, qn),
												infile_meta_at(ref_result, rn));
	return ani_opt->s < 0 ? outrow.selected_ani : outrow.metric;
}

static int ani_choose_dense_matrix_block_size(int ref_n, int qry_n)
{
	if (ref_n <= 0 || qry_n <= 0)
		return 0;
	const size_t bytes_per_query =
		(size_t)ref_n * (sizeof(ctx_mut2_t) + sizeof(obj_section_t));
	if (bytes_per_query == 0)
		return 0;
	size_t budget = GetAvailableMemory() / 4;
	if (budget < bytes_per_query * 64)
		budget = bytes_per_query * 64;
	size_t block = budget / bytes_per_query;
	if (block < 64)
		block = 64;
	if (block > 64)
		block = 64;
	if (block > (size_t)qry_n)
		block = (size_t)qry_n;

	const char *env = getenv("KSSD3A_ANI_MATRIX_BLOCK_SIZE");
	if (!env || env[0] == '\0')
		env = getenv("KSSD3A_MATRIX_BLOCK_SIZE");
	if (env && env[0] != '\0') {
		char *end = NULL;
		errno = 0;
		long requested = strtol(env, &end, 10);
		if (errno == 0 && end != env && *end == '\0' && requested > 0) {
			if (requested > qry_n)
				requested = qry_n;
			block = (size_t)requested;
		}
	}
	return (int)block;
}

static uint32_t *ani_ctx_counts_for_sketch(const unify_sketch_t *sketch, bool ignoreconflict)
{
	uint32_t *counts = calloc((size_t)sketch->infile_num, sizeof(counts[0]));
	if (!counts)
		err(EXIT_FAILURE, "%s(): calloc ctx counts", __func__);
	for (int i = 0; i < sketch->infile_num; ++i) {
		const uint64_t begin = sketch->sketch_index[i];
		const uint64_t end = sketch->sketch_index[i + 1];
		const uint64_t *arr = sketch->comb_sketch + begin;
		const size_t len = (size_t)(end - begin);
		counts[i] = sketch->conflict
						? count_ctx_runs_sorted_ctxobj64_local(arr, len, ignoreconflict)
						: (uint32_t)len;
	}
	return counts;
}

static ctxgidobj_t *load_self_sorted_index_or_build(const unify_sketch_t *sketch,
													const char *sketch_dir,
													size_t *index_bytes,
													bool *from_file,
													bool *is_mmap)
{
	*index_bytes = 0;
	*from_file = false;
	*is_mmap = false;
	const size_t total_entries = (size_t)sketch->sketch_index[sketch->infile_num];
	if (sketch_dir && sketch_dir[0] != '\0' &&
		file_exists_in_folder((char *)sketch_dir, (char *)sorted_comb_ctxgid64obj32)) {
		char *index_path = test_get_fullpath(sketch_dir, sorted_comb_ctxgid64obj32);
		ctxgidobj_t *index = read_reference_sorted_index(index_path, index_bytes, is_mmap);
		free(index_path);
		if (*index_bytes != total_entries * sizeof(index[0]))
			errx(EINVAL, "%s(): sorted index size mismatch", __func__);
		*from_file = true;
		return index;
	}
	if (force_ref_index_requested())
		errx(EXIT_FAILURE, "KSSD3A_FORCE_REF_INDEX requires reference index '%s/%s'",
			 sketch_dir, sorted_comb_ctxgid64obj32);
	*index_bytes = total_entries * sizeof(ctxgidobj_t);
	return comb_sortedsketch64_2sortedcomb_ctxgid64obj32((unify_sketch_t *)sketch);
}

static void free_self_sorted_index(ctxgidobj_t *index, size_t index_bytes,
								   bool from_file, bool is_mmap)
{
	if (!index)
		return;
	if (from_file)
		free_reference_sorted_index(index, index_bytes, is_mmap);
	else
		free(index);
}

static inline size_t ani_lower_offset(int row, int col)
{
	return (size_t)row * (size_t)(row - 1) / 2u + (size_t)col;
}

static int32_t ani_scale_fixed6(double value)
{
	if (!isfinite(value))
		value = 1.0;
	if (value > 2147.0)
		value = 2147.0;
	if (value < -2147.0)
		value = -2147.0;
	return (int32_t)(value >= 0.0
					 ? value * 1000000.0 + 0.5
					 : value * 1000000.0 - 0.5);
}

static void ani_kput_scaled6(kstring_t *s, int32_t scaled)
{
	uint32_t mag;
	if (scaled < 0) {
		kputc('-', s);
		mag = (uint32_t)(-scaled);
	} else {
		mag = (uint32_t)scaled;
	}
	const uint32_t whole = mag / 1000000U;
	uint32_t frac = mag % 1000000U;
	kputuw(whole, s);
	kputc('.', s);
	char frac_buf[6];
	for (int i = 5; i >= 0; --i) {
		frac_buf[i] = (char)('0' + (frac % 10));
		frac /= 10;
	}
	kputsn(frac_buf, 6, s);
}

static double ani_indexed_self_pair_metric(const unify_sketch_t *sketch,
										   const ani_opt_t *ani_opt,
										   int qgid, int rgid,
										   const ani_features_t *features,
										   const uint32_t *qry_ctx_count,
										   const uint32_t *ref_ctx_count)
{
	const uint32_t qry_ctx = qry_ctx_count[qgid];
	const uint32_t ref_ctx = ref_ctx_count[rgid];
	double metric = ani_matrix_exception_value(ani_opt);
	if (qry_ctx > 0 && ref_ctx > 0) {
		const double af_qry = (double)features->XnY_ctx / (double)qry_ctx;
		const double af_ref = (double)features->XnY_ctx / (double)ref_ctx;
		ani_features_t tmp = *features;
		tmp.X_ctx = qry_ctx;
		const double blastn_af_qry = lm3ways_af_ANIb_from_features(&tmp);
		tmp = *features;
		tmp.X_ctx = ref_ctx;
		const double blastn_af_ref = lm3ways_af_ANIb_from_features(&tmp);
		ani_row_t outrow = make_selected_output_row(
			(uint32_t)rgid, features, ani_opt,
			qry_ctx, ref_ctx, af_qry, blastn_af_qry,
			af_ref, blastn_af_ref,
			infile_meta_at(sketch, (uint32_t)qgid),
			infile_meta_at(sketch, (uint32_t)rgid));
		metric = ani_opt->s < 0 ? outrow.selected_ani : outrow.metric;
	}
	return metric;
}

bool comb_sortedsketch64_indexed_self_full(ani_opt_t *ani_opt)
{
	unify_sketch_t *sketch = generic_sketch_parse(ani_opt->qrydir, SKETCH_PARSE_NONE);
	if (sketch->stat_type != 2) {
		free_unify_sketch(sketch);
		return false;
	}

	const int n = sketch->infile_num;
	size_t index_bytes = 0;
	bool index_from_file = false;
	bool index_is_mmap = false;
	ctxgidobj_t *sorted_index =
		load_self_sorted_index_or_build(sketch, ani_opt->qrydir,
										&index_bytes, &index_from_file,
										&index_is_mmap);
	const size_t ref_sksize = (size_t)sketch->sketch_index[n];
	uint32_t *qry_ctx_count = ani_ctx_counts_for_sketch(sketch, false);
	uint32_t *ref_ctx_count = ani_opt->ignoreconflict
								  ? ani_ctx_counts_for_sketch(sketch, true)
								  : qry_ctx_count;

	const int block_size = ani_choose_dense_matrix_block_size(n, n);
	if (block_size <= 0)
		errx(EINVAL, "%s(): invalid ANI full-matrix block size", __func__);
	ctx_mut2_t *ctx = malloc((size_t)n * (size_t)block_size * sizeof(ctx[0]));
	obj_section_t *obj = malloc((size_t)n * (size_t)block_size * sizeof(obj[0]));
	if (!ctx || !obj)
		err(EXIT_FAILURE, "%s(): OOM indexed ANI full-matrix block", __func__);

	const size_t lower_n = (size_t)n * (size_t)(n - 1) / 2u;
	int32_t *lower = lower_n > 0 ? malloc(lower_n * sizeof(lower[0])) : NULL;
	if (lower_n > 0 && !lower)
		err(EXIT_FAILURE, "%s(): OOM indexed ANI full lower matrix", __func__);

	ani_opt_t scan_opt = *ani_opt;
	scan_opt.p = ani_opt->p > 0 ? ani_opt->p : 1;
	for (int offset = 0; offset < n; offset += block_size) {
		const int this_block = (offset + block_size > n) ? n - offset : block_size;
		memset(ctx, 0, (size_t)n * (size_t)this_block * sizeof(ctx[0]));
		memset(obj, 0, (size_t)n * (size_t)this_block * sizeof(obj[0]));
		count_ctx_obj_frm_comb_sketch_section_lower(
			ctx, obj, sorted_index, ref_sksize, n, offset, this_block,
			sketch->comb_sketch + sketch->sketch_index[offset],
			sketch->sketch_index + offset, &scan_opt);

		for (int i = 0; i < this_block; ++i) {
			const int qgid = offset + i;
			for (int rgid = 0; rgid < qgid; ++rgid) {
				ani_features_t features = {
					.XnY_ctx = MCTX(n, i, rgid).num_ctx,
					.N_diff_obj = MOBJ(n, i, rgid).diff_obj,
					.N_diff_obj_section = MOBJ(n, i, rgid).diff_obj_section,
					.N_mut2_ctx = MCTX(n, i, rgid).num_mut2_ctx,
				};
				const double metric = ani_indexed_self_pair_metric(
					sketch, ani_opt, qgid, rgid, &features,
					qry_ctx_count, ref_ctx_count);
				lower[ani_lower_offset(qgid, rgid)] = ani_scale_fixed6(metric);
			}
		}
	}

	FILE *outfp = ani_opt->outf[0] == '\0' ? stdout : fopen(ani_opt->outf, "w");
	if (outfp == NULL)
		err(errno, "%s", ani_opt->outf);
	kstring_t row = {0, 0, NULL};
	for (int rn = 0; rn < n; ++rn) {
		kputc('\t', &row);
		kputs(sketch->gname[rn], &row);
	}
	kputc('\n', &row);
	fwrite(row.s, 1, row.l, outfp);

	const int32_t diag_scaled = ani_scale_fixed6(ani_matrix_diagonal_value(ani_opt));
	for (int qgid = 0; qgid < n; ++qgid) {
		row.l = 0;
		kputs(sketch->gname[qgid], &row);
		for (int rgid = 0; rgid < n; ++rgid) {
			int32_t scaled = diag_scaled;
			if (qgid > rgid)
				scaled = lower[ani_lower_offset(qgid, rgid)];
			else if (qgid < rgid)
				scaled = lower[ani_lower_offset(rgid, qgid)];
			kputc('\t', &row);
			ani_kput_scaled6(&row, scaled);
		}
		kputc('\n', &row);
		fwrite(row.s, 1, row.l, outfp);
	}

	if (outfp != stdout)
		fclose(outfp);
	free(row.s);
	free(lower);
	free(ctx);
	free(obj);
	free(qry_ctx_count);
	if (ref_ctx_count != qry_ctx_count)
		free(ref_ctx_count);
	free_self_sorted_index(sorted_index, index_bytes, index_from_file, index_is_mmap);
	free_unify_sketch(sketch);
	return true;
}

bool comb_sortedsketch64_indexed_self_triangle(ani_opt_t *ani_opt)
{
	unify_sketch_t *sketch = generic_sketch_parse(ani_opt->qrydir, SKETCH_PARSE_NONE);
	if (sketch->stat_type != 2) {
		free_unify_sketch(sketch);
		return false;
	}

	const int n = sketch->infile_num;
	size_t index_bytes = 0;
	bool index_from_file = false;
	bool index_is_mmap = false;
	ctxgidobj_t *sorted_index =
		load_self_sorted_index_or_build(sketch, ani_opt->qrydir,
										&index_bytes, &index_from_file,
										&index_is_mmap);
	const size_t ref_sksize = (size_t)sketch->sketch_index[n];
	uint32_t *qry_ctx_count = ani_ctx_counts_for_sketch(sketch, false);
	uint32_t *ref_ctx_count = ani_opt->ignoreconflict
								  ? ani_ctx_counts_for_sketch(sketch, true)
								  : qry_ctx_count;

	const int block_size = ani_choose_dense_matrix_block_size(n, n);
	if (block_size <= 0)
		errx(EINVAL, "%s(): invalid ANI triangle block size", __func__);
	ctx_mut2_t *ctx = malloc((size_t)n * (size_t)block_size * sizeof(ctx[0]));
	obj_section_t *obj = malloc((size_t)n * (size_t)block_size * sizeof(obj[0]));
	if (!ctx || !obj)
		err(EXIT_FAILURE, "%s(): OOM indexed ANI triangle block", __func__);

	ani_opt_t scan_opt = *ani_opt;
	scan_opt.p = ani_opt->p > 0 ? ani_opt->p : 1;
	kstring_t row = {0, 0, NULL};
	FILE *outfp = ani_opt->outf[0] == '\0' ? stdout : fopen(ani_opt->outf, "w");
	if (outfp == NULL)
		err(errno, "%s", ani_opt->outf);

	for (int offset = 0; offset < n; offset += block_size) {
		const int this_block = (offset + block_size > n) ? n - offset : block_size;
		memset(ctx, 0, (size_t)n * (size_t)this_block * sizeof(ctx[0]));
		memset(obj, 0, (size_t)n * (size_t)this_block * sizeof(obj[0]));
		count_ctx_obj_frm_comb_sketch_section_lower(
			ctx, obj, sorted_index, ref_sksize, n, offset, this_block,
			sketch->comb_sketch + sketch->sketch_index[offset],
			sketch->sketch_index + offset, &scan_opt);

		for (int i = 0; i < this_block; ++i) {
			const int qgid = offset + i;
			row.l = 0;
			kputs(sketch->gname[qgid], &row);
			for (int rgid = 0; rgid < qgid; ++rgid) {
				ani_features_t features = {
					.XnY_ctx = MCTX(n, i, rgid).num_ctx,
					.N_diff_obj = MOBJ(n, i, rgid).diff_obj,
					.N_diff_obj_section = MOBJ(n, i, rgid).diff_obj_section,
					.N_mut2_ctx = MCTX(n, i, rgid).num_mut2_ctx,
				};
				const double metric = ani_indexed_self_pair_metric(
					sketch, ani_opt, qgid, rgid, &features,
					qry_ctx_count, ref_ctx_count);
				ksprintf(&row, "\t%lf", metric);
			}
			if (ani_opt->d)
				ksprintf(&row, "\t%lf", ani_matrix_diagonal_value(ani_opt));
			kputc('\n', &row);
			fwrite(row.s, 1, row.l, outfp);
		}
	}

	if (outfp != stdout)
		fclose(outfp);
	free(row.s);
	free(ctx);
	free(obj);
	free(qry_ctx_count);
	if (ref_ctx_count != qry_ctx_count)
		free(ref_ctx_count);
	free_self_sorted_index(sorted_index, index_bytes, index_from_file, index_is_mmap);
	free_unify_sketch(sketch);
	return true;
}

static inline void kv_append_rows(kv_ani_row_t *dst, const kv_ani_row_t *src)
{
    const size_t n = kv_size(*src);
    if (!n) return;
    const size_t old = kv_size(*dst);
    kv_resize(ani_row_t, *dst, old + n);
    memcpy(&kv_A(*dst, old), &kv_A(*src, 0), n * sizeof(ani_row_t));
    kv_size(*dst) = old + n;
}

static inline void format_rows_to_kstr(const unify_sketch_t *qry, const unify_sketch_t *ref,
                                       uint32_t qn, const kv_ani_row_t *rows,
                                       const ani_opt_t *ani_opt, kstring_t *ks_out)
{
    for (size_t i = 0; i < kv_size(*rows); ++i) {
        const ani_row_t *r = &kv_A(*rows, i);
        append_unified_detail_row(ks_out, ani_opt, qry->gname[qn], ref->gname[r->rn],
                                  r, unify_annotation_at(ref, r->rn));
    }
}


/* Compute (qn,rn); return 1 if it passes filters and fill *row_out */
// called by void comb_sortedsketch64Xcomb_sortedsketch64_filter_and_sort_survivors(ani_opt_t *ani_opt);
static inline int compute_row_if_survivor(
    const unify_sketch_t *qry, const unify_sketch_t *ref,
    uint32_t qn, uint32_t rn, const ani_opt_t *opt,
    ani_row_t *row_out)
{
    uint64_t *arr_qry = qry->comb_sketch + qry->sketch_index[qn];
    size_t    len_qry = qry->sketch_index[qn + 1] - qry->sketch_index[qn];

    uint64_t *arr_ref = ref->comb_sketch + ref->sketch_index[rn];
    size_t    len_ref = ref->sketch_index[rn + 1] - ref->sketch_index[rn];

    ani_features_t f;
    /* Compare one shared context once, choosing the closest object pair when either side keeps conflicts. */
    get_ani_features_ctx_min_over_conflicts_both_filtered(arr_qry, len_qry, arr_ref, len_ref, opt->ignoreconflict && ref->conflict, &f);

    const uint32_t qry_ctx = qry->conflict ? count_ctx_runs_sorted_ctxobj64_local(arr_qry, len_qry, false) : (uint32_t)len_qry;
    const uint32_t ref_ctx = ref->conflict ? count_ctx_runs_sorted_ctxobj64_local(arr_ref, len_ref, opt->ignoreconflict) : (uint32_t)len_ref;
    if (qry_ctx == 0 || ref_ctx == 0) return 0;

    const double af_q = (double)f.XnY_ctx / (double)qry_ctx;
    const double af_r = (double)f.XnY_ctx / (double)ref_ctx;

	if (!ani_report_af_pass(opt, af_q, af_r)) return 0;

    ani_features_t tmp = f;
    tmp.X_ctx = qry_ctx;
    const double blastn_af_q = lm3ways_af_ANIb_from_features(&tmp);
    tmp = f;
    tmp.X_ctx = ref_ctx;
    const double blastn_af_r = lm3ways_af_ANIb_from_features(&tmp);

    ani_row_t r = make_selected_output_row(rn, &f, opt, qry_ctx, ref_ctx,
                                           af_q, blastn_af_q, af_r, blastn_af_r,
                                           infile_meta_at(qry, qn), infile_meta_at(ref, rn));
    if (r.selected_ani <= opt->anicut) return 0;

    *row_out = r;
    return 1;
}

typedef struct {
    uint64_t ctx;
    size_t beg;
    size_t end;
} qry_ctx_run_t;

typedef struct {
    qry_ctx_run_t *runs;
    size_t n_runs;
    size_t *fence;
    int fence_k;
} qry_ctx_lookup_t;

typedef struct {
    uint64_t *keys;  /* ctx + 1; 0 means empty */
    uint64_t *vals;  /* packed beg/end */
    uint32_t mask;
    uint32_t n_runs;
} qry_ctx_hash_lookup_t;

typedef enum {
    QRAW_LOOKUP_SORTED = 0,
    QRAW_LOOKUP_HASH = 1
} qraw_lookup_mode_t;

static inline uint32_t next_pow2_u32_local(uint32_t x)
{
    if (x <= 1)
        return 1;
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return x + 1;
}

static inline uint64_t pack_qry_run_be(size_t beg, size_t end)
{
    return ((uint64_t)(uint32_t)beg << 32) | (uint64_t)(uint32_t)end;
}

static inline size_t unpack_qry_run_beg(uint64_t be) { return (size_t)(uint32_t)(be >> 32); }
static inline size_t unpack_qry_run_end(uint64_t be) { return (size_t)(uint32_t)be; }

static inline uint32_t top_k_bits_ctx_local(uint64_t ctx, int k)
{
    if (k <= 0)
        return 0;
    const int ctx_bits = Bitslen.ctx;
    if (k >= ctx_bits)
        return (uint32_t)ctx;
    return (uint32_t)(ctx >> (ctx_bits - k));
}

static int choose_qry_fence_k(size_t n_runs)
{
    int max_k = Bitslen.ctx < 20 ? Bitslen.ctx : 20;
    if (max_k < 0)
        max_k = 0;
    int k = max_k < 12 ? max_k : 12;
    while (k > 4 && ((size_t)1u << k) > n_runs * 4)
        k--;
    return k;
}

static void build_qry_ctx_lookup(const uint64_t *qry, size_t n, qry_ctx_lookup_t *lookup)
{
    memset(lookup, 0, sizeof(*lookup));
    if (n == 0)
        return;

    lookup->runs = malloc(n * sizeof(*lookup->runs));
    if (!lookup->runs)
        err(EXIT_FAILURE, "%s(): OOM query context runs", __func__);

    const uint8_t nobjbits = Bitslen.obj;
    size_t i = 0;
    while (i < n) {
        const uint64_t ctx = qry[i] >> nobjbits;
        const size_t beg = i;
        do { ++i; } while (i < n && (qry[i] >> nobjbits) == ctx);
        lookup->runs[lookup->n_runs++] = (qry_ctx_run_t){ctx, beg, i};
    }

    lookup->fence_k = choose_qry_fence_k(lookup->n_runs);
    const size_t buckets = (size_t)1u << lookup->fence_k;
    lookup->fence = malloc((buckets + 1) * sizeof(*lookup->fence));
    if (!lookup->fence)
        err(EXIT_FAILURE, "%s(): OOM query context fence", __func__);

    for (size_t t = 0; t <= buckets; ++t)
        lookup->fence[t] = lookup->n_runs;

    size_t next = 0;
    for (size_t r = 0; r < lookup->n_runs; ++r) {
        const uint32_t topk = top_k_bits_ctx_local(lookup->runs[r].ctx, lookup->fence_k);
        while (next <= topk && next <= buckets)
            lookup->fence[next++] = r;
        if (next > buckets)
            break;
    }
    while (next <= buckets)
        lookup->fence[next++] = lookup->n_runs;
}

static void free_qry_ctx_lookup(qry_ctx_lookup_t *lookup)
{
    if (!lookup)
        return;
    free(lookup->runs);
    free(lookup->fence);
    memset(lookup, 0, sizeof(*lookup));
}

static void build_qry_ctx_hash_lookup(const uint64_t *qry, size_t n, qry_ctx_hash_lookup_t *lookup)
{
    memset(lookup, 0, sizeof(*lookup));
    if (n == 0)
        return;
    if (n > UINT32_MAX)
        errx(EXIT_FAILURE, "%s(): query sketch has too many entries for packed hash lookup", __func__);

    const uint32_t cap = next_pow2_u32_local((uint32_t)n * 2u);
    lookup->keys = calloc(cap, sizeof(*lookup->keys));
    lookup->vals = malloc((size_t)cap * sizeof(*lookup->vals));
    if (!lookup->keys || !lookup->vals)
        err(EXIT_FAILURE, "%s(): OOM query context hash", __func__);
    lookup->mask = cap - 1;

    const uint8_t nobjbits = Bitslen.obj;
    size_t i = 0;
    while (i < n) {
        const uint64_t ctx = qry[i] >> nobjbits;
        const size_t beg = i;
        do { ++i; } while (i < n && (qry[i] >> nobjbits) == ctx);
        const uint64_t key = ctx + 1;
        uint32_t pos = (uint32_t)mix64(key) & lookup->mask;
        while (lookup->keys[pos] && lookup->keys[pos] != key)
            pos = (pos + 1) & lookup->mask;
        lookup->keys[pos] = key;
        lookup->vals[pos] = pack_qry_run_be(beg, i);
        lookup->n_runs++;
    }
}

static void free_qry_ctx_hash_lookup(qry_ctx_hash_lookup_t *lookup)
{
    if (!lookup)
        return;
    free(lookup->keys);
    free(lookup->vals);
    memset(lookup, 0, sizeof(*lookup));
}

static inline int lookup_qry_ctx_hash_run(const qry_ctx_hash_lookup_t *lookup, uint64_t ctx,
                                          size_t *beg, size_t *end)
{
    if (!lookup || lookup->n_runs == 0)
        return 0;
    const uint64_t key = ctx + 1;
    uint32_t pos = (uint32_t)mix64(key) & lookup->mask;
    for (;;) {
        const uint64_t cur = lookup->keys[pos];
        if (!cur)
            return 0;
        if (cur == key) {
            const uint64_t be = lookup->vals[pos];
            *beg = unpack_qry_run_beg(be);
            *end = unpack_qry_run_end(be);
            return 1;
        }
        pos = (pos + 1) & lookup->mask;
    }
}

static inline const qry_ctx_run_t *lookup_qry_ctx_run(const qry_ctx_lookup_t *lookup, uint64_t ctx)
{
    if (!lookup || lookup->n_runs == 0)
        return NULL;

    const size_t buckets = (size_t)1u << lookup->fence_k;
    uint32_t topk = top_k_bits_ctx_local(ctx, lookup->fence_k);
    if (topk >= buckets)
        topk = (uint32_t)buckets - 1;

    size_t lo = lookup->fence[topk];
    size_t hi = lookup->fence[topk + 1];
    while (lo < hi) {
        const size_t mid = lo + ((hi - lo) >> 1);
        if (lookup->runs[mid].ctx < ctx)
            lo = mid + 1;
        else
            hi = mid;
    }

    if (lo < lookup->n_runs && lookup->runs[lo].ctx == ctx)
        return &lookup->runs[lo];
    return NULL;
}

static int compute_streamed_ref_row_one_qry(
    const uint64_t *qry, const qry_ctx_lookup_t *qry_lookup,
    const qry_ctx_hash_lookup_t *qry_hash_lookup, qraw_lookup_mode_t lookup_mode,
    const uint64_t *ref, size_t ref_len,
    uint32_t rn, const ani_opt_t *opt, bool ref_conflict,
    const infile_meta_t *qry_asm,
    const infile_meta_t *ref_asm,
    ani_row_t *row_out)
{
    const uint32_t qry_ctx = lookup_mode == QRAW_LOOKUP_HASH
        ? qry_hash_lookup->n_runs
        : (uint32_t)qry_lookup->n_runs;
    if (qry_ctx == 0 || ref_len == 0)
        return 0;
    if (ref_len > UINT32_MAX)
        errx(EXIT_FAILURE, "%s(): reference sketch for genome %u has too many entries", __func__, rn);
    const uint32_t ref_ctx_total = ref_conflict
        ? count_ctx_runs_sorted_ctxobj64_local(ref, ref_len, opt->ignoreconflict)
        : (uint32_t)ref_len;
    if (ref_ctx_total == 0)
        return 0;
    const uint32_t need_X = ani_report_af_needed_ctx(opt, qry_ctx, ref_ctx_total);

    const uint8_t nobjbits = Bitslen.obj;
    const uint64_t objmask = (nobjbits == 64) ? UINT64_MAX : ((1ULL << nobjbits) - 1ULL);
    ani_features_t f;
    memset(&f, 0, sizeof(f));

    uint32_t ref_ctx = 0;
    for (size_t j = 0; j < ref_len; ) {
        const uint64_t ctx = ref[j] >> nobjbits;
        const size_t ref_beg = j;
        do { ++j; } while (j < ref_len && (ref[j] >> nobjbits) == ctx);
        const size_t ref_end = j;

        if (opt->ignoreconflict && ref_end - ref_beg > 1)
            continue;
        ref_ctx++;
        if ((uint32_t)f.XnY_ctx + (ref_ctx_total - ref_ctx) < need_X)
            return 0;

        size_t qry_beg = 0, qry_end = 0;
        if (lookup_mode == QRAW_LOOKUP_HASH) {
            if (!lookup_qry_ctx_hash_run(qry_hash_lookup, ctx, &qry_beg, &qry_end))
                continue;
        } else {
            const qry_ctx_run_t *qry_run = lookup_qry_ctx_run(qry_lookup, ctx);
            if (!qry_run)
                continue;
            qry_beg = qry_run->beg;
            qry_end = qry_run->end;
        }
        f.XnY_ctx++;
        const int min_diff = min_diff_sections_ctxobj64_runs(qry, qry_beg, qry_end,
                                                             ref, ref_beg, ref_end, objmask);
        if (min_diff > 0) {
            f.N_diff_obj++;
            f.N_diff_obj_section += min_diff;
            if (min_diff > 1)
                f.N_mut2_ctx++;
        }
    }

    if (ref_ctx == 0 || (uint32_t)f.XnY_ctx < (uint32_t)opt->ctxcut ||
        (uint32_t)f.XnY_ctx < need_X)
        return 0;

    const double af_q = (double)f.XnY_ctx / (double)qry_ctx;
    const double af_r = (double)f.XnY_ctx / (double)ref_ctx;
	if (!ani_report_af_pass(opt, af_q, af_r))
		return 0;

    ani_features_t tmp = f;
    tmp.X_ctx = qry_ctx;
    const double blastn_af_q = lm3ways_af_ANIb_from_features(&tmp);
    tmp = f;
    tmp.X_ctx = ref_ctx;
    const double blastn_af_r = lm3ways_af_ANIb_from_features(&tmp);

    ani_row_t row = make_selected_output_row(rn, &f, opt, qry_ctx, ref_ctx,
                                             af_q, blastn_af_q, af_r, blastn_af_r,
                                             qry_asm, ref_asm);
    if (row.selected_ani < opt->anicut)
        return 0;
    *row_out = row;
    return 1;
}

int stream_ref_sketches_one_qraw_lookup(ani_opt_t *ani_opt)
{
    unify_sketch_t *qry = generic_sketch_parse(ani_opt->qrydir, ani_query_parse_flags(ani_opt));
    load_infile_meta_for_best_guard(qry, ani_opt->qrydir, ani_opt);
    if (qry->stat_type != 2)
        errx(EXIT_FAILURE, "%s(): small-query streaming requires a 64-bit query sketch", __func__);
    if (qry->infile_num != 1)
        errx(EXIT_FAILURE, "%s(): expected exactly one query sample, found %d",
             __func__, qry->infile_num);
    const_comask_init(&qry->stats.lco_stat_val);

    size_t ref_stat_size = 0;
    char *ref_stat_path = test_get_fullpath(ani_opt->refdir, sketch_stat);
    dim_sketch_stat_t *ref_stat = read_from_file(ref_stat_path, &ref_stat_size);
    free(ref_stat_path);
    if (ref_stat->hash_id != qry->hash_id)
        errx(EXIT_FAILURE, "%s(): hash_id mismatch between query and reference sketches", __func__);
    const uint32_t ref_n = (uint32_t)ref_stat->infile_num;
    char (*refname)[PATHLEN] = (char (*)[PATHLEN])(ref_stat + 1);

    size_t ref_idx_size = 0;
    char *ref_idx_path = test_get_fullpath(ani_opt->refdir, idx_sketch_suffix);
    uint64_t *ref_idx = read_from_file(ref_idx_path, &ref_idx_size);
    free(ref_idx_path);
    if (ref_idx_size != ((size_t)ref_n + 1) * sizeof(ref_idx[0]))
        errx(EXIT_FAILURE, "%s(): %s/%s has %zu bytes, expected %zu",
             __func__, ani_opt->refdir, idx_sketch_suffix, ref_idx_size,
             ((size_t)ref_n + 1) * sizeof(ref_idx[0]));

    char (*refanno)[PATHLEN] = read_optional_sketch_annotations(ani_opt->refdir, (int)ref_n);
    infile_meta_t *ref_infile_meta =
        ani_best_guard_enabled(ani_opt) ? read_optional_sketch_infile_meta_stats(ani_opt->refdir, (int)ref_n) : NULL;

    const uint64_t *qry_arr = qry->comb_sketch + qry->sketch_index[0];
    const size_t qry_len = (size_t)(qry->sketch_index[1] - qry->sketch_index[0]);
    qraw_lookup_mode_t lookup_mode = QRAW_LOOKUP_HASH;
    const char *lookup_env = getenv("KSSD3A_QRAW_LOOKUP");
    if (lookup_env && lookup_env[0]) {
        if (strcmp(lookup_env, "sorted") == 0 || strcmp(lookup_env, "fence") == 0 ||
            strcmp(lookup_env, "fencepost") == 0) {
            lookup_mode = QRAW_LOOKUP_SORTED;
        } else if (strcmp(lookup_env, "hash") == 0) {
            lookup_mode = QRAW_LOOKUP_HASH;
        } else {
            errx(EXIT_FAILURE, "%s(): invalid KSSD3A_QRAW_LOOKUP='%s' (use hash or sorted)",
                 __func__, lookup_env);
        }
    }
    qry_ctx_lookup_t qry_lookup;
    qry_ctx_hash_lookup_t qry_hash_lookup;
    memset(&qry_lookup, 0, sizeof(qry_lookup));
    memset(&qry_hash_lookup, 0, sizeof(qry_hash_lookup));
    if (lookup_mode == QRAW_LOOKUP_HASH)
        build_qry_ctx_hash_lookup(qry_arr, qry_len, &qry_hash_lookup);
    else
        build_qry_ctx_lookup(qry_arr, qry_len, &qry_lookup);

    FILE *outfp = ani_opt->outf[0] == '\0' ? stdout : fopen(ani_opt->outf, "w");
    if (!outfp)
        err(errno, "%s", ani_opt->outf);
    print_ani_detail_header(outfp, ani_opt, true);

    char *ref_comb_path = test_get_fullpath(ani_opt->refdir, combined_sketch_suffix);
    FILE *ref_fp = fopen(ref_comb_path, "rb");
    if (!ref_fp)
        err(errno, "%s", ref_comb_path);
    free(ref_comb_path);
    if (ref_idx[0] != 0) {
        const off_t offset = (off_t)(ref_idx[0] * sizeof(uint64_t));
        if (fseeko(ref_fp, offset, SEEK_SET) != 0)
            err(errno, "%s(): failed to seek reference sketch", __func__);
    }

    const int P = ani_opt->p > 0 ? ani_opt->p : 1;
    const size_t max_block_entries = (size_t)1u << 22; /* 32 MiB of uint64_t sketch data */
    uint64_t *ref_block = NULL;
    size_t ref_block_cap = 0;
    kv_ani_row_t survivors;
    kv_init(survivors);

    for (uint32_t block_start = 0; block_start < ref_n; ) {
        const uint64_t block_begin = ref_idx[block_start];
        uint32_t block_end_gid = block_start;
        uint64_t block_end = block_begin;
        while (block_end_gid < ref_n) {
            const uint64_t next_end = ref_idx[block_end_gid + 1];
            if (block_end_gid > block_start && next_end - block_begin > max_block_entries)
                break;
            block_end_gid++;
            block_end = next_end;
        }
        const size_t block_entries = (size_t)(block_end - block_begin);
        if (block_entries > ref_block_cap) {
            uint64_t *new_block = realloc(ref_block, block_entries * sizeof(*ref_block));
            if (!new_block)
                err(EXIT_FAILURE, "%s(): OOM reference stream block", __func__);
            ref_block = new_block;
            ref_block_cap = block_entries;
        }
        if (block_entries > 0) {
            const size_t got = fread(ref_block, sizeof(*ref_block), block_entries, ref_fp);
            if (got != block_entries)
                errx(EXIT_FAILURE, "%s(): short read from reference sketch: got %zu, expected %zu",
                     __func__, got, block_entries);
        }

        kv_ani_row_t *tls = calloc((size_t)P, sizeof(*tls));
        if (!tls)
            err(EXIT_FAILURE, "%s(): OOM survivor buffers", __func__);
        for (int t = 0; t < P; ++t)
            kv_init(tls[t]);

        const uint32_t block_count = block_end_gid - block_start;
#pragma omp parallel num_threads(P)
        {
            const int tid = omp_get_thread_num();
            ani_row_t row;
#pragma omp for schedule(dynamic, 64)
            for (uint32_t bi = 0; bi < block_count; ++bi) {
                const uint32_t rn = block_start + bi;
                const size_t ref_off = (size_t)(ref_idx[rn] - block_begin);
                const size_t ref_len = (size_t)(ref_idx[rn + 1] - ref_idx[rn]);
                const uint64_t *ref_arr = ref_len ? ref_block + ref_off : NULL;
                if (compute_streamed_ref_row_one_qry(qry_arr, &qry_lookup, &qry_hash_lookup,
                                                     lookup_mode, ref_arr, ref_len,
                                                     rn, ani_opt, ref_stat->conflict,
                                                     infile_meta_at(qry, 0),
                                                     ref_infile_meta ? &ref_infile_meta[rn] : NULL,
                                                     &row))
                    kv_push(ani_row_t, tls[tid], row);
            }
        }

        for (int t = 0; t < P; ++t) {
            if (kv_size(tls[t]))
                kv_append_rows(&survivors, &tls[t]);
            kv_destroy(tls[t]);
        }
        free(tls);
        block_start = block_end_gid;
    }

    if (kv_size(survivors))
        qsort(&kv_A(survivors, 0), kv_size(survivors), sizeof(ani_row_t), cmp_ani_desc);

    size_t out_n = kv_size(survivors);
    if (ani_opt->ntop > 0 && (size_t)ani_opt->ntop < out_n)
        out_n = (size_t)ani_opt->ntop;
    for (size_t i = 0; i < out_n; ++i) {
        const ani_row_t *r = &kv_A(survivors, i);
        print_unified_detail_row(outfp, ani_opt, qry->gname[0], refname[r->rn],
                                 r, annotation_at(refanno, r->rn));
    }

    kv_destroy(survivors);
    free(ref_block);
    fclose(ref_fp);
    if (outfp != stdout)
        fclose(outfp);
    free_qry_ctx_lookup(&qry_lookup);
    free_qry_ctx_hash_lookup(&qry_hash_lookup);
    if (refanno)
        free_read_from_file(refanno, (size_t)ref_n * PATHLEN);
    if (ref_infile_meta)
        free_read_from_file(ref_infile_meta, (size_t)ref_n * sizeof(ref_infile_meta[0]));
    free_read_from_file(ref_idx, ref_idx_size);
    free_read_from_file(ref_stat, ref_stat_size);
    free_unify_sketch(qry);
    return 0;
}

typedef struct {
    uint64_t ctx;
    uint32_t qid;
    uint32_t beg;
    uint32_t end;
} qraw_multi_qctx_run_t;

typedef struct {
    qraw_multi_qctx_run_t *runs;
    size_t n_runs;
    size_t *fence;
    int fence_k;
    uint32_t *qry_ctx_counts;
    uint32_t qry_n;
} qraw_multi_lookup_t;

static int cmp_qraw_multi_qctx_run(const void *pa, const void *pb)
{
    const qraw_multi_qctx_run_t *a = (const qraw_multi_qctx_run_t *)pa;
    const qraw_multi_qctx_run_t *b = (const qraw_multi_qctx_run_t *)pb;
    if (a->ctx < b->ctx) return -1;
    if (a->ctx > b->ctx) return 1;
    return (a->qid > b->qid) - (a->qid < b->qid);
}

static void build_qraw_multi_lookup(const unify_sketch_t *qry, qraw_multi_lookup_t *lookup)
{
    memset(lookup, 0, sizeof(*lookup));
    const uint32_t Q = (uint32_t)qry->infile_num;
    lookup->qry_n = Q;
    lookup->qry_ctx_counts = calloc(Q, sizeof(*lookup->qry_ctx_counts));
    if (!lookup->qry_ctx_counts)
        err(EXIT_FAILURE, "%s(): OOM query context counts", __func__);

    const uint64_t total_entries = qry->sketch_index[Q];
    if (total_entries == 0)
        return;
    if (total_entries > UINT32_MAX)
        errx(EXIT_FAILURE, "%s(): query sketch has too many entries for packed multi-query lookup", __func__);

    lookup->runs = malloc((size_t)total_entries * sizeof(*lookup->runs));
    if (!lookup->runs)
        err(EXIT_FAILURE, "%s(): OOM multi-query context runs", __func__);

    const uint8_t nobjbits = Bitslen.obj;
    for (uint32_t q = 0; q < Q; ++q) {
        size_t i = (size_t)qry->sketch_index[q];
        const size_t end_q = (size_t)qry->sketch_index[q + 1];
        while (i < end_q) {
            const uint64_t ctx = qry->comb_sketch[i] >> nobjbits;
            const size_t beg = i;
            do { ++i; } while (i < end_q && (qry->comb_sketch[i] >> nobjbits) == ctx);
            lookup->runs[lookup->n_runs++] = (qraw_multi_qctx_run_t){
                .ctx = ctx, .qid = q, .beg = (uint32_t)beg, .end = (uint32_t)i
            };
            lookup->qry_ctx_counts[q]++;
        }
    }

    if (lookup->n_runs > 1)
        qsort(lookup->runs, lookup->n_runs, sizeof(*lookup->runs), cmp_qraw_multi_qctx_run);

    lookup->fence_k = choose_qry_fence_k(lookup->n_runs);
    const size_t buckets = (size_t)1u << lookup->fence_k;
    lookup->fence = malloc((buckets + 1) * sizeof(*lookup->fence));
    if (!lookup->fence)
        err(EXIT_FAILURE, "%s(): OOM multi-query context fence", __func__);

    for (size_t t = 0; t <= buckets; ++t)
        lookup->fence[t] = lookup->n_runs;

    size_t next = 0;
    for (size_t r = 0; r < lookup->n_runs; ++r) {
        const uint32_t topk = top_k_bits_ctx_local(lookup->runs[r].ctx, lookup->fence_k);
        while (next <= topk && next <= buckets)
            lookup->fence[next++] = r;
        if (next > buckets)
            break;
    }
    while (next <= buckets)
        lookup->fence[next++] = lookup->n_runs;
}

static void free_qraw_multi_lookup(qraw_multi_lookup_t *lookup)
{
    if (!lookup)
        return;
    free(lookup->runs);
    free(lookup->fence);
    free(lookup->qry_ctx_counts);
    memset(lookup, 0, sizeof(*lookup));
}

static inline void qraw_multi_lookup_ctx_range(const qraw_multi_lookup_t *lookup, uint64_t ctx,
                                               size_t *lo_out, size_t *hi_out)
{
    *lo_out = *hi_out = 0;
    if (!lookup || lookup->n_runs == 0)
        return;

    const size_t buckets = (size_t)1u << lookup->fence_k;
    uint32_t topk = top_k_bits_ctx_local(ctx, lookup->fence_k);
    if (topk >= buckets)
        topk = (uint32_t)buckets - 1;

    size_t lo = lookup->fence[topk];
    size_t hi = lookup->fence[topk + 1];
    while (lo < hi) {
        const size_t mid = lo + ((hi - lo) >> 1);
        if (lookup->runs[mid].ctx < ctx)
            lo = mid + 1;
        else
            hi = mid;
    }
    if (lo >= lookup->n_runs || lookup->runs[lo].ctx != ctx)
        return;

    size_t end = lo + 1;
    while (end < lookup->n_runs && lookup->runs[end].ctx == ctx)
        end++;
    *lo_out = lo;
    *hi_out = end;
}

int stream_ref_sketches_multi_qraw_sortedindex(ani_opt_t *ani_opt)
{
    unify_sketch_t *qry = generic_sketch_parse(ani_opt->qrydir, ani_query_parse_flags(ani_opt));
    load_infile_meta_for_best_guard(qry, ani_opt->qrydir, ani_opt);
    if (qry->stat_type != 2)
        errx(EXIT_FAILURE, "%s(): small-query streaming requires 64-bit query sketches", __func__);
    if (qry->infile_num <= 1)
        errx(EXIT_FAILURE, "%s(): expected multiple query samples, found %d",
             __func__, qry->infile_num);
    const_comask_init(&qry->stats.lco_stat_val);
    const uint32_t Q = (uint32_t)qry->infile_num;

    size_t ref_stat_size = 0;
    char *ref_stat_path = test_get_fullpath(ani_opt->refdir, sketch_stat);
    dim_sketch_stat_t *ref_stat = read_from_file(ref_stat_path, &ref_stat_size);
    free(ref_stat_path);
    if (ref_stat->hash_id != qry->hash_id)
        errx(EXIT_FAILURE, "%s(): hash_id mismatch between query and reference sketches", __func__);
    const uint32_t ref_n = (uint32_t)ref_stat->infile_num;
    char (*refname)[PATHLEN] = (char (*)[PATHLEN])(ref_stat + 1);

    size_t ref_idx_size = 0;
    char *ref_idx_path = test_get_fullpath(ani_opt->refdir, idx_sketch_suffix);
    uint64_t *ref_idx = read_from_file(ref_idx_path, &ref_idx_size);
    free(ref_idx_path);
    if (ref_idx_size != ((size_t)ref_n + 1) * sizeof(ref_idx[0]))
        errx(EXIT_FAILURE, "%s(): %s/%s has %zu bytes, expected %zu",
             __func__, ani_opt->refdir, idx_sketch_suffix, ref_idx_size,
             ((size_t)ref_n + 1) * sizeof(ref_idx[0]));

    char (*refanno)[PATHLEN] = read_optional_sketch_annotations(ani_opt->refdir, (int)ref_n);
    infile_meta_t *ref_infile_meta =
        ani_best_guard_enabled(ani_opt) ? read_optional_sketch_infile_meta_stats(ani_opt->refdir, (int)ref_n) : NULL;
    qraw_multi_lookup_t lookup;
    build_qraw_multi_lookup(qry, &lookup);

    FILE *outfp = ani_opt->outf[0] == '\0' ? stdout : fopen(ani_opt->outf, "w");
    if (!outfp)
        err(errno, "%s", ani_opt->outf);
    print_ani_detail_header(outfp, ani_opt, true);

    char *ref_comb_path = test_get_fullpath(ani_opt->refdir, combined_sketch_suffix);
    FILE *ref_fp = fopen(ref_comb_path, "rb");
    if (!ref_fp)
        err(errno, "%s", ref_comb_path);
    free(ref_comb_path);
    if (ref_idx[0] != 0) {
        const off_t offset = (off_t)(ref_idx[0] * sizeof(uint64_t));
        if (fseeko(ref_fp, offset, SEEK_SET) != 0)
            err(errno, "%s(): failed to seek reference sketch", __func__);
    }

    const int P = ani_opt->p > 0 ? ani_opt->p : 1;
    const size_t max_block_entries = (size_t)1u << 22;
    uint64_t *ref_block = NULL;
    size_t ref_block_cap = 0;
    kv_ani_row_t *survivors = calloc(Q, sizeof(*survivors));
    if (!survivors)
        err(EXIT_FAILURE, "%s(): OOM survivor vectors", __func__);
    for (uint32_t q = 0; q < Q; ++q)
        kv_init(survivors[q]);

    const uint8_t nobjbits = Bitslen.obj;
    const uint64_t objmask = (nobjbits == 64) ? UINT64_MAX : ((1ULL << nobjbits) - 1ULL);

    for (uint32_t block_start = 0; block_start < ref_n; ) {
        const uint64_t block_begin = ref_idx[block_start];
        uint32_t block_end_gid = block_start;
        uint64_t block_end = block_begin;
        while (block_end_gid < ref_n) {
            const uint64_t next_end = ref_idx[block_end_gid + 1];
            if (block_end_gid > block_start && next_end - block_begin > max_block_entries)
                break;
            block_end_gid++;
            block_end = next_end;
        }
        const size_t block_entries = (size_t)(block_end - block_begin);
        if (block_entries > ref_block_cap) {
            uint64_t *new_block = realloc(ref_block, block_entries * sizeof(*ref_block));
            if (!new_block)
                err(EXIT_FAILURE, "%s(): OOM reference stream block", __func__);
            ref_block = new_block;
            ref_block_cap = block_entries;
        }
        if (block_entries > 0) {
            const size_t got = fread(ref_block, sizeof(*ref_block), block_entries, ref_fp);
            if (got != block_entries)
                errx(EXIT_FAILURE, "%s(): short read from reference sketch: got %zu, expected %zu",
                     __func__, got, block_entries);
        }

        kv_ani_row_t *tls = calloc((size_t)P * Q, sizeof(*tls));
        if (!tls)
            err(EXIT_FAILURE, "%s(): OOM thread-local survivor vectors", __func__);
        for (size_t i = 0; i < (size_t)P * Q; ++i)
            kv_init(tls[i]);

        const uint32_t block_count = block_end_gid - block_start;
#pragma omp parallel num_threads(P)
        {
            const int tid = omp_get_thread_num();
            ani_features_t *features = calloc(Q, sizeof(*features));
            if (!features)
                err(EXIT_FAILURE, "%s(): OOM per-reference features", __func__);

#pragma omp for schedule(dynamic, 64)
            for (uint32_t bi = 0; bi < block_count; ++bi) {
                const uint32_t rn = block_start + bi;
                const size_t ref_off = (size_t)(ref_idx[rn] - block_begin);
                const size_t ref_len = (size_t)(ref_idx[rn + 1] - ref_idx[rn]);
                const uint64_t *ref = ref_len ? ref_block + ref_off : NULL;
                memset(features, 0, (size_t)Q * sizeof(*features));
                uint32_t ref_ctx = 0;

                for (size_t j = 0; j < ref_len; ) {
                    const uint64_t ctx = ref[j] >> nobjbits;
                    const size_t ref_beg = j;
                    do { ++j; } while (j < ref_len && (ref[j] >> nobjbits) == ctx);
                    const size_t ref_end = j;

                    if (ani_opt->ignoreconflict && ref_end - ref_beg > 1)
                        continue;
                    ref_ctx++;

                    size_t qlo = 0, qhi = 0;
                    qraw_multi_lookup_ctx_range(&lookup, ctx, &qlo, &qhi);
                    for (size_t qi = qlo; qi < qhi; ++qi) {
                        const qraw_multi_qctx_run_t *qr = &lookup.runs[qi];
                        ani_features_t *f = &features[qr->qid];
                        f->XnY_ctx++;
                        const int min_diff = min_diff_sections_ctxobj64_runs(
                            qry->comb_sketch, qr->beg, qr->end,
                            ref, ref_beg, ref_end, objmask);
                        if (min_diff > 0) {
                            f->N_diff_obj++;
                            f->N_diff_obj_section += min_diff;
                            if (min_diff > 1)
                                f->N_mut2_ctx++;
                        }
                    }
                }

                if (ref_ctx == 0)
                    continue;

                for (uint32_t q = 0; q < Q; ++q) {
                    ani_features_t f = features[q];
                    const uint32_t qry_ctx = lookup.qry_ctx_counts[q];
                    if (qry_ctx == 0 || (uint32_t)f.XnY_ctx < (uint32_t)ani_opt->ctxcut)
                        continue;

                    const double af_q = (double)f.XnY_ctx / (double)qry_ctx;
                    const double af_r = (double)f.XnY_ctx / (double)ref_ctx;
					if (!ani_report_af_pass(ani_opt, af_q, af_r))
						continue;

                    ani_features_t tmp = f;
                    tmp.X_ctx = qry_ctx;
                    const double blastn_af_q = lm3ways_af_ANIb_from_features(&tmp);
                    tmp = f;
                    tmp.X_ctx = ref_ctx;
                    const double blastn_af_r = lm3ways_af_ANIb_from_features(&tmp);

                    ani_row_t row = make_selected_output_row(rn, &f, ani_opt, qry_ctx, ref_ctx,
                                                             af_q, blastn_af_q, af_r, blastn_af_r,
                                                             infile_meta_at(qry, q),
                                                             ref_infile_meta ? &ref_infile_meta[rn] : NULL);
                    if (row.selected_ani < ani_opt->anicut)
                        continue;
                    kv_push(ani_row_t, tls[(size_t)tid * Q + q], row);
                }
            }

            free(features);
        }

        for (int t = 0; t < P; ++t) {
            for (uint32_t q = 0; q < Q; ++q) {
                kv_ani_row_t *src = &tls[(size_t)t * Q + q];
                if (kv_size(*src))
                    kv_append_rows(&survivors[q], src);
                kv_destroy(*src);
            }
        }
        free(tls);
        block_start = block_end_gid;
    }

    for (uint32_t q = 0; q < Q; ++q) {
        if (kv_size(survivors[q]))
            qsort(&kv_A(survivors[q], 0), kv_size(survivors[q]), sizeof(ani_row_t), cmp_ani_desc);
        size_t out_n = kv_size(survivors[q]);
        if (ani_opt->ntop > 0 && (size_t)ani_opt->ntop < out_n)
            out_n = (size_t)ani_opt->ntop;
        for (size_t i = 0; i < out_n; ++i) {
            const ani_row_t *r = &kv_A(survivors[q], i);
            print_unified_detail_row(outfp, ani_opt, qry->gname[q], refname[r->rn],
                                     r, annotation_at(refanno, r->rn));
        }
        kv_destroy(survivors[q]);
    }

    free(survivors);
    free(ref_block);
    fclose(ref_fp);
    if (outfp != stdout)
        fclose(outfp);
    free_qraw_multi_lookup(&lookup);
    if (refanno)
        free_read_from_file(refanno, (size_t)ref_n * PATHLEN);
    if (ref_infile_meta)
        free_read_from_file(ref_infile_meta, (size_t)ref_n * sizeof(ref_infile_meta[0]));
    free_read_from_file(ref_idx, ref_idx_size);
    free_read_from_file(ref_stat, ref_stat_size);
    free_unify_sketch(qry);
    return 0;
}

/* --- main --- */
/** 
void comb_sortedsketch64Xcomb_sortedsketch64_filter_and_sort_survivors(ani_opt_t *ani_opt)
{
    unify_sketch_t *qry = generic_sketch_parse(ani_opt->qrydir, ani_query_parse_flags(ani_opt));
    unify_sketch_t *ref = generic_sketch_parse(ani_opt->refdir, ani_ref_parse_flags(ani_opt));
    load_infile_meta_for_best_guard(qry, ani_opt->qrydir, ani_opt);
    load_infile_meta_for_best_guard(ref, ani_opt->refdir, ani_opt);
    if (ref->conflict)
        errx(EXIT_FAILURE, "%s(): ref '%s' contains conflicting objects!", __func__, ani_opt->refdir);

    const uint32_t Q = qry->infile_num;
    const uint32_t R = ref->infile_num;

    FILE *outfp = (ani_opt->outf[0] == '\0') ? stdout : fopen(ani_opt->outf, "w");
    if (!outfp) err(errno, "%s", ani_opt->outf);
    print_ani_detail_header(outfp, ani_opt, false);

    const int P = (ani_opt->p > 0) ? ani_opt->p : 1;

    // ---- P==1: fast serial path (logic unchanged) ---- //
    if (P == 1) {
        for (uint32_t qn = 0; qn < Q; ++qn) {
            kv_ani_row_t surv; kv_init(surv);
            kstring_t ks_out = (kstring_t){0,0,0};

            ani_row_t tmp;
            for (uint32_t rn = 0; rn < R; ++rn) {
                if (compute_row_if_survivor(qry, ref, qn, rn, ani_opt, &tmp))
                    kv_push(ani_row_t, surv, tmp);
            }

            if (kv_size(surv))
                qsort(&kv_A(surv,0), kv_size(surv), sizeof(ani_row_t), cmp_ani_desc);

            format_rows_to_kstr(qry, ref, qn, &surv, ani_opt, &ks_out);
            if (ks_out.l) fwrite(ks_out.s, 1, ks_out.l, outfp);

            kv_destroy(surv);
            free(ks_out.s);
        }
        if (outfp != stdout) fclose(outfp);
        return;
    }

#ifdef PRINT_UNORDERED
    // ---- Parallel, unordered print (as-ready) ---- //
    #pragma omp parallel num_threads(P)
    {
        #pragma omp single nowait
        {
            const int many_queries = (Q >= (uint32_t)P);

            for (uint32_t qn = 0; qn < Q; ++qn) {
                #pragma omp task firstprivate(qn, many_queries) shared(qry,ref,ani_opt,outfp)
                {
                    kv_ani_row_t surv; kv_init(surv);

                    if (many_queries) {
                        ani_row_t tmp;
                        for (uint32_t rn = 0; rn < R; ++rn) {
                            if (compute_row_if_survivor(qry, ref, qn, rn, ani_opt, &tmp))
                                kv_push(ani_row_t, surv, tmp);
                        }
                    } else {
                        #pragma omp taskgroup
                        {
                            const uint32_t CHUNK = 256;
                            for (uint32_t start = 0; start < R; start += CHUNK) {
                                const uint32_t end = (start + CHUNK < R) ? start + CHUNK : R;
                                #pragma omp task firstprivate(start,end,qn) shared(qry,ref,ani_opt,surv)
                                {
                                    kv_ani_row_t local; kv_init(local);
                                    ani_row_t tmp;
                                    for (uint32_t rn = start; rn < end; ++rn) {
                                        if (compute_row_if_survivor(qry, ref, qn, rn, ani_opt, &tmp))
                                            kv_push(ani_row_t, local, tmp);
                                    }
                                    if (kv_size(local)) {
                                        #pragma omp critical (append_survivors)
                                        kv_append_rows(&surv, &local);
                                    }
                                    kv_destroy(local);
                                }
                            }
                        }
                    }

                    if (kv_size(surv))
                        qsort(&kv_A(surv,0), kv_size(surv), sizeof(ani_row_t), cmp_ani_desc);

                    kstring_t ks_out = (kstring_t){0,0,0};
                    format_rows_to_kstr(qry, ref, qn, &surv, ani_opt, &ks_out);

                    if (ks_out.l) {
                        #pragma omp critical (print)
                        fwrite(ks_out.s, 1, ks_out.l, outfp);
                    }

                    free(ks_out.s);
                    kv_destroy(surv);
                } // compute+print task //
            }
        } // single //
    } // parallel //
#else
    // ---- Parallel, ordered print: compute in parallel; only print is serialized ---- //

    // tokens: comp[qn] signals compute done; order[i] forms a print chain i->i+1 //
    char *comp  = (char*)calloc(Q,   1);
    char *order = (char*)calloc(Q+1, 1);
    if (!comp || !order) err(EXIT_FAILURE, "calloc tokens");

    // per-qn output buffers produced by compute, consumed by print //
    kstring_t *ks_arr = (kstring_t*)calloc(Q, sizeof(kstring_t));
    if (!ks_arr) err(EXIT_FAILURE, "calloc ks_arr");

    #pragma omp parallel num_threads(P)
    {
        #pragma omp single nowait
        {
            const int many_queries = (Q >= (uint32_t)P);

            for (uint32_t qn = 0; qn < Q; ++qn) {

                // 1) compute task: fully parallel across qn //
                #pragma omp task firstprivate(qn, many_queries) shared(qry,ref,ani_opt,ks_arr,comp) depend(out: comp[qn])
                {
                    kv_ani_row_t surv; kv_init(surv);

                    if (many_queries) {
                        ani_row_t tmp;
                        for (uint32_t rn = 0; rn < R; ++rn) {
                            if (compute_row_if_survivor(qry, ref, qn, rn, ani_opt, &tmp))
                                kv_push(ani_row_t, surv, tmp);
                        }
                    } else {
                        #pragma omp taskgroup
                        {
                            const uint32_t CHUNK = 256;
                            for (uint32_t start = 0; start < R; start += CHUNK) {
                                const uint32_t end = (start + CHUNK < R) ? start + CHUNK : R;
                                #pragma omp task firstprivate(start,end,qn) shared(qry,ref,ani_opt,surv)
                                {
                                    kv_ani_row_t local; kv_init(local);
                                    ani_row_t tmp;
                                    for (uint32_t rn = start; rn < end; ++rn) {
                                        if (compute_row_if_survivor(qry, ref, qn, rn, ani_opt, &tmp))
                                            kv_push(ani_row_t, local, tmp);
                                    }
                                    if (kv_size(local)) {
                                        #pragma omp critical (append_survivors)
                                        kv_append_rows(&surv, &local);
                                    }
                                    kv_destroy(local);
                                }
                            }
                        }
                    }

                    if (kv_size(surv))
                        qsort(&kv_A(surv,0), kv_size(surv), sizeof(ani_row_t), cmp_ani_desc);

                    kstring_t ks_out = (kstring_t){0,0,0};
                    format_rows_to_kstr(qry, ref, qn, &surv, ani_opt, &ks_out);

                    ks_arr[qn] = ks_out; // move: keep buffer for print task //
                    kv_destroy(surv);
                }

                // 2) print task: serialized by order[] but can run as soon as qn computed //
                #pragma omp task firstprivate(qn) shared(outfp, ks_arr, comp, order) depend(in: comp[qn]) depend(in: order[qn]) depend(out: order[qn+1])
                {
                    kstring_t *ksp = &ks_arr[qn];
                    if (ksp->l)
                        fwrite(ksp->s, 1, ksp->l, outfp);
                    free(ksp->s);
                    ksp->s = NULL; ksp->l = ksp->m = 0;
                }
            } // for qn //
        } // single //
    } // parallel //

    free(ks_arr);
    free(order);
    free(comp);
#endif // PRINT_UNORDERED //

    if (outfp != stdout) fclose(outfp);
}
*/

void comb_sortedsketch64Xcomb_sortedsketch64(ani_opt_t *ani_opt)
{
	unify_sketch_t *qry_result = generic_sketch_parse(ani_opt->qrydir, ani_query_parse_flags(ani_opt));
	unify_sketch_t *ref_result = generic_sketch_parse(ani_opt->refdir, ani_ref_parse_flags(ani_opt));
	load_infile_meta_for_best_guard(qry_result, ani_opt->qrydir, ani_opt);
	load_infile_meta_for_best_guard(ref_result, ani_opt->refdir, ani_opt);
	const bool same_sketch = strcmp(ani_opt->qrydir, ani_opt->refdir) == 0;

	FILE *outfp = ani_opt->outf[0] == '\0' ? stdout : fopen(ani_opt->outf, "w");
	if (outfp == NULL)
		err(errno, "%s", ani_opt->outf);
	if (ani_opt->fmt == 1)
	{
		for (uint32_t rn = 0; rn < ref_result->infile_num; rn++)
			fprintf(outfp, "\t%s", ref_result->gname[rn]);
		fprintf(outfp, "\n");
		for (uint32_t qn = 0; qn < qry_result->infile_num; qn++)
		{
			fprintf(outfp, "%s", qry_result->gname[qn]);
			for (uint32_t rn = 0; rn < ref_result->infile_num; rn++)
			{
				double metric;
				if (same_sketch && qn == rn)
					metric = ani_matrix_diagonal_value(ani_opt);
				else if (same_sketch)
				{
					const uint32_t q = qn > rn ? qn : rn;
					const uint32_t r = qn > rn ? rn : qn;
					metric = ani_matrix_pair_value(qry_result, q, ref_result, r, ani_opt);
				}
				else
					metric = ani_matrix_pair_value(qry_result, qn, ref_result, rn, ani_opt);
				fprintf(outfp, "\t%lf", metric);
			}
			fprintf(outfp, "\n");
		}
		if (outfp != stdout)
			fclose(outfp);
		free_unify_sketch(qry_result);
		free_unify_sketch(ref_result);
		return;
	}
	if (ani_opt->fmt == 2)
	{
		if (!same_sketch)
			errx(EXIT_FAILURE, "%s(): triangle output requires one sketch or identical -r/-q sketches", __func__);
		for (uint32_t qn = 0; qn < qry_result->infile_num; qn++)
		{
			fprintf(outfp, "%s", qry_result->gname[qn]);
			for (uint32_t rn = 0; rn < qn; rn++)
				fprintf(outfp, "\t%lf", ani_matrix_pair_value(qry_result, qn, ref_result, rn, ani_opt));
			if (ani_opt->d)
				fprintf(outfp, "\t%lf", ani_matrix_diagonal_value(ani_opt));
			fprintf(outfp, "\n");
		}
		if (outfp != stdout)
			fclose(outfp);
		free_unify_sketch(qry_result);
		free_unify_sketch(ref_result);
		return;
	}

	print_ani_detail_header(outfp, ani_opt, false);
	// #pragma omp parallel for num_threads(32) schedule(guided)
	for (uint32_t rn = 0; rn < ref_result->infile_num; rn++)
	{

		uint64_t *arr_ref = ref_result->comb_sketch + ref_result->sketch_index[rn];
		size_t len_ref = ref_result->sketch_index[rn + 1] - ref_result->sketch_index[rn];
		// printf(outfp,"%s", ref_result->gname[rn]);

		for (uint32_t qn = 0; qn < qry_result->infile_num; qn++)
		{
			ani_features_t ani_features;
			uint64_t *arr_qry = qry_result->comb_sketch + qry_result->sketch_index[qn];
			size_t len_qry = qry_result->sketch_index[qn + 1] - qry_result->sketch_index[qn];
			get_ani_features_from_two_sorted_ctxobj64(arr_ref, len_ref, arr_qry, len_qry, &ani_features);
			double af_qry = (double)ani_features.XnY_ctx / len_qry;
			double af_ref = (double)ani_features.XnY_ctx / len_ref;

			if (!ani_report_af_pass(ani_opt, af_qry, af_ref))
				continue;
			ani_features.X_ctx = len_qry;
			double blastn_af_qry = lm3ways_af_ANIb_from_features(&ani_features);
			ani_features.X_ctx = len_ref;
			double blastn_af_ref = lm3ways_af_ANIb_from_features(&ani_features);
			ani_row_t outrow = make_selected_output_row(rn, &ani_features, ani_opt,
														(uint32_t)len_qry,
														(uint32_t)len_ref,
														af_qry, blastn_af_qry,
														af_ref, blastn_af_ref,
														infile_meta_at(qry_result, qn),
														infile_meta_at(ref_result, rn));
			if (outrow.selected_ani < ani_opt->anicut)
				continue;
			print_unified_detail_row(outfp, ani_opt, qry_result->gname[qn], ref_result->gname[rn],
									 &outrow, unify_annotation_at(ref_result, rn));
		}
	}
	if (outfp != stdout)
		fclose(outfp);
	free_unify_sketch(qry_result);
	free_unify_sketch(ref_result);
}

void comb_sortedsketch64_self_matrix(ani_opt_t *ani_opt)
{
	unify_sketch_t *sketch = generic_sketch_parse(ani_opt->qrydir, SKETCH_PARSE_NONE);
	load_infile_meta_for_best_guard(sketch, ani_opt->qrydir, ani_opt);

	FILE *outfp = ani_opt->outf[0] == '\0' ? stdout : fopen(ani_opt->outf, "w");
	if (outfp == NULL)
		err(errno, "%s", ani_opt->outf);

	if (ani_opt->fmt == 1)
	{
		for (uint32_t rn = 0; rn < sketch->infile_num; rn++)
			fprintf(outfp, "\t%s", sketch->gname[rn]);
		fprintf(outfp, "\n");
		for (uint32_t qn = 0; qn < sketch->infile_num; qn++)
		{
			fprintf(outfp, "%s", sketch->gname[qn]);
			for (uint32_t rn = 0; rn < sketch->infile_num; rn++)
			{
				double metric;
				if (qn == rn)
					metric = ani_matrix_diagonal_value(ani_opt);
				else
				{
					const uint32_t q = qn > rn ? qn : rn;
					const uint32_t r = qn > rn ? rn : qn;
					metric = ani_matrix_pair_value(sketch, q, sketch, r, ani_opt);
				}
				fprintf(outfp, "\t%lf", metric);
			}
			fprintf(outfp, "\n");
		}
	}
	else if (ani_opt->fmt == 2)
	{
		for (uint32_t qn = 0; qn < sketch->infile_num; qn++)
		{
			fprintf(outfp, "%s", sketch->gname[qn]);
			for (uint32_t rn = 0; rn < qn; rn++)
				fprintf(outfp, "\t%lf", ani_matrix_pair_value(sketch, qn, sketch, rn, ani_opt));
			if (ani_opt->d)
				fprintf(outfp, "\t%lf", ani_matrix_diagonal_value(ani_opt));
			fprintf(outfp, "\n");
		}
	}
	else
		errx(EXIT_FAILURE, "%s(): one-sketch ANI requires -m1 full matrix or -m2 triangle", __func__);

	if (outfp != stdout)
		fclose(outfp);
	free_unify_sketch(sketch);
}

void check_comb_sortedsketch64(unify_sketch_t *result)
{
	for (uint32_t rn = 0; rn < result->infile_num; rn++)
	{
		uint64_t *arr = result->comb_sketch + result->sketch_index[rn];
		size_t len = result->sketch_index[rn + 1] - result->sketch_index[rn];
		for (uint32_t i = 1; i < len; i++)
		{
			if (arr[i] < arr[i - 1])
				err(EXIT_FAILURE, "%s(): %dth genome %dth kmer < %dth kmer (%lx<%lx)", __func__, rn, i, i - 1, arr[i], arr[i - 1]);
		}
	}
}

//

/**
 * Finds the first occurrence index in large sorted array `b` for each element in a small sorted array `a`.
 *
 * @param a       Sorted array of uint64_t elements (ascending order)
 * @param a_size  Number of elements in array `a`
 * @param b       Sorted array of uint64_t elements (ascending order)
 * @param b_size  Number of elements in array `b`
 * @return        Array of indices (size_t*) where indices[i] = first occurrence of a[i] in b,
 *                SIZE_MAX if not found. Caller must free() the returned array.
 */
size_t *find_first_occurrences_AT_ctxgidobj_arr(const uint64_t *a, size_t a_size,
												const ctxgidobj_t *b, size_t b_size)
{
	size_t *indices = malloc(a_size * sizeof(size_t));
	if (!indices)
		return NULL;
	size_t low = 0; // Track lower bound for binary search
	int nobjbits = Bitslen.obj;

	for (size_t i = 0; i < a_size; ++i)
	{
		const uint64_t target = a[i] >> nobjbits; // (2*(2*holen + iolen ));
		// when conflict objects are kept, skip searching if target[i+1] == target[i].
		if (i > 0 && target == a[i - 1] >> nobjbits)
		{
			indices[i] = indices[i - 1]; // Use the previous index if the current target is the same
			continue;
		}

		size_t high = b_size;
		// Leftmost binary search within [low, high)
		while (low < high)
		{
			size_t mid = low + (high - low) / 2;
			if ((b[mid].ctxgid >> GID_NBITS) < target)
			{
				low = mid + 1;
			}
			else
			{
				high = mid;
			}
		}

		// Check if target was found: modified from: if (low < b_size && b[low] == target) {
		if (low < b_size && (b[low].ctxgid >> GID_NBITS == target))
			indices[i] = low;
		else
			indices[i] = SIZE_MAX; // Not found
	}

	return indices;
}

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
	FILE *outfp, ani_opt_t *ani_opt, int matrix_mode)
{
	ani_features_t ani_features;

	for (int i = 0; i < this_block_size; i++)
	{
		int qry_gid = qry_gid_offset + i;
		int qry_sketch_size = qry_ctx_count[qry_gid];

		if (matrix_mode)
		{
			fprintf(outfp, "%s", qryfname[qry_gid]);
		}

		int loop_j = matrix_mode ? ref_infile_num : num_passid_block[i];
		for (int n = 0; n < loop_j; n++)
		{
			int j = matrix_mode ? n : sort_idani_block[i][n].id;
			ani_features.XnY_ctx = MCTX(ref_infile_num, i, j).num_ctx;
			ani_features.N_diff_obj_section = MOBJ(ref_infile_num, i, j).diff_obj_section;
			ani_features.N_mut2_ctx = MCTX(ref_infile_num, i, j).num_mut2_ctx;
			ani_features.N_diff_obj = MOBJ(ref_infile_num, i, j).diff_obj;

			int ref_sketch_size = ref_ctx_count[j];
				if (qry_sketch_size == 0 || ref_sketch_size == 0)
				{
					if (matrix_mode)
					{
						double metric = ani_opt->s < 0 ? 1.0 - ani_opt->e : ani_opt->e;
						fprintf(outfp, "\t%lf", metric);
					}
					continue;
				}
			float af_qry = (float)ani_features.XnY_ctx / qry_sketch_size;
			float af_ref = (float)ani_features.XnY_ctx / ref_sketch_size;
			ani_features.X_ctx = qry_sketch_size;
			float blastn_af_qry = lm3ways_af_ANIb_from_features(&ani_features);
			ani_features.X_ctx = ref_sketch_size;
			float blastn_af_ref = lm3ways_af_ANIb_from_features(&ani_features);

				ani_row_t outrow = make_selected_output_row((uint32_t)j, &ani_features, ani_opt,
															(uint32_t)qry_sketch_size,
															(uint32_t)ref_sketch_size,
															af_qry, blastn_af_qry,
															af_ref, blastn_af_ref,
															qry_infile_meta ? &qry_infile_meta[qry_gid] : NULL,
															ref_infile_meta ? &ref_infile_meta[j] : NULL);
				double metric = ani_opt->s < 0 ? outrow.selected_ani : outrow.metric;

				if (matrix_mode)
				{
					fprintf(outfp, "\t%lf", metric);
				}
				else if (outrow.selected_ani >= ani_opt->anicut)
				{
					print_unified_detail_row(outfp, ani_opt, qryfname[qry_gid], refname[j],
											 &outrow, annotation_at(refanno, (uint32_t)j));
				}
		}
		if (matrix_mode)
		{
			fprintf(outfp, "\n");
		}
	}
}

void simple_sortedsketch64Xcomb_sortedsketch64(simple_sketch_t *simple_sketch, infile_tab_t *genomes_infiletab, ani_opt_t *ani_opt)
{
	uint64_t *arr_ref = simple_sketch->comb_sketch;
	if (arr_ref == NULL)
		err(EXIT_FAILURE, "%s(): simple_sketch->comb_sketch is NULL", __func__);

	if (simple_sketch->infile_num != genomes_infiletab->infile_num)
		err(EXIT_FAILURE, "%s(): infile_num mismatch: %d vs %d", __func__, simple_sketch->infile_num, genomes_infiletab->infile_num);
	size_t ref_sketch_size = simple_sketch->sketch_index[1] - simple_sketch->sketch_index[0];
	assert(simple_sketch->infile_num == genomes_infiletab->infile_num);

	FILE *outfp = ani_opt->outf[0] == '\0' ? stdout : fopen(ani_opt->outf, "w");
	if (outfp == NULL)
		err(errno, "%s", ani_opt->outf);
	print_ani_detail_header(outfp, ani_opt, true);

	for (uint32_t qn = 1; qn < simple_sketch->infile_num; qn++)
	{
		ani_features_t ani_features;
		uint64_t *arr_qry = simple_sketch->comb_sketch + simple_sketch->sketch_index[qn];
		size_t qry_sketch_size = simple_sketch->sketch_index[qn + 1] - simple_sketch->sketch_index[qn];
		if (qry_sketch_size == 0 || ref_sketch_size == 0)
			continue;

		get_ani_features_from_two_sorted_ctxobj64(arr_ref, ref_sketch_size, arr_qry, qry_sketch_size, &ani_features);

		float af_qry = (double)ani_features.XnY_ctx / qry_sketch_size;

		float af_ref = (double)ani_features.XnY_ctx / ref_sketch_size;
		ani_features.X_ctx = qry_sketch_size;
		float blastn_af_qry = lm3ways_af_ANIb_from_features(&ani_features);
		ani_features.X_ctx = ref_sketch_size;
		float blastn_af_ref = lm3ways_af_ANIb_from_features(&ani_features);

		if (!ani_report_af_pass(ani_opt, af_qry, af_ref))
			continue;
		ani_row_t outrow = make_selected_output_row(0, &ani_features, ani_opt,
													(uint32_t)qry_sketch_size,
													(uint32_t)ref_sketch_size,
													af_qry, blastn_af_qry,
													af_ref, blastn_af_ref,
													NULL, NULL);
		if (outrow.selected_ani >= ani_opt->anicut)
		{
			print_unified_detail_row(outfp, ani_opt,
									 genomes_infiletab->organized_infile_tab[qn].fpath,
									 genomes_infiletab->organized_infile_tab[0].fpath,
									 &outrow, NULL);
		}
	}
	if (outfp != stdout)
		fclose(outfp);
}

// fencepost search method
#include <limits.h>

static inline uint32_t top_k_bits_u64(uint64_t x, int k)
{
	if (k <= 0)
		return 0u;
	if (k >= 64)
		return (uint32_t)x; /* defensive; we clamp k well below 64 */
	return (uint32_t)(x >> (64 - k));
}

int kssd_choose_k_fenceposts(size_t b_size, size_t a_size)
{
	if (a_size == 0 || b_size == 0)
		return 14;
	double ratio = (double)b_size / (double)a_size;			   /* ≈ m = b/a */
	int k = (ratio > 1.0) ? (int)floor(log2(ratio) + 0.5) : 0; /* ~log2(m) */
	if (k < 12)
		k = 12; /* 4K buckets  (~32 KB of fenceposts) */
	if (k > 20)
		k = 20; /* 1M buckets  (~8  MB of fenceposts, 64-bit size_t) */
	return k;
}

int kssd_build_fenceposts_ctxgid(const ctxgidobj_t *b, size_t b_size, int k, size_t *F)
{
	if (!b || !F || k < 0 || k > 32)
		return -1;
	const size_t buckets = (size_t)1u << k;

	for (size_t t = 0; t <= buckets; ++t)
		F[t] = b_size;

	size_t next = 0;
	for (size_t i = 0; i < b_size; ++i)
	{
		uint64_t key = (uint64_t)b[i].ctxgid >> GID_NBITS; /* effective key */
		uint32_t topk = top_k_bits_u64(key, k);
		while (next <= topk)
			F[next++] = i;
		if (next > buckets)
			break;
	}
	while (next <= buckets)
		F[next++] = b_size;
	return 0;
}

static inline size_t lb_in_bucket_ctxgid(const ctxgidobj_t *b, const size_t *F, int k,
										 uint64_t target_key)
{
	const size_t buckets = (size_t)1u << k;
	uint32_t t = top_k_bits_u64(target_key, k);
	if (t > buckets)
		t = (uint32_t)buckets; /* defensive */

	size_t lo = F[t];
	size_t hi = F[t + 1];

	while (lo < hi)
	{
		size_t mid = lo + ((hi - lo) >> 1);
		uint64_t mid_key = (uint64_t)b[mid].ctxgid >> GID_NBITS;
		if (mid_key < target_key)
			lo = mid + 1;
		else
			hi = mid;
	}
	return lo; /* first position with key >= target_key (may be hi) */
}

size_t *kssd_find_first_occurrences_fenceposts(const uint64_t *a, size_t a_size,
											   const ctxgidobj_t *b, size_t b_size,
											   const size_t *F, int k, unsigned nobjbits)
{
	if (!a || !b || !F)
		return NULL;

	size_t *idx = (size_t *)malloc(a_size * sizeof *idx);
	if (!idx)
		return NULL;

	uint64_t prev_key = UINT64_MAX;
	size_t prev_idx = SIZE_MAX;

	for (size_t i = 0; i < a_size; ++i)
	{
		uint64_t target_key = a[i] >> nobjbits;

		if (i && target_key == prev_key)
		{ /* reuse for duplicates in a */
			idx[i] = prev_idx;
			continue;
		}

		size_t pos = lb_in_bucket_ctxgid(b, F, k, target_key);

		if (pos < b_size && (((uint64_t)b[pos].ctxgid >> GID_NBITS) == target_key))
			idx[i] = pos;
		else
			idx[i] = SIZE_MAX;

		prev_key = target_key;
		prev_idx = idx[i];
	}
	return idx;
}



/* comb_manysmall_sortedsketch64Xcomb_fewlarge_sortedsketch64_filter_and_sort_survivors using hash table  */
//helpers: ctx→(beg,end) hash table for query a

typedef struct {
    uint64_t *keys;   // store (ctx+1); 0 means empty
    uint64_t *vals;   // pack(beg,end)
    uint32_t  mask;   // cap-1 (cap power of 2)
    uint32_t  n_runs;
} ctxrun_ht_t;

static inline uint32_t next_pow2_u32(uint32_t x)
{
    if (x <= 1) return 1;
    --x;
    x |= x >> 1; x |= x >> 2; x |= x >> 4; x |= x >> 8; x |= x >> 16;
    return x + 1;
}

static inline uint64_t pack_be(uint32_t beg, uint32_t end)
{
    return ((uint64_t)beg << 32) | (uint64_t)end;
}
static inline uint32_t unpack_beg(uint64_t be) { return (uint32_t)(be >> 32); }
static inline uint32_t unpack_end(uint64_t be) { return (uint32_t)(be); }

static inline void ctxrun_ht_init_from_n(ctxrun_ht_t *ht, uint32_t n)
{
    // Over-allocate from n (len_qry). Very fast lookups; simple build.
    uint32_t cap = next_pow2_u32(n * 2u);         // ~0.5 load factor if nrun~n
    ht->keys = (uint64_t*)calloc(cap, sizeof(uint64_t));
    ht->vals = (uint64_t*)malloc((size_t)cap * sizeof(uint64_t));
    if (!ht->keys || !ht->vals) err(EXIT_FAILURE, "malloc ctxrun_ht");
    ht->mask = cap - 1;
    ht->n_runs = 0;
}

static inline void ctxrun_ht_destroy(ctxrun_ht_t *ht)
{
    free(ht->keys); free(ht->vals);
    ht->keys = NULL; ht->vals = NULL; ht->mask = 0; ht->n_runs = 0;
}

static inline void ctxrun_ht_put(ctxrun_ht_t *ht, uint64_t ctx, uint64_t val_be)
{
    const uint64_t k = ctx + 1; // reserve 0
    uint32_t pos = (uint32_t)mix64(k) & ht->mask;
    while (ht->keys[pos] && ht->keys[pos] != k) pos = (pos + 1) & ht->mask;
    ht->keys[pos] = k;
    ht->vals[pos] = val_be;
}

static inline int ctxrun_ht_get(const ctxrun_ht_t *ht, uint64_t ctx, uint64_t *out_be)
{
    const uint64_t k = ctx + 1;
    uint32_t pos = (uint32_t)mix64(k) & ht->mask;
    for (;;) {
        const uint64_t cur = ht->keys[pos];
        if (!cur) return 0;
        if (cur == k) { *out_be = ht->vals[pos]; return 1; }
        pos = (pos + 1) & ht->mask;
    }
}

// Build ctx->(beg,end) from sorted a[] in ONE pass.
// Must be called once per query qn; afterwards ht is read-only and thread-safe.
static inline void build_ctxrun_ht_from_sorted_query(const uint64_t *a, uint32_t n, uint8_t nobjbits, ctxrun_ht_t *ht)
{
    ctxrun_ht_init_from_n(ht, n);

    for (uint32_t i = 0; i < n; ) {
        const uint64_t ctx = a[i] >> nobjbits;
        const uint32_t beg = i;
        do { ++i; } while (i < n && (a[i] >> nobjbits) == ctx);
        const uint32_t end = i;
        ctxrun_ht_put(ht, ctx, pack_be(beg, end));
        ht->n_runs++;
    }
}

static inline uint32_t count_ctx_runs_sorted_ctxobj64(const uint64_t *a, uint32_t n, uint8_t nobjbits,
													  bool ignoreconflict)
{
    uint32_t runs = 0;
    for (uint32_t i = 0; i < n; ) {
        const uint64_t ctx = a[i] >> nobjbits;
        const uint32_t begin = i;
        do { ++i; } while (i < n && (a[i] >> nobjbits) == ctx);
        if (!ignoreconflict || i - begin == 1)
            runs++;
    }
    return runs;
}

// Return 1 if enough shared contexts remain for the report AF cutoff, else 0.
// Fills features regardless; caller can compute dist/ani.
static inline int get_features_scan_b_hash_a(
    const uint64_t *a, uint32_t n,
    const uint64_t *b, uint32_t m,
    const ctxrun_ht_t *ht,                    // built from a
	uint32_t need_X,                          // ceil(afcut * min(query ctx count, reference ctx count))
    bool ignore_ref_conflict,
    ani_features_t *f)
{
    const uint8_t  nobjbits = Bitslen.obj;
    const uint64_t objmask  = ((nobjbits == 64) ? ~0ULL : ((1ULL << nobjbits) - 1ULL));

    memset(f, 0, sizeof *f);

    for (uint32_t j = 0; j < m; ) {
        // Early fail: even if all remaining b entries hit, cannot reach need_X
        if ((uint32_t)f->XnY_ctx + (m - j) < need_X)
            return 0;

        const uint64_t ctxB = b[j] >> nobjbits;
        const uint32_t b_begin = j;
        do { ++j; } while (j < m && (b[j] >> nobjbits) == ctxB);
        const uint32_t b_end = j;

        if (ignore_ref_conflict && b_end - b_begin > 1)
            continue;

        uint64_t be;
        if (!ctxrun_ht_get(ht, ctxB, &be))
            continue;

        const uint32_t beg  = unpack_beg(be);
        const uint32_t end  = unpack_end(be);

        f->XnY_ctx++;
        const int min_diff = min_diff_sections_ctxobj64_runs(a, beg, end, b, b_begin, b_end, objmask);

        if (min_diff > 0) {
            f->N_diff_obj++;
            f->N_diff_obj_section += min_diff;
            if (min_diff > 1) f->N_mut2_ctx++;
        }
    }

    return ((uint32_t)f->XnY_ctx >= need_X);
}


void comb_manysmall_sortedsketch64Xcomb_fewlarge_sortedsketch64_filter_and_sort_survivors(ani_opt_t *ani_opt)
{
    unify_sketch_t *qry = generic_sketch_parse(ani_opt->qrydir, ani_query_parse_flags(ani_opt));
    unify_sketch_t *ref = generic_sketch_parse(ani_opt->refdir, ani_ref_parse_flags(ani_opt));
    load_infile_meta_for_best_guard(qry, ani_opt->qrydir, ani_opt);
    load_infile_meta_for_best_guard(ref, ani_opt->refdir, ani_opt);

    const uint32_t Q = qry->infile_num;
    const uint32_t R = ref->infile_num;

    FILE *outfp = (ani_opt->outf[0] == '\0') ? stdout : fopen(ani_opt->outf, "w");
    if (!outfp) err(errno, "%s", ani_opt->outf);
    print_ani_detail_header(outfp, ani_opt, false);

    const int P = (ani_opt->p > 0) ? ani_opt->p : 1;

    // store formatted output per qn for ordered printing
    kstring_t *ks_arr = (kstring_t*)calloc(Q, sizeof(kstring_t));
    if (!ks_arr) err(EXIT_FAILURE, "calloc ks_arr");

    const uint8_t  nobjbits = Bitslen.obj;

    for (uint32_t qn = 0; qn < Q; ++qn) {

        // ---- query a (large) ----
        const uint64_t *a = qry->comb_sketch + qry->sketch_index[qn];
        const uint32_t  n = (uint32_t)(qry->sketch_index[qn + 1] - qry->sketch_index[qn]);

        // Build ctx->(beg,end) once for this qn (key speedup)
        ctxrun_ht_t ht;
        build_ctxrun_ht_from_sorted_query(a, n, nobjbits, &ht);

        // Thread-local survivor buffers
        kv_ani_row_t *tls = (kv_ani_row_t*)calloc((size_t)P, sizeof(*tls));
        if (!tls) err(EXIT_FAILURE, "calloc tls");
        for (int t = 0; t < P; ++t) kv_init(tls[t]);

        #pragma omp parallel num_threads(P)
        {
            const int tid = omp_get_thread_num();
            ani_row_t row;

            #pragma omp for schedule(dynamic, 256)
            for (uint32_t rn = 0; rn < R; ++rn) {

                // ---- ref b (small, conflict-free) ----
                const uint64_t *b = ref->comb_sketch + ref->sketch_index[rn];
                const uint32_t  m = (uint32_t)(ref->sketch_index[rn + 1] - ref->sketch_index[rn]);
                const uint32_t  m_ctx = ref->conflict ? count_ctx_runs_sorted_ctxobj64(b, m, nobjbits, ani_opt->ignoreconflict) : m;
                const uint32_t  n_ctx = qry->conflict ? ht.n_runs : n;
                if (m_ctx == 0 || n_ctx == 0)
                    continue;

				const uint32_t need_X = ani_report_af_needed_ctx(ani_opt, n_ctx, m_ctx);

                ani_features_t f;
                get_features_scan_b_hash_a(a, n, b, m, &ht, need_X, ani_opt->ignoreconflict && ref->conflict, &f);
                if (f.XnY_ctx < need_X) continue;

                const double af_q = (double)f.XnY_ctx / (double)n_ctx;
                const double af_r = (double)f.XnY_ctx / (double)m_ctx;

				if (!ani_report_af_pass(ani_opt, af_q, af_r)) continue;

	                ani_features_t tmp = f;
	                tmp.X_ctx = n_ctx;
	                const double blastn_af_q = lm3ways_af_ANIb_from_features(&tmp);
	                tmp = f;
	                tmp.X_ctx = m_ctx;
	                const double blastn_af_r = lm3ways_af_ANIb_from_features(&tmp);

	                row = make_selected_output_row(rn, &f, ani_opt, n_ctx, m_ctx,
	                                               af_q, blastn_af_q, af_r, blastn_af_r,
	                                               infile_meta_at(qry, qn), infile_meta_at(ref, rn));
	                if (row.selected_ani <= ani_opt->anicut) continue;

	                kv_push(ani_row_t, tls[tid], row);
	            }
        } // omp parallel

        // merge tls to surv
        kv_ani_row_t surv; kv_init(surv);
        for (int t = 0; t < P; ++t) {
            if (kv_size(tls[t])) kv_append_rows(&surv, &tls[t]);
            kv_destroy(tls[t]);
        }
        free(tls);

        ctxrun_ht_destroy(&ht);

        // sort + format
        if (kv_size(surv))
            qsort(&kv_A(surv,0), kv_size(surv), sizeof(ani_row_t), cmp_ani_desc);

        kstring_t ks_out = (kstring_t){0,0,0};
        format_rows_to_kstr(qry, ref, qn, &surv, ani_opt, &ks_out);
        ks_arr[qn] = ks_out;

        kv_destroy(surv);
    }

    // ordered print
    for (uint32_t qn = 0; qn < Q; ++qn) {
        if (ks_arr[qn].l) fwrite(ks_arr[qn].s, 1, ks_arr[qn].l, outfp);
        free(ks_arr[qn].s);
    }
    free(ks_arr);

    if (outfp != stdout) fclose(outfp);
}
