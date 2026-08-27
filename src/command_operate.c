/*command_operate.c extent the legency kssd set functions to support long sketch (*.comblco) */
#include "global_basic.h"
#include "command_ani.h"
#include "command_operate.h"
#include "sketch_inspect.h"
#include "../klib/khash.h"
#include <stdint.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

const char lpan_prefix[] = "lpan"; // uint64_t pan
const char luniq_pan_prefix[] = "luniq_pan";
// common vars
static size_t file_size;
static dim_sketch_stat_t lco_stat_readin, lco_stat_pan, lco_stat_origin;
static struct stat s;
static char outfpath[PATHLEN + 20];
static int ret;
extern const char sorted_comb_ctxgid64obj32[];

static void copy_lsketch_annotations(const char *indir, const char *outdir, int infile_num)
{
	if (infile_num <= 0 || !file_exists_in_folder(indir, sketch_anno_stat))
		return;

	size_t anno_file_size = 0;
	char *anno_path = test_get_fullpath(indir, sketch_anno_stat);
	char (*annotations)[PATHLEN] = read_from_file(anno_path, &anno_file_size);
	free(anno_path);

	const size_t expected_size = (size_t)infile_num * PATHLEN;
	if (anno_file_size != expected_size)
		err(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
			__func__, indir, sketch_anno_stat, anno_file_size, expected_size);

	write_to_file(test_create_fullpath(outdir, sketch_anno_stat), annotations, anno_file_size);
	free_read_from_file(annotations, anno_file_size);
}

void sketch_inspect_print_samples(const char *sketch_path)
{
	void *mem_stat = read_from_file(test_get_fullpath(sketch_path, sketch_stat), &file_size);
	memcpy(&lco_stat_readin, mem_stat, sizeof(lco_stat_readin));
	char (*tmpname)[PATHLEN] = mem_stat + sizeof(dim_sketch_stat_t);
	uint64_t *mem_index = (uint64_t *)read_from_file(test_get_fullpath(sketch_path, idx_sketch_suffix), &file_size);
	for (int i = 0; i < lco_stat_readin.infile_num; i++)
		printf("%lu\t%s\n", mem_index[i + 1] - mem_index[i], tmpname[i]);
	free_all(mem_stat, mem_index, NULL);
}

void print_lco_gnames(set_opt_t *set_opt)
{
	sketch_inspect_print_samples(set_opt->insketchpath);
}

void sketch_inspect_print_content(const char *sketch_path, int show_mode)
{
	FILE *fp = NULL, *fab = NULL, *fpos = NULL;
	uint64_t kmer;
	uint64_t pos;
	uint32_t abun;
	ctxgidobj_t ctxgidobj;
	bool abundance = 0;
	if (show_mode == 1)
	{
		uint64_t *sketch_index = read_from_file(test_get_fullpath(sketch_path, idx_sketch_suffix), &file_size);
		int infile_num = file_size / sizeof(uint64_t) - 1;
		if ((fp = fopen(test_get_fullpath(sketch_path, combined_sketch_suffix), "rb")) == NULL)
			err(EXIT_FAILURE, "%s(): Failed to open file '%s/%s'", __func__, sketch_path, combined_sketch_suffix);

		if (file_exists_in_folder(sketch_path, combined_ab_suffix))
		{
			if ((fab = fopen(test_get_fullpath(sketch_path, combined_ab_suffix), "rb")) == NULL)
				err(EXIT_FAILURE, "%s(): Failed to open file '%s/%s'", __func__, sketch_path, combined_ab_suffix);
			abundance = 1;
		}

		for (int i = 0; i < infile_num; i++)
		{
			for (uint64_t j = sketch_index[i]; j < sketch_index[i + 1]; j++)
			{
				fread(&kmer, sizeof(kmer), 1, fp);
				if (abundance)
				{
					fread(&abun, sizeof(abun), 1, fab);
					printf("%d\t%lx\t%u\n", i, kmer, abun);
				}
				else
					printf("%d\t%lx\n", i, kmer);
			}
		}
		free(sketch_index);
	}
	else if (show_mode == 2)
	{
		if ((fp = fopen(test_get_fullpath(sketch_path, sorted_comb_ctxgid64obj32), "rb")) == NULL)
			err(EXIT_FAILURE, "%s(): Failed to open file '%s/%s'", __func__, sketch_path, sorted_comb_ctxgid64obj32);
		while (!feof(fp))
		{
			fread(&ctxgidobj, sizeof(ctxgidobj), 1, fp);
			printf("%lx\t%x\n", ctxgidobj.ctxgid, ctxgidobj.obj);
		}
	}
	else if (show_mode == 3)
	{
		uint64_t *sketch_index = read_from_file(test_get_fullpath(sketch_path, idx_sketch_suffix), &file_size);
		int infile_num = file_size / sizeof(uint64_t) - 1;
		const uint64_t sketch_entries = sketch_index[infile_num];
		if ((fp = fopen(test_get_fullpath(sketch_path, combined_sketch_suffix), "rb")) == NULL)
			err(EXIT_FAILURE, "%s(): Failed to open file '%s/%s'", __func__, sketch_path, combined_sketch_suffix);
		if ((fpos = fopen(test_get_fullpath(sketch_path, sketch_position_suffix), "rb")) == NULL)
			err(EXIT_FAILURE, "%s(): Failed to open file '%s/%s'; build this sketch with --position",
				__func__, sketch_path, sketch_position_suffix);

		struct stat comb_st;
		struct stat pos_st;
		if (fstat(fileno(fp), &comb_st) != 0)
			err(errno, "%s(): fstat %s/%s", __func__, sketch_path, combined_sketch_suffix);
		if (fstat(fileno(fpos), &pos_st) != 0)
			err(errno, "%s(): fstat %s/%s", __func__, sketch_path, sketch_position_suffix);
		const off_t expected_size = (off_t)(sketch_entries * sizeof(uint64_t));
		if (comb_st.st_size != expected_size)
			err(EINVAL, "%s(): %s/%s has %ld bytes, expected %ld",
				__func__, sketch_path, combined_sketch_suffix,
				(long)comb_st.st_size, (long)expected_size);
		if (pos_st.st_size != expected_size)
			err(EINVAL, "%s(): %s/%s has %ld bytes, expected %ld",
				__func__, sketch_path, sketch_position_suffix,
				(long)pos_st.st_size, (long)expected_size);

		for (int i = 0; i < infile_num; i++)
		{
			for (uint64_t j = sketch_index[i]; j < sketch_index[i + 1]; j++)
			{
				if (fread(&kmer, sizeof(kmer), 1, fp) != 1)
					err(EXIT_FAILURE, "%s(): Failed to read %s", __func__, combined_sketch_suffix);
				if (fread(&pos, sizeof(pos), 1, fpos) != 1)
					err(EXIT_FAILURE, "%s(): Failed to read %s", __func__, sketch_position_suffix);
				printf("%d\t%lx\t%lu\n", i, kmer, pos);
			}
		}
		free(sketch_index);
	}
	else
		err(EXIT_FAILURE, "%s(): only show modes 1, 2 and 3 are supported, show_mode =%d", __func__, show_mode);
	if (fp)
		fclose(fp);
	if (abundance)
		fclose(fab);
	if (fpos)
		fclose(fpos);
}

void show_content(set_opt_t *set_opt)
{
	sketch_inspect_print_content(set_opt->insketchpath, set_opt->show);
}

KHASH_SET_INIT_INT64(kmer_set)

static inline void u64_swap(uint64_t *a, uint64_t *b)
{
	uint64_t tmp = *a;
	*a = *b;
	*b = tmp;
}

static void u64_insertion_sort(uint64_t *arr, size_t n)
{
	for (size_t i = 1; i < n; ++i)
	{
		uint64_t key = arr[i];
		size_t j = i;
		while (j > 0 && key < arr[j - 1])
		{
			arr[j] = arr[j - 1];
			--j;
		}
		arr[j] = key;
	}
}

static size_t u64_partition(uint64_t *arr, size_t n)
{
	uint64_t *left = arr;
	uint64_t *mid = arr + n / 2;
	uint64_t *right = arr + n - 1;
	if (*mid < *left)
		u64_swap(mid, left);
	if (*right < *left)
		u64_swap(right, left);
	if (*mid < *right)
		u64_swap(mid, right);
	const uint64_t pivot = *right;

	size_t i = 0;
	size_t j = n - 1;
	while (1)
	{
		while (arr[i] < pivot)
			++i;
		while (pivot < arr[j])
			--j;
		if (i >= j)
			return j;
		u64_swap(&arr[i], &arr[j]);
		++i;
		--j;
	}
}

static void u64_parallel_quicksort_task(uint64_t *arr, size_t n)
{
	const size_t insertion_threshold = 32;
	const size_t parallel_threshold = (size_t)1 << 20;
	if (n <= insertion_threshold)
	{
		u64_insertion_sort(arr, n);
		return;
	}

	size_t p = u64_partition(arr, n);
#ifdef _OPENMP
#pragma omp task shared(arr) if (p + 1 > parallel_threshold)
#endif
	u64_parallel_quicksort_task(arr, p + 1);
#ifdef _OPENMP
#pragma omp task shared(arr) if (n - p - 1 > parallel_threshold)
#endif
	u64_parallel_quicksort_task(arr + p + 1, n - p - 1);
#ifdef _OPENMP
#pragma omp taskwait
#endif
}

static void sort_uint64_values(uint64_t *arr, size_t n, int threads)
{
#ifdef _OPENMP
	if (threads < 1)
		threads = 1;
#pragma omp parallel num_threads(threads)
	{
#pragma omp single nowait
		u64_parallel_quicksort_task(arr, n);
	}
#else
	u64_parallel_quicksort_task(arr, n);
#endif
}

static size_t compact_union_or_unique_uint64(uint64_t *sorted, size_t n, bool unique_only)
{
	size_t out = 0;
	for (size_t i = 0; i < n;)
	{
		size_t j = i + 1;
		while (j < n && sorted[j] == sorted[i])
			++j;
		if (!unique_only || j == i + 1)
			sorted[out++] = sorted[i];
		i = j;
	}
	return out;
}

static void kh_resize_for_size(khash_t(kmer_set) *h, size_t n)
{
	if (n == 0)
		return;
	const size_t requested = n + n / 3 + 1024;
	if (requested <= (size_t)UINT32_MAX)
		kh_resize(kmer_set, h, (khint_t)requested);
}

static int lsketch_ctx_obj_bits(const dim_sketch_stat_t *stat)
{
	const int ctx_bits = stat->coden_len > 0 ? 4 * stat->coden_len : 4 * stat->hclen;
	const int obj_bits = 2 * stat->klen - ctx_bits;
	if (obj_bits < 0 || obj_bits >= 64)
		errx(EXIT_FAILURE, "invalid context-object layout: klen=%d coden_len=%d hclen=%d gives obj_bits=%d",
			 stat->klen, stat->coden_len, stat->hclen, obj_bits);
	return obj_bits;
}

static uint64_t lsketch_set_key(uint64_t ctxobj, set_key_mode_t key_mode, int obj_bits)
{
	return key_mode == SET_KEY_CTX ? ctxobj >> obj_bits : ctxobj;
}

static size_t *build_u64_highbit_bucket_offsets(const uint64_t *sorted, size_t n, unsigned bits)
{
	const size_t bucket_count = (size_t)1 << bits;
	size_t *offsets = calloc(bucket_count + 1, sizeof(offsets[0]));
	if (!offsets)
		err(errno, "%s(): OOM markerdb lookup buckets", __func__);

	size_t bucket = 0;
	for (size_t i = 0; i < n; ++i)
	{
		const size_t value_bucket = sorted[i] >> (64 - bits);
		while (bucket < value_bucket)
			offsets[++bucket] = i;
	}
	while (bucket < bucket_count)
		offsets[++bucket] = n;
	return offsets;
}

static inline bool u64_bucketed_contains(
	const uint64_t *sorted,
	const size_t *bucket_offsets,
	unsigned bits,
	uint64_t value)
{
	const size_t bucket = value >> (64 - bits);
	size_t lo = bucket_offsets[bucket];
	size_t hi = bucket_offsets[bucket + 1];
	while (lo < hi)
	{
		const size_t mid = lo + (hi - lo) / 2;
		if (sorted[mid] < value)
			lo = mid + 1;
		else
			hi = mid;
	}
	return lo < bucket_offsets[bucket + 1] && sorted[lo] == value;
}

static void pwrite_all(int fd, const void *data, size_t bytes, off_t offset, const char *path)
{
	const char *cursor = data;
	while (bytes > 0)
	{
		const ssize_t written = pwrite(fd, cursor, bytes, offset);
		if (written < 0)
		{
			if (errno == EINTR)
				continue;
			err(errno, "%s(): write %s", __func__, path);
		}
		if (written == 0)
			err(EIO, "%s(): short write %s", __func__, path);
		cursor += written;
		offset += written;
		bytes -= (size_t)written;
	}
}

int lsketch_union(set_opt_t *set_opt)
{ // for both union and uniq union

	void *mem_stat = read_from_file(test_get_fullpath(set_opt->insketchpath, sketch_stat), &file_size);
	const size_t stat_file_size = file_size;
	memcpy(&lco_stat_readin, mem_stat, sizeof(lco_stat_readin));
	if (lco_stat_readin.infile_num == 1)
	{ // no need create
		printf("only 1 sketch, use %s as pan-sketch?(Y/N)\n", set_opt->insketchpath);
		char inpbuff;
		scanf(" %c", &inpbuff);
		if ((inpbuff == 'Y') || (inpbuff == 'y'))
		{
			chdir(set_opt->insketchpath);
			if (rename(combined_sketch_suffix, lpan_prefix) != 0)
				err(errno, "lsketch_union()");
			printf("the union directory: %s created successfully\n", set_opt->insketchpath);
			return 1;
		}
	}
	// union operation
	uint64_t *mem_comblco = (uint64_t *)read_from_file(test_get_fullpath(set_opt->insketchpath, combined_sketch_suffix), &file_size);
	const size_t in_kmer_ct = file_size / sizeof(mem_comblco[0]);
	uint64_t *sorted_kmers = (uint64_t *)malloc(file_size);
	if (!sorted_kmers)
		err(errno, "%s(): OOM copying %s/%s for set operation",
			__func__, set_opt->insketchpath, combined_sketch_suffix);
	memcpy(sorted_kmers, mem_comblco, file_size);
	free_read_from_file(mem_comblco, file_size);
	mem_comblco = NULL;

	sort_uint64_values(sorted_kmers, in_kmer_ct, set_opt->p);

	if (set_opt->operation == 2)
	{ // -u: normal union mode
		const size_t union_ct = compact_union_or_unique_uint64(sorted_kmers, in_kmer_ct, false);
		write_to_file(test_create_fullpath(set_opt->outdir, lpan_prefix),
					  sorted_kmers, union_ct * sizeof(sorted_kmers[0]));
	}
	else if (set_opt->operation == 3)
	{ // -q: uniq union mode
		const size_t unique_ct = compact_union_or_unique_uint64(sorted_kmers, in_kmer_ct, true);
		if (!set_opt->q2markerdb)
		{
			write_to_file(test_create_fullpath(set_opt->outdir, luniq_pan_prefix),
						  sorted_kmers, unique_ct * sizeof(sorted_kmers[0]));
		}
		else
		{ // genereate markerdb directly instead of uniq union
			const size_t unique_bytes = unique_ct * sizeof(sorted_kmers[0]);
			uint64_t *shrunk = realloc(sorted_kmers, unique_bytes);
			if (shrunk || unique_bytes == 0)
				sorted_kmers = shrunk;

			const unsigned bucket_bits = 24;
			size_t *bucket_offsets = build_u64_highbit_bucket_offsets(
				sorted_kmers,
				unique_ct,
				bucket_bits);

			mem_comblco = (uint64_t *)read_from_file(test_get_fullpath(set_opt->insketchpath, combined_sketch_suffix), &file_size);
			const size_t comb_file_size = file_size;
			uint64_t *fco_pos = (uint64_t *)read_from_file(test_get_fullpath(set_opt->insketchpath, idx_sketch_suffix), &file_size);
			uint64_t *post_fco_pos = calloc((lco_stat_readin.infile_num + 1), sizeof(uint64_t));
			if (!post_fco_pos)
				err(errno, "%s(): OOM output index", __func__);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(set_opt->p)
#endif
			for (uint32_t i = 0; i < lco_stat_readin.infile_num; i++)
			{
				uint64_t count = 0;
				for (uint64_t n = fco_pos[i]; n < fco_pos[i + 1]; n++)
				{
					if (u64_bucketed_contains(sorted_kmers, bucket_offsets, bucket_bits, mem_comblco[n]))
						count++;
				}
				post_fco_pos[i + 1] = count;
			}

			for (uint32_t i = 0; i < lco_stat_readin.infile_num; i++)
				post_fco_pos[i + 1] += post_fco_pos[i];
			if (post_fco_pos[lco_stat_readin.infile_num] != unique_ct)
				err(EINVAL, "%s(): markerdb output count %lu != unique count %zu",
					__func__, post_fco_pos[lco_stat_readin.infile_num], unique_ct);

			char *out_comb_path = test_create_fullpath(set_opt->outdir, combined_sketch_suffix);
			int outfd = open(out_comb_path, O_CREAT | O_TRUNC | O_WRONLY, 0644);
			if (outfd < 0)
				err(errno, "%s(): cannot create %s/%s", __func__, set_opt->outdir, combined_sketch_suffix);
			if (ftruncate(outfd, (off_t)unique_bytes) != 0)
				err(errno, "%s(): resize %s", __func__, out_comb_path);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(set_opt->p)
#endif
			for (uint32_t i = 0; i < lco_stat_readin.infile_num; i++)
			{
				const uint64_t sample_count = post_fco_pos[i + 1] - post_fco_pos[i];
				if (sample_count == 0)
					continue;
				uint64_t *sample_markers = malloc(sample_count * sizeof(sample_markers[0]));
				if (!sample_markers)
					err(errno, "%s(): OOM sample marker buffer", __func__);
				uint64_t out = 0;
				for (uint64_t n = fco_pos[i]; n < fco_pos[i + 1]; n++)
				{
					if (u64_bucketed_contains(sorted_kmers, bucket_offsets, bucket_bits, mem_comblco[n]))
						sample_markers[out++] = mem_comblco[n];
				}
				if (out != sample_count)
					err(EINVAL, "%s(): marker count changed for sample %u", __func__, i);
				pwrite_all(
					outfd,
					sample_markers,
					sample_count * sizeof(sample_markers[0]),
					(off_t)(post_fco_pos[i] * sizeof(sample_markers[0])),
					out_comb_path);
				free(sample_markers);
			}
			if (close(outfd) != 0)
				err(errno, "%s(): close %s", __func__, out_comb_path);
			write_to_file(format_string("%s/%s", set_opt->outdir, idx_sketch_suffix), post_fco_pos, (lco_stat_readin.infile_num + 1) * sizeof(post_fco_pos[0]));
			free(out_comb_path);
			free(bucket_offsets);
			free_read_from_file(mem_comblco, comb_file_size);
			mem_comblco = NULL;
			free_all(fco_pos, post_fco_pos, NULL);
		}
	}
	else
		err(EINVAL, "operation value %d neither 2 (-u: union) nor 3 (-q :uniq uion )", set_opt->operation);
	write_to_file(test_create_fullpath(set_opt->outdir, sketch_stat), mem_stat, stat_file_size);
	copy_lsketch_annotations(set_opt->insketchpath, set_opt->outdir, lco_stat_readin.infile_num);
	free(sorted_kmers);
	free_all(mem_stat, NULL);
	return 1;
}

int lsketch_operate(set_opt_t *set_opt)
{
	clock_t start_time = clock();
	void *mem_stat_pan = read_from_file(test_get_fullpath(set_opt->pansketchpath, sketch_stat), &file_size);
	memcpy(&lco_stat_pan, mem_stat_pan, sizeof(lco_stat_pan));
	void *mem_stat_lco = read_from_file(test_get_fullpath(set_opt->insketchpath, sketch_stat), &file_size);
	memcpy(&lco_stat_origin, mem_stat_lco, sizeof(lco_stat_origin));
	if (lco_stat_pan.hash_id != lco_stat_origin.hash_id)
		errx(EXIT_FAILURE, "%s(): %s sketching id %u != %s id %u", __func__, set_opt->pansketchpath, lco_stat_pan.hash_id, set_opt->insketchpath, lco_stat_origin.hash_id);
	int key_obj_bits = 0;
	if (set_opt->key_mode == SET_KEY_CTX)
	{
		if (lco_stat_pan.coden_len <= 0 || lco_stat_origin.coden_len <= 0)
			errx(EXIT_FAILURE, "%s(): --key ctx currently supports only -T coden long sketches", __func__);
		key_obj_bits = lsketch_ctx_obj_bits(&lco_stat_origin);
		const int pan_obj_bits = lsketch_ctx_obj_bits(&lco_stat_pan);
		if (pan_obj_bits != key_obj_bits)
			errx(EXIT_FAILURE, "%s(): --key ctx requires matching object-bit widths: pan=%d input=%d",
				 __func__, pan_obj_bits, key_obj_bits);
	}
	// copy sketch stat file to result sketch
	write_to_file(test_create_fullpath(set_opt->outdir, sketch_stat), mem_stat_lco, file_size);
	copy_lsketch_annotations(set_opt->insketchpath, set_opt->outdir, lco_stat_origin.infile_num);

	char *lco_fpath;
	if (file_exists_in_folder(set_opt->pansketchpath, lpan_prefix))
		lco_fpath = test_get_fullpath(set_opt->pansketchpath, lpan_prefix);
	else if (file_exists_in_folder(set_opt->pansketchpath, luniq_pan_prefix))
		lco_fpath = test_get_fullpath(set_opt->pansketchpath, luniq_pan_prefix);
	else
		err(EXIT_FAILURE, "%s():cannot find %s or %s under %s ", __func__, lpan_prefix, luniq_pan_prefix, set_opt->pansketchpath);

	uint64_t *mem_pan = (uint64_t *)read_from_file(lco_fpath, &file_size);
	khash_t(kmer_set) *h = kh_init(kmer_set);
	size_t kmer_ct = file_size / sizeof(uint64_t);
	kh_resize_for_size(h, kmer_ct);
	for (size_t i = 0; i < kmer_ct; i++)
		kh_put(kmer_set, h, lsketch_set_key(mem_pan[i], set_opt->key_mode, key_obj_bits), &ret);
	// read comblco.index to mem
	uint64_t *fco_pos = (uint64_t *)read_from_file(test_get_fullpath(set_opt->insketchpath, idx_sketch_suffix), &file_size);
	// post operation index	calloc
	uint64_t *post_fco_pos = calloc((lco_stat_origin.infile_num + 1), sizeof(uint64_t));
	// read comblco to mem
	uint64_t *tmp_comblco_mem = (uint64_t *)read_from_file(test_get_fullpath(set_opt->insketchpath, combined_sketch_suffix), &file_size);
	const size_t input_comb_bytes = file_size;
	if (input_comb_bytes % sizeof(tmp_comblco_mem[0]) != 0)
		errx(EXIT_FAILURE, "%s(): %s/%s size %zu is not a multiple of %zu",
			 __func__, set_opt->insketchpath, combined_sketch_suffix,
			 input_comb_bytes, sizeof(tmp_comblco_mem[0]));
	const size_t input_kmer_ct = input_comb_bytes / sizeof(tmp_comblco_mem[0]);
	const bool preserve_abundance = lco_stat_origin.koc;
	uint32_t *tmp_abundance_mem = NULL;
	uint32_t *post_abundance_mem = NULL;
	if (preserve_abundance)
	{
		size_t abundance_file_size = 0;
		char *abundance_path = test_get_fullpath(set_opt->insketchpath, combined_ab_suffix);
		tmp_abundance_mem = (uint32_t *)read_from_file(abundance_path, &abundance_file_size);
		free(abundance_path);
		const size_t expected_abundance_bytes = input_kmer_ct * sizeof(tmp_abundance_mem[0]);
		if (abundance_file_size != expected_abundance_bytes)
			errx(EXIT_FAILURE, "%s(): %s/%s has %zu bytes, expected %zu",
				 __func__, set_opt->insketchpath, combined_ab_suffix,
				 abundance_file_size, expected_abundance_bytes);
		post_abundance_mem = (uint32_t *)malloc(expected_abundance_bytes);
		if (expected_abundance_bytes > 0 && !post_abundance_mem)
			err(errno, "%s(): OOM output abundance buffer", __func__);
	}
	uint64_t *post_comblco_mem = (uint64_t *)malloc(input_comb_bytes);
	if (input_comb_bytes > 0 && !post_comblco_mem)
		err(errno, "%s(): OOM output sketch buffer", __func__);
	size_t post_kmer_ct = 0;

	// sketch operation
	for (uint32_t i = 0; i < lco_stat_origin.infile_num; i++)
	{
		for (uint64_t n = fco_pos[i]; n < fco_pos[i + 1]; n++)
		{
			// make sure set_opt->operation == 0 if subtract, == 1 if intersect
			uint64_t key = lsketch_set_key(tmp_comblco_mem[n], set_opt->key_mode, key_obj_bits);
			if (set_opt->operation == (kh_get(kmer_set, h, key) != kh_end(h)))
			{
				post_comblco_mem[post_kmer_ct++] = tmp_comblco_mem[n];
				if (post_abundance_mem)
					post_abundance_mem[post_kmer_ct - 1] = tmp_abundance_mem[n];
			}
		}
		post_fco_pos[i + 1] = post_kmer_ct;
	}
	kh_destroy(kmer_set, h);
	// write to result comblco
	sprintf(outfpath, "%s/%s", set_opt->outdir, combined_sketch_suffix);
	write_to_file(outfpath, post_comblco_mem, post_kmer_ct * sizeof(post_comblco_mem[0]));
	if (preserve_abundance)
	{
		sprintf(outfpath, "%s/%s", set_opt->outdir, combined_ab_suffix);
		write_to_file(outfpath, post_abundance_mem, post_kmer_ct * sizeof(post_abundance_mem[0]));
	}
	// write index
	sprintf(outfpath, "%s/%s", set_opt->outdir, idx_sketch_suffix);
	write_to_file(outfpath, post_fco_pos, (lco_stat_origin.infile_num + 1) * sizeof(post_fco_pos[0]));

	free_all(mem_stat_pan, mem_stat_lco, lco_fpath, mem_pan, fco_pos, post_fco_pos,
			 tmp_comblco_mem, tmp_abundance_mem, post_comblco_mem, post_abundance_mem, NULL);
	return 1;
}

int lgrouping_genomes(set_opt_t *set_opt)
{ // for sorted lcombco only !!! old unsorted one not suportted

	void *mem_stat = read_from_file(test_get_fullpath(set_opt->insketchpath, sketch_stat), &file_size);
	memcpy(&lco_stat_readin, mem_stat, sizeof(lco_stat_readin));

	compan_t *subset = organize_taxf(set_opt->subsetf);
	if (lco_stat_readin.infile_num != subset->gn)
		err(EXIT_FAILURE, "%s(): %s's genome number %d != %s's line number %d", __func__, set_opt->insketchpath, lco_stat_readin.infile_num, set_opt->subsetf, subset->gn);

	// read index and comblco
	uint64_t *tmp_idx = (uint64_t *)read_from_file(test_get_fullpath(set_opt->insketchpath, idx_sketch_suffix), &file_size);
	uint64_t *mem_comblco = (uint64_t *)read_from_file(test_get_fullpath(set_opt->insketchpath, combined_sketch_suffix), &file_size);

	if (tmp_idx[subset->gn] * sizeof(uint64_t) != file_size)
		err(EXIT_FAILURE, "%s(): %s last(%u) index(%lu) * sizeof(uint64_t) != %s file size (%lu) ",
			__func__, idx_sketch_suffix, subset->gn, tmp_idx[subset->gn], combined_sketch_suffix, file_size);
	// out index and comblco
	char *lcombco_f = test_create_fullpath(set_opt->outdir, combined_sketch_suffix);
	int fd = open(lcombco_f, O_CREAT | O_RDWR, 0644);
	ftruncate(fd, file_size); // allowing much larger grouped_comblco than using malloc
	uint64_t *grouped_comblco = (uint64_t *)mmap(NULL, file_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);

	uint64_t *out_idx = calloc((subset->taxn + 1), sizeof(uint64_t));
	int outfn = 0; // uint64_t grouped_kmer_ct = 0;

	for (int t = 0; t < subset->taxn; t++)
	{
		if (subset->tax[t].taxid == 0)
			continue; // ignore taxid 0
		outfn++;
		uint64_t start_offset = out_idx[outfn] = out_idx[outfn - 1];
		for (int n = 1; n <= subset->tax[t].gids[0]; n++)
		{
			int gid = subset->tax[t].gids[n];
			memcpy(grouped_comblco + start_offset, mem_comblco + tmp_idx[gid], (tmp_idx[gid + 1] - tmp_idx[gid]) * sizeof(grouped_comblco[0]));
			start_offset += (tmp_idx[gid + 1] - tmp_idx[gid]);
		}

		if (subset->tax[t].gids[0] > 1)
		{
			qsort(grouped_comblco + out_idx[outfn], start_offset - out_idx[outfn], sizeof(grouped_comblco[0]), qsort_comparator_uint64);
			size_t len = dedup_sorted_uint64(grouped_comblco + out_idx[outfn], start_offset - out_idx[outfn]);
			out_idx[outfn] += len;
		}
		else
			out_idx[outfn] = start_offset;
	}
	// write grouped kmer and index to result
	// write_to_file(test_create_fullpath(set_opt->outdir,combined_sketch_suffix), grouped_comblco, grouped_kmer_ct*sizeof(grouped_comblco[0]));
	if (msync(grouped_comblco, out_idx[outfn] * sizeof(grouped_comblco[0]), MS_SYNC) == -1)
		err(EXIT_FAILURE, "%s(): msync error", __func__);
	if (ftruncate(fd, out_idx[outfn] * sizeof(grouped_comblco[0])) == -1)
		err(EXIT_FAILURE, "%s(): ftrucate resize failed", __func__);
	if (munmap(grouped_comblco, file_size) == -1)
		err(EXIT_FAILURE, "%s(): munmap", __func__);
	close(fd);

	write_to_file(test_create_fullpath(set_opt->outdir, idx_sketch_suffix), out_idx, (outfn + 1) * sizeof(out_idx[0]));
	// write stat file
	lco_stat_readin.infile_num = outfn;
	lco_stat_readin.koc = 0;
	char (*input_anno)[PATHLEN] = NULL;
	size_t input_anno_size = 0;
	if (file_exists_in_folder(set_opt->insketchpath, sketch_anno_stat))
	{
		char *anno_path = test_get_fullpath(set_opt->insketchpath, sketch_anno_stat);
		input_anno = read_from_file(anno_path, &input_anno_size);
		free(anno_path);
		const size_t expected_size = (size_t)subset->gn * PATHLEN;
		if (input_anno_size != expected_size)
			err(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
				__func__, set_opt->insketchpath, sketch_anno_stat,
				input_anno_size, expected_size);
	}

	char (*tmpfname)[PATHLEN] = calloc(outfn, PATHLEN);
	char (*tmpanno)[PATHLEN] = input_anno ? calloc(outfn, PATHLEN) : NULL;
	if (!tmpfname || (input_anno && !tmpanno))
		err(errno, "%s(): OOM output names/annotations", __func__);
	int idx = 0;
	for (int t = 0; t < subset->taxn; t++)
	{
		if (subset->tax[t].taxid != 0)
		{
			if (subset->tax[t].taxname != NULL)
				snprintf(tmpfname[idx], PATHLEN, "%d_%s", subset->tax[t].taxid, subset->tax[t].taxname);
			else
				snprintf(tmpfname[idx], PATHLEN, "%d", subset->tax[t].taxid);
			if (tmpanno)
			{
				if (subset->tax[t].gids[0] == 1)
					memcpy(tmpanno[idx], input_anno[subset->tax[t].gids[1]], PATHLEN);
				else
					snprintf(tmpanno[idx], PATHLEN, "%s", tmpfname[idx]);
			}
			idx++;
		}
		free_all(subset->tax[t].gids, subset->tax[t].taxname, NULL);
	}
	concat_and_write_to_file(test_create_fullpath(set_opt->outdir, sketch_stat), &lco_stat_readin, sizeof(lco_stat_readin), tmpfname, PATHLEN * outfn);
	if (tmpanno)
		write_to_file(test_create_fullpath(set_opt->outdir, sketch_anno_stat), tmpanno, (size_t)outfn * PATHLEN);
	if (munmap(mem_comblco, file_size) == -1)
	{
		if (mem_comblco != NULL)
			free(mem_comblco);
		else
			err(EXIT_FAILURE, "%s(): munmap", __func__);
	}

	if (input_anno)
		free_read_from_file(input_anno, input_anno_size);
	free_all(mem_stat, tmp_idx, out_idx, subset->tax, subset, tmpfname, tmpanno, NULL);
	return outfn;
}
