#define _GNU_SOURCE   // must come before any system header
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include "global_basic.h"
#include "command_sketch.h"
#include "command_ani.h"
#include "command_progress.h"
#include "pairwise_graph.h"
#include "sketch_rearrange.h"
#include "kssdlib_sort.h"
#include <time.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <limits.h>

#ifndef likely
#define likely(x) __builtin_expect(!!(x), 1)
#endif
#ifndef unlikely
#define unlikely(x) __builtin_expect(!!(x), 0)
#endif

// shared global vars
const char sketch_stat[] = "lcofiles.stat";
const char sketch_qc_stat[] = "lcofiles.qc";
const char sketch_anno_stat[] = "lcofiles.anno";
const char sketch_infile_meta_stat[] = "lcofiles.infilemeta";
const char sketch_position_suffix[] = "comblco.position";
const char sketch_suffix[] = "lco"; // long co, 64bits
const char combined_sketch_suffix[] = "comblco";
const char idx_sketch_suffix[] = "comblco.index";
const char combined_ab_suffix[] = "comblco.a";
const char sorted_comb_ctxgid64obj32[] = "sortedcomb_ctxgid64obj32";
// public vars shared across files
uint32_t FILTER, hash_id;
dim_sketch_stat_t comblco_stat_one, comblco_stat_it;
// tmp container in this scope
static size_t file_size;
static char tmp_fname[PATHLEN + 20];
static struct stat tmpstat;

static void remove_sketch_positions(const char *outdir);
static void sketch_few_files_with_intrafile_parallel_pos(sketch_opt_t *opt, infile_tab_t *tab, int BATCH_READS);

static inline bool sketch_needs_kmer_counts(const sketch_opt_t *opt)
{
    return opt->abundance || opt->kmerocrs > 1 || opt->npercentile > 0.0 || opt->reads_qc;
}

static void radix_sort_u32_inplace(uint32_t *a, size_t n)
{
    if (n < 2) return;
    uint32_t *tmp = (uint32_t*)malloc(n * sizeof(uint32_t));
    if (!tmp) {
        qsort(a, n, sizeof(uint32_t), qsort_comparator_uint32);
        return;
    }

    size_t cnt[256];
    for (unsigned pass = 0; pass < 4; ++pass) {
        for (int i = 0; i < 256; ++i) cnt[i] = 0;
        const unsigned sh = pass * 8;
        for (size_t i = 0; i < n; ++i) ++cnt[(a[i] >> sh) & 0xFFu];
        size_t sum = 0;
        for (int i = 0; i < 256; ++i) {
            const size_t c = cnt[i];
            cnt[i] = sum;
            sum += c;
        }
        for (size_t i = 0; i < n; ++i) {
            const unsigned b = (unsigned)((a[i] >> sh) & 0xFFu);
            tmp[cnt[b]++] = a[i];
        }
        uint32_t *swap = a;
        a = tmp;
        tmp = swap;
    }

    free(tmp);
}

static uint32_t kmer_count_percentile_threshold(const SortedKV_Arrays_t *kv, double percentile,
                                                uint32_t base_cutoff)
{
    if (percentile <= 0.0 || kv->len == 0 || !kv->values)
        return base_cutoff;

    uint32_t *counts = (uint32_t*)malloc(kv->len * sizeof(uint32_t));
    if (!counts)
        err(errno, "%s(): OOM count percentile buffer", __func__);

    size_t selected = 0;
    uint64_t total_weight = 0;
    for (size_t i = 0; i < kv->len; ++i) {
        if (kv->values[i] >= base_cutoff) {
            const uint32_t count = kv->values[i];
            counts[selected++] = count;
            if (UINT64_MAX - total_weight < count)
                total_weight = UINT64_MAX;
            else
                total_weight += count;
        }
    }
    if (selected == 0) {
        free(counts);
        return base_cutoff;
    }

    radix_sort_u32_inplace(counts, selected);

    long double rank_ld = ceill((long double)percentile * (long double)total_weight);
    uint64_t rank = rank_ld >= (long double)UINT64_MAX ? UINT64_MAX : (uint64_t)rank_ld;
    if (rank < 1) rank = 1;

    uint64_t cumulative = 0;
    uint32_t threshold = counts[selected - 1];
    for (size_t i = 0; i < selected; ++i) {
        const uint32_t count = counts[i];
        if (UINT64_MAX - cumulative < count)
            cumulative = UINT64_MAX;
        else
            cumulative += count;
        if (cumulative >= rank) {
            threshold = count;
            break;
        }
    }
    free(counts);
    return threshold;
}

typedef struct kmer_count_filter_result
{
    uint32_t raw_cutoff;
    uint32_t effective_cutoff;
    uint32_t upper_cutoff;
    uint32_t reads_qc_mode;
    dim_sketch_qc_stat_t sample_qc;
    bool used_reads_qc;
} kmer_count_filter_result_t;

typedef struct count_histogram
{
    uint32_t *counts;
    uint64_t *freqs;
    size_t len;
    uint64_t total;
} count_histogram_t;

typedef struct reads_qc_range
{
    bool valid;
    uint32_t lower;
    uint32_t upper;
    uint32_t mode;
    long double mean;
    long double variance;
    long double nb_size;
} reads_qc_range_t;

typedef struct reads_qc_peak
{
    size_t mode_idx;
    size_t lower_limit_idx;
    size_t valley_idx;
    uint32_t weighted_q25;
    uint32_t weighted_median;
    uint32_t weighted_q75;
    bool has_error_valley;
    bool used_weighted_anchor;
} reads_qc_peak_t;

static inline bool sketch_reports_count_filter(const sketch_opt_t *opt)
{
    return opt->npercentile > 0.0 || opt->reads_qc;
}

static inline bool sketch_records_sample_qc(const sketch_opt_t *opt)
{
    return opt->abundance || opt->reads_qc;
}

static inline dim_sketch_qc_stat_t sample_qc_from_range(const reads_qc_range_t *range,
                                                        bool applied)
{
    dim_sketch_qc_stat_t qc = {0};
    if (!range || !range->valid)
        return qc;

    qc.flags = DIM_SKETCH_QC_RANGE_VALID;
    if (applied)
        qc.flags |= DIM_SKETCH_QC_RANGE_APPLIED;
    qc.reads_qc_lower = range->lower;
    qc.reads_qc_upper = range->upper;
    qc.reads_qc_mode = range->mode;
    return qc;
}

static long double histogram_smooth_freq(const count_histogram_t *hist, size_t idx);

static void free_count_histogram(count_histogram_t *hist)
{
    free(hist->counts);
    free(hist->freqs);
    *hist = (count_histogram_t){0};
}

static bool reads_qc_debug_enabled(void)
{
    const char *value = getenv("KSSD3_READSQC_DEBUG");
    return value && value[0] != '\0' && strcmp(value, "0") != 0;
}

static void debug_print_count_histogram(const count_histogram_t *hist, const reads_qc_peak_t *peak)
{
    if (!reads_qc_debug_enabled())
        return;

    fprintf(stderr, "\n--readsQC histogram bins=%zu total_distinct=%lu selected_mode=%u",
            hist->len, (unsigned long)hist->total,
            hist->len ? hist->counts[peak->mode_idx] : 0);
    if (peak->weighted_median) {
        fprintf(stderr, " weighted_q25=%u weighted_median=%u weighted_q75=%u",
                peak->weighted_q25, peak->weighted_median, peak->weighted_q75);
    }
    if (peak->has_error_valley) {
        fprintf(stderr, " valley=%u lower_limit=%u",
                hist->counts[peak->valley_idx], hist->counts[peak->lower_limit_idx]);
    }
    fputc('\n', stderr);

    const size_t max_bins = hist->len < 160 ? hist->len : 160;
    for (size_t i = 0; i < max_bins; ++i) {
        fprintf(stderr, "--readsQC hist count=%u distinct=%lu smooth=%.0Lf%s%s\n",
                hist->counts[i], (unsigned long)hist->freqs[i],
                histogram_smooth_freq(hist, i),
                i == peak->mode_idx ? " mode" : "",
                peak->has_error_valley && i == peak->valley_idx ? " valley" : "");
    }
    if (hist->len > max_bins)
        fprintf(stderr, "--readsQC hist truncated after %zu bins\n", max_bins);
}

static count_histogram_t build_count_histogram(const SortedKV_Arrays_t *kv, uint32_t min_count)
{
    count_histogram_t hist = {0};
    if (!kv->len || !kv->values)
        return hist;

    uint32_t *counts = (uint32_t *)malloc(kv->len * sizeof(uint32_t));
    if (!counts)
        err(errno, "%s(): OOM count histogram buffer", __func__);

    size_t selected = 0;
    for (size_t i = 0; i < kv->len; ++i) {
        if (kv->values[i] >= min_count)
            counts[selected++] = kv->values[i];
    }

    if (!selected) {
        free(counts);
        return hist;
    }

    radix_sort_u32_inplace(counts, selected);

    uint64_t *freqs = (uint64_t *)malloc(selected * sizeof(uint64_t));
    if (!freqs)
        err(errno, "%s(): OOM count histogram frequencies", __func__);

    size_t bins = 0;
    counts[bins] = counts[0];
    freqs[bins] = 1;
    for (size_t i = 1; i < selected; ++i) {
        if (counts[i] == counts[bins]) {
            ++freqs[bins];
        } else {
            ++bins;
            counts[bins] = counts[i];
            freqs[bins] = 1;
        }
    }
    ++bins;

    hist.counts = counts;
    hist.freqs = freqs;
    hist.len = bins;
    hist.total = (uint64_t)selected;
    return hist;
}

static bool histogram_weighted_stats(const count_histogram_t *hist, size_t first, size_t last,
                                     uint64_t *mass_out, long double *mean_out,
                                     long double *variance_out)
{
    if (!hist->len || first > last || last >= hist->len)
        return false;

    uint64_t mass = 0;
    long double sum = 0.0L;
    long double sumsq = 0.0L;
    for (size_t i = first; i <= last; ++i) {
        const long double freq = (long double)hist->freqs[i];
        const long double count = (long double)hist->counts[i];
        mass += hist->freqs[i];
        sum += count * freq;
        sumsq += count * count * freq;
    }

    if (!mass)
        return false;

    const long double mean = sum / (long double)mass;
    long double variance = sumsq / (long double)mass - mean * mean;
    if (variance < 0.0L)
        variance = 0.0L;

    if (mass_out) *mass_out = mass;
    if (mean_out) *mean_out = mean;
    if (variance_out) *variance_out = variance;
    return true;
}

static long double negative_binomial_size(long double mean, long double variance)
{
    if (mean <= 0.0L || variance <= mean)
        return 0.0L;
    return (mean * mean) / (variance - mean);
}

static long double log_count_pmf(uint32_t count, long double mean, long double nb_size)
{
    const long double x = (long double)count;
    if (mean <= 0.0L)
        return -INFINITY;
    if (nb_size <= 0.0L) {
        return -mean + x * logl(mean) - lgammal(x + 1.0L);
    }

    const long double p = nb_size / (nb_size + mean);
    return lgammal(x + nb_size) - lgammal(nb_size) - lgammal(x + 1.0L)
           + nb_size * logl(p) + x * logl(1.0L - p);
}

static long double histogram_smooth_freq(const count_histogram_t *hist, size_t idx)
{
    if (!hist->len)
        return 0.0L;
    if (hist->len == 1)
        return (long double)hist->freqs[0];

    long double total = 2.0L * (long double)hist->freqs[idx];
    long double weight = 2.0L;
    if (idx > 0 && hist->counts[idx - 1] + 1 == hist->counts[idx]) {
        total += (long double)hist->freqs[idx - 1];
        weight += 1.0L;
    }
    if (idx + 1 < hist->len && hist->counts[idx] + 1 == hist->counts[idx + 1]) {
        total += (long double)hist->freqs[idx + 1];
        weight += 1.0L;
    }
    return total / weight;
}

static bool histogram_local_minimum(const count_histogram_t *hist, size_t idx)
{
    if (idx == 0 || idx + 1 >= hist->len)
        return false;
    if (hist->counts[idx - 1] + 1 != hist->counts[idx]
        || hist->counts[idx] + 1 != hist->counts[idx + 1])
        return false;

    const long double left = histogram_smooth_freq(hist, idx - 1);
    const long double mid = histogram_smooth_freq(hist, idx);
    const long double right = histogram_smooth_freq(hist, idx + 1);
    return mid <= left && mid <= right;
}

static bool histogram_local_maximum(const count_histogram_t *hist, size_t idx)
{
    const long double mid = histogram_smooth_freq(hist, idx);
    const long double left = idx > 0 ? histogram_smooth_freq(hist, idx - 1) : -1.0L;
    const long double right = idx + 1 < hist->len ? histogram_smooth_freq(hist, idx + 1) : -1.0L;

    if (idx > 0 && hist->counts[idx - 1] + 1 != hist->counts[idx])
        return false;
    if (idx + 1 < hist->len && hist->counts[idx] + 1 != hist->counts[idx + 1])
        return false;
    return mid >= left && mid >= right;
}

static size_t histogram_best_peak_between(const count_histogram_t *hist, size_t first, size_t last)
{
    size_t best = first;
    long double best_freq = -1.0L;
    for (size_t i = first; i <= last && i < hist->len; ++i) {
        const long double freq = histogram_smooth_freq(hist, i);
        if (freq > best_freq) {
            best_freq = freq;
            best = i;
        }
    }
    return best;
}

static uint32_t histogram_weighted_quantile_count(const count_histogram_t *hist, long double quantile)
{
    if (!hist->len)
        return 0;

    long double total = 0.0L;
    for (size_t i = 0; i < hist->len; ++i)
        total += (long double)hist->counts[i] * (long double)hist->freqs[i];

    if (total <= 0.0L)
        return hist->counts[0];

    const long double target = total * quantile;
    long double cumulative = 0.0L;
    for (size_t i = 0; i < hist->len; ++i) {
        cumulative += (long double)hist->counts[i] * (long double)hist->freqs[i];
        if (cumulative >= target)
            return hist->counts[i];
    }

    return hist->counts[hist->len - 1];
}

static size_t histogram_first_count_at_least(const count_histogram_t *hist, uint32_t count)
{
    for (size_t i = 0; i < hist->len; ++i) {
        if (hist->counts[i] >= count)
            return i;
    }
    return hist->len - 1;
}

static size_t histogram_last_count_at_most(const count_histogram_t *hist, uint32_t count)
{
    for (size_t i = hist->len; i > 0; --i) {
        if (hist->counts[i - 1] <= count)
            return i - 1;
    }
    return 0;
}

static size_t histogram_lowest_smooth_between(const count_histogram_t *hist, size_t first, size_t last)
{
    size_t best = first;
    long double best_freq = histogram_smooth_freq(hist, first);
    for (size_t i = first + 1; i <= last && i < hist->len; ++i) {
        const long double freq = histogram_smooth_freq(hist, i);
        if (freq < best_freq) {
            best_freq = freq;
            best = i;
        }
    }
    return best;
}

static bool select_weighted_anchor_peak(const count_histogram_t *hist, reads_qc_peak_t *peak)
{
    peak->weighted_q25 = histogram_weighted_quantile_count(hist, 0.25L);
    peak->weighted_median = histogram_weighted_quantile_count(hist, 0.50L);
    peak->weighted_q75 = histogram_weighted_quantile_count(hist, 0.75L);

    if (!peak->weighted_median || peak->weighted_q25 > peak->weighted_q75)
        return false;

    size_t first = histogram_first_count_at_least(hist, peak->weighted_q25);
    size_t last = histogram_last_count_at_most(hist, peak->weighted_q75);
    if (first > last)
        return false;

    size_t best = first;
    long double best_freq = -1.0L;
    bool found_local_peak = false;
    for (size_t i = first; i <= last && i < hist->len; ++i) {
        if (!histogram_local_maximum(hist, i))
            continue;
        const long double freq = histogram_smooth_freq(hist, i);
        if (freq > best_freq) {
            best_freq = freq;
            best = i;
            found_local_peak = true;
        }
    }

    if (!found_local_peak)
        best = histogram_best_peak_between(hist, first, last);

    peak->mode_idx = best;
    peak->lower_limit_idx = 0;
    peak->valley_idx = 0;
    peak->has_error_valley = false;
    peak->used_weighted_anchor = true;

    if (best > 1) {
        const size_t valley = histogram_lowest_smooth_between(hist, 0, best - 1);
        if (valley + 1 <= best) {
            peak->valley_idx = valley;
            peak->lower_limit_idx = valley + 1;
            peak->has_error_valley = true;
        }
    }

    return true;
}

static reads_qc_peak_t select_reads_qc_peak(const count_histogram_t *hist)
{
    reads_qc_peak_t peak = {0};
    if (!hist->len)
        return peak;

    peak.mode_idx = histogram_best_peak_between(hist, 0, hist->len - 1);
    peak.lower_limit_idx = 0;

    if (select_weighted_anchor_peak(hist, &peak))
        return peak;

    const long double global_peak_freq = histogram_smooth_freq(hist, peak.mode_idx);
    long double best_score = -1.0L;

    for (size_t valley = 1; valley + 1 < hist->len; ++valley) {
        if (!histogram_local_minimum(hist, valley))
            continue;

        const size_t right_peak = histogram_best_peak_between(hist, valley + 1, hist->len - 1);
        if (right_peak <= valley)
            continue;

        const long double valley_freq = histogram_smooth_freq(hist, valley);
        const long double right_peak_freq = histogram_smooth_freq(hist, right_peak);

        if (right_peak_freq < 25.0L)
            continue;
        if (right_peak_freq < global_peak_freq * 0.001L)
            continue;
        if (right_peak_freq < valley_freq * 1.5L)
            continue;
        if (!histogram_local_maximum(hist, right_peak) && right_peak + 1 < hist->len)
            continue;

        const long double score = right_peak_freq * log1pl((long double)hist->counts[right_peak]);
        if (score > best_score) {
            best_score = score;
            peak.mode_idx = right_peak;
            peak.lower_limit_idx = valley + 1;
            peak.valley_idx = valley;
            peak.has_error_valley = true;
        }
    }

    return peak;
}

static bool histogram_bin_is_nb_like(const count_histogram_t *hist, size_t idx, size_t mode_idx,
                                     long double mean, long double nb_size)
{
    const long double observed_ratio =
        (long double)hist->freqs[idx] / (long double)hist->freqs[mode_idx];
    const long double log_expected_ratio =
        log_count_pmf(hist->counts[idx], mean, nb_size)
        - log_count_pmf(hist->counts[mode_idx], mean, nb_size);

    if (!isfinite((double)log_expected_ratio))
        return false;

    const long double expected_ratio = expl(log_expected_ratio);
    if (expected_ratio < 0.005L && observed_ratio < 0.005L)
        return false;

    const long double tolerance = 3.0L;
    return observed_ratio >= expected_ratio / tolerance
           && observed_ratio <= expected_ratio * tolerance;
}

static reads_qc_range_t infer_reads_qc_count_range(const SortedKV_Arrays_t *kv, uint32_t min_count)
{
    reads_qc_range_t range = {0};
    count_histogram_t hist = build_count_histogram(kv, min_count);
    if (!hist.len) {
        free_count_histogram(&hist);
        return range;
    }

    const reads_qc_peak_t peak_info = select_reads_qc_peak(&hist);
    debug_print_count_histogram(&hist, &peak_info);
    const size_t mode_idx = peak_info.mode_idx;
    const size_t lower_limit_idx = peak_info.lower_limit_idx;

    const uint64_t peak = hist.freqs[mode_idx];
    const uint64_t min_fit_freq = peak / 20 > 2 ? peak / 20 : 2;

    size_t fit_left = mode_idx;
    while (fit_left > lower_limit_idx
           && hist.counts[fit_left - 1] + 1 == hist.counts[fit_left]
           && hist.freqs[fit_left - 1] >= min_fit_freq) {
        --fit_left;
    }

    size_t fit_right = mode_idx;
    while (fit_right + 1 < hist.len
           && hist.counts[fit_right] + 1 == hist.counts[fit_right + 1]
           && hist.freqs[fit_right + 1] >= min_fit_freq) {
        ++fit_right;
    }

    while (fit_right - fit_left + 1 < 3 && (fit_left > lower_limit_idx || fit_right + 1 < hist.len)) {
        if (fit_left > lower_limit_idx)
            --fit_left;
        if (fit_right + 1 < hist.len)
            ++fit_right;
    }

    uint64_t mass = 0;
    long double mean = 0.0L, variance = 0.0L;
    if (!histogram_weighted_stats(&hist, fit_left, fit_right, &mass, &mean, &variance)) {
        free_count_histogram(&hist);
        return range;
    }

    const long double nb_size = negative_binomial_size(mean, variance);

    size_t range_left = mode_idx;
    while (range_left > lower_limit_idx && hist.counts[range_left - 1] + 1 == hist.counts[range_left]) {
        if (!histogram_bin_is_nb_like(&hist, range_left - 1, mode_idx, mean, nb_size))
            break;
        --range_left;
    }

    size_t range_right = mode_idx;
    while (range_right + 1 < hist.len && hist.counts[range_right] + 1 == hist.counts[range_right + 1]) {
        if (!histogram_bin_is_nb_like(&hist, range_right + 1, mode_idx, mean, nb_size))
            break;
        ++range_right;
    }

    if (range_left == range_right) {
        range_left = fit_left;
        range_right = fit_right;
    }

    if (peak_info.used_weighted_anchor) {
        const size_t weighted_idx = histogram_first_count_at_least(&hist, peak_info.weighted_median);
        if (weighted_idx < range_left)
            range_left = weighted_idx < lower_limit_idx ? lower_limit_idx : weighted_idx;
        if (weighted_idx > range_right)
            range_right = weighted_idx;
    }

    range.valid = true;
    range.lower = hist.counts[range_left];
    range.upper = hist.counts[range_right];
    range.mode = hist.counts[mode_idx];
    range.mean = mean;
    range.variance = variance;
    range.nb_size = nb_size;

    free_count_histogram(&hist);
    return range;
}

static void filter_count_range_SortedKV_Arrays(SortedKV_Arrays_t *kv, uint32_t lower, uint32_t upper)
{
    if (!kv->len || !kv->values)
        return;

    size_t write_idx = 0;
    for (size_t i = 0; i < kv->len; ++i) {
        const uint32_t count = kv->values[i];
        if (count >= lower && (!upper || count <= upper)) {
            kv->keys[write_idx] = kv->keys[i];
            kv->values[write_idx] = kv->values[i];
            if (kv->positions)
                kv->positions[write_idx] = kv->positions[i];
            ++write_idx;
        }
    }
    kv->len = write_idx;
}

static void print_kmer_count_filter_result(const sketch_opt_t *opt, const char *name,
                                           const kmer_count_filter_result_t *result,
                                           bool leading_newline)
{
    if (leading_newline)
        fputc('\n', stderr);

    if (result->used_reads_qc) {
        fprintf(stderr, "--readsQC count range for %s: mode=%u lower=%u upper=%u",
                name, result->reads_qc_mode, result->effective_cutoff,
                result->upper_cutoff);
        if (opt->ncap > 0)
            fprintf(stderr, " ncap=%u", (uint32_t)opt->ncap);
        fputc('\n', stderr);
    } else if (opt->reads_qc && opt->npercentile <= 0.0) {
        fprintf(stderr, "--readsQC count range for %s: insufficient count histogram; lower=%u\n",
                name, result->effective_cutoff);
    } else if (opt->ncap > 0) {
        fprintf(stderr, "--npercentile %.6g count cutoff for %s: raw=%u capped=%u effective=%u\n",
                opt->npercentile, name, result->raw_cutoff, (uint32_t)opt->ncap,
                result->effective_cutoff);
    } else {
        fprintf(stderr, "--npercentile %.6g count cutoff for %s: %u\n",
                opt->npercentile, name, result->effective_cutoff);
    }
}

static inline kmer_count_filter_result_t apply_kmer_count_filters(SortedKV_Arrays_t *kv,
                                                                  const sketch_opt_t *opt)
{
    if (!kv->len)
        return (kmer_count_filter_result_t){0};

    uint32_t cutoff = opt->kmerocrs > 1 ? (uint32_t)opt->kmerocrs : 1u;
    if (opt->npercentile > 0.0) {
        const uint32_t pct_cutoff = kmer_count_percentile_threshold(kv, opt->npercentile, cutoff);
        if (pct_cutoff > cutoff)
            cutoff = pct_cutoff;
    }

    uint32_t upper = 0;
    uint32_t reads_qc_mode = 0;
    dim_sketch_qc_stat_t sample_qc = {0};
    bool used_reads_qc = false;
    if (sketch_records_sample_qc(opt)) {
        const uint32_t range_min_count = cutoff < 2 ? 2u : cutoff;
        if (opt->reads_qc && cutoff < 2)
            cutoff = 2;

        const reads_qc_range_t range = infer_reads_qc_count_range(kv, range_min_count);
        if (range.valid) {
            reads_qc_mode = range.mode;
            sample_qc = sample_qc_from_range(&range, opt->reads_qc);
            if (opt->reads_qc) {
                if (range.lower > cutoff)
                    cutoff = range.lower;
                upper = range.upper;
                used_reads_qc = true;
            }
        }
    }

    const uint32_t raw_cutoff = cutoff;
    if (opt->ncap > 0 && cutoff > (uint32_t)opt->ncap)
        cutoff = (uint32_t)opt->ncap;

    if (cutoff > 1 || upper)
        filter_count_range_SortedKV_Arrays(kv, cutoff, upper);

    return (kmer_count_filter_result_t){
        .raw_cutoff = raw_cutoff,
        .effective_cutoff = cutoff,
        .upper_cutoff = upper,
        .reads_qc_mode = reads_qc_mode,
        .sample_qc = sample_qc,
        .used_reads_qc = used_reads_qc,
    };
}

uint32_t get_sketching_id(uint32_t hclen, uint32_t holen, uint32_t iolen, uint32_t drfold, uint32_t FILTER)
{
    return GET_SKETCHING_ID(hclen, holen, iolen, drfold, FILTER);
}

void compute_sketch(sketch_opt_t *sketch_opt_val, infile_tab_t *infile_stat)
{
    if (!sketch_opt_val->position)
        remove_sketch_positions(sketch_opt_val->outdir);
    if (sketch_opt_val->split_mfa)
    { // mfa files parse
        mfa2sortedctxobj64_v2(sketch_opt_val, infile_stat);
        return;
    }
    // Decide mode based on #files vs. threads
    if ((!sketch_opt_val->asone) && infile_stat->infile_num >= sketch_opt_val->p)         // Mode-A: per-file parallel
        sketch_many_files_in_parallel(sketch_opt_val, infile_stat, 1024);
    else         // Mode-B: in-file parallel
        sketch_few_files_with_intrafile_parallel(sketch_opt_val, infile_stat, 4096);

    return;

};

void combine_lco(sketch_opt_t *sketch_opt_val, infile_tab_t *infile_stat)
{
    // combine *.lco to comblco
    int i = 0;
    char *comb_ab_fn, *comb_sketch_fn;
    FILE *comb_sketch_fp, *comb_ab_fp;
    comb_sketch_fn = format_string("%s/%s", sketch_opt_val->outdir, combined_sketch_suffix);
    sprintf(tmp_fname, "%s/%d.%s", sketch_opt_val->outdir, i, sketch_suffix);
    if (rename(tmp_fname, comb_sketch_fn))
        err(errno, "%s():%s rename error", __func__, tmp_fname);
    if ((comb_sketch_fp = fopen(comb_sketch_fn, "ab")) == NULL)
        err(errno, "%s() open file error: %s", __func__, comb_sketch_fn);
    if (sketch_opt_val->abundance)
    {
        comb_ab_fn = format_string("%s.a", comb_sketch_fn);
        sprintf(tmp_fname, "%s/%d.%s.a", sketch_opt_val->outdir, i, sketch_suffix);
        if (rename(tmp_fname, comb_ab_fn))
            err(errno, "%s():%s rename error", __func__, tmp_fname);
        if ((comb_ab_fp = fopen(comb_ab_fn, "ab")) == NULL)
            err(errno, "%s() open file error: %s", __func__, comb_ab_fn);
    }

    for (int i = 1; i < infile_stat->infile_num; i++)
    {
        sprintf(tmp_fname, "%s/%d.%s", sketch_opt_val->outdir, i, sketch_suffix);
        uint64_t *mem_lco = read_from_file(tmp_fname, &file_size);
        fwrite(mem_lco, file_size, 1, comb_sketch_fp);
        remove(tmp_fname);
        free(mem_lco);
        // abundance
        if (!sketch_opt_val->abundance)
            continue;
        sprintf(tmp_fname, "%s/%d.%s.a", sketch_opt_val->outdir, i, sketch_suffix);
        uint32_t *mem_ab = read_from_file(tmp_fname, &file_size);
        fwrite(mem_ab, file_size, 1, comb_ab_fp);
        remove(tmp_fname);
        free(mem_ab);
    }
    fclose(comb_sketch_fp);
    if (sketch_opt_val->abundance)
        fclose(comb_ab_fp);
}

void gen_inverted_index4comblco(const char *refdir)
{

    unify_sketch_t *ref_result = generic_sketch_parse(refdir, SKETCH_PARSE_NONE);
    const_comask_init(&ref_result->stats.lco_stat_val);

    uint64_t sketch_size = ref_result->sketch_index[ref_result->infile_num];
    if (sketch_size >= UINT32_MAX)
        err(EXIT_FAILURE, "%s():sketch_index maximun %lu exceed UINT32_MAX %u", __func__, sketch_size, UINT32_MAX);
    if (ref_result->infile_num >= (1 << GID_NBITS))
        err(EXIT_FAILURE, "%s(): genome numer %d exceed maximum:%u", __func__, ref_result->infile_num, 1 << GID_NBITS);
    if (GID_NBITS + 4 * hclen > 64)
        err(EXIT_FAILURE, "%s(): context_bits_len(%d)+gid_bits_len(%d) exceed 64", __func__, 4 * hclen, GID_NBITS);
    ctxgidobj_t *ctxgidobj = ctxobj64_2ctxgidobj(ref_result->sketch_index, ref_result->comb_sketch, ref_result->infile_num, sketch_size);
    free_unify_sketch(ref_result);
    ctxgidobj_sort_array(ctxgidobj, sketch_size);
    // printf("sketch_size=%lu\t%d\t%d\n",sizeof(ctxgidobj[0]),sizeof(ctxgidobj_t),sketch_size);
    write_to_file(format_string("%s/%s", refdir, sorted_comb_ctxgid64obj32), ctxgidobj, sizeof(ctxgidobj[0]) * sketch_size);
    free(ctxgidobj);
}

static void write_sketch_input_annotations(const char *outdir, infile_tab_t *infile_stat);
static void remove_sketch_annotations(const char *outdir);
static void write_sketch_input_infile_meta(const char *outdir, infile_tab_t *infile_stat);
static void write_sketch_asone_infile_meta(const char *outdir, infile_tab_t *infile_stat);
static void write_sketch_infile_meta_stats(const char *outdir, const infile_meta_t *stats,
                                   size_t sample_count);

static void write_sketch_stat_ex(const char *outdir, infile_tab_t *infile_stat,
                                 bool write_annotations, bool write_infile_meta)
{
    char (*tmpname)[PATHLEN] = malloc(infile_stat->infile_num * PATHLEN);

    for (int i = 0; i < infile_stat->infile_num; i++)
            snprintf(tmpname[i], PATHLEN, "%s", infile_stat->organized_infile_tab[i].fpath);
//        memcpy(tmpname[i], (infile_stat->organized_infile_tab)[i].fpath, PATHLEN);
    concat_and_write_to_file(test_create_fullpath(outdir, sketch_stat), &comblco_stat_one, sizeof(comblco_stat_one), tmpname, infile_stat->infile_num * PATHLEN);
    free(tmpname);
    if (write_annotations)
        write_sketch_input_annotations(outdir, infile_stat);
    else
        remove_sketch_annotations(outdir);
    if (write_infile_meta)
        write_sketch_input_infile_meta(outdir, infile_stat);
}

void write_sketch_stat(const char *outdir, infile_tab_t *infile_stat, bool write_annotations)
{
    write_sketch_stat_ex(outdir, infile_stat, write_annotations, true);
}

static void sketch_annotation_copy(char dest[PATHLEN], const char *name, const char *comment)
{
    memset(dest, 0, PATHLEN);
    size_t pos = 0;
    const char *parts[2] = {name, comment};
    for (int part = 0; part < 2; ++part) {
        const char *src = parts[part];
        if (!src || !src[0])
            continue;
        if (part == 1 && pos > 0 && pos < PATHLEN - 1)
            dest[pos++] = ' ';
        for (; *src && pos < PATHLEN - 1; ++src) {
            const unsigned char c = (unsigned char)*src;
            dest[pos++] = (c < 32 || c == 127 || c == '\t') ? ' ' : (char)c;
        }
    }
    dest[PATHLEN - 1] = '\0';
}

static void write_sketch_annotations(const char *outdir, char (*annotations)[PATHLEN],
                                     size_t sample_count)
{
    if (!annotations || sample_count == 0)
        return;
    write_to_file(test_create_fullpath(outdir, sketch_anno_stat), annotations,
                  sample_count * PATHLEN);
}

static void remove_sketch_annotations(const char *outdir)
{
    char *anno_path = format_string("%s/%s", outdir, sketch_anno_stat);
    if (!anno_path)
        err(errno, "%s(): OOM annotation path", __func__);
    if (unlink(anno_path) != 0 && errno != ENOENT)
        err(errno, "%s(): cannot remove stale %s", __func__, anno_path);
    free(anno_path);
}

static void remove_sketch_positions(const char *outdir)
{
    char *pos_path = format_string("%s/%s", outdir, sketch_position_suffix);
    if (!pos_path)
        err(errno, "%s(): OOM position path", __func__);
    if (unlink(pos_path) != 0 && errno != ENOENT)
        err(errno, "%s(): cannot remove stale %s", __func__, pos_path);
    free(pos_path);
}

static void write_sketch_qc_stats(const char *outdir, const dim_sketch_qc_stat_t *stats,
                                  size_t sample_count)
{
    if (!stats || sample_count == 0)
        return;
    write_to_file(test_create_fullpath(outdir, sketch_qc_stat), stats,
                  sample_count * sizeof(stats[0]));
}

static void append_sketch_qc_stat(dim_sketch_qc_stat_t **stats, size_t *len, size_t *cap,
                                  dim_sketch_qc_stat_t value)
{
    if (!stats || !len || !cap)
        return;
    if (*len == *cap) {
        size_t new_cap = *cap ? (*cap * 2) : 1024;
        dim_sketch_qc_stat_t *new_stats =
            (dim_sketch_qc_stat_t *)realloc(*stats, new_cap * sizeof((*stats)[0]));
        if (!new_stats)
            err(errno, "%s(): OOM sample QC stats", __func__);
        *stats = new_stats;
        *cap = new_cap;
    }
    (*stats)[(*len)++] = value;
}

KSEQ_INIT(gzFile, gzread)

static int cmp_u64_asc_local(const void *a, const void *b)
{
    const uint64_t va = *(const uint64_t *)a;
    const uint64_t vb = *(const uint64_t *)b;
    return (va > vb) - (va < vb);
}

static int8_t infile_fmt_from_path(const char *path)
{
    if (!path)
        return KSSD_INFILE_FMT_UNKNOWN;
    if (isOK_fmt_infile((char *)path, fastq_fmt, FQ_FMT_SZ))
        return KSSD_INFILE_FMT_FASTQ;
    if (isOK_fmt_infile((char *)path, fasta_fmt, FAS_FMT_SZ))
        return KSSD_INFILE_FMT_FASTA;
    return KSSD_INFILE_FMT_UNKNOWN;
}

static uint8_t infile_flags_from_path(const char *path, const sketch_opt_t *opt)
{
    uint8_t flags = 0;
    if (opt && opt->pipecmd && opt->pipecmd[0] != '\0')
        flags |= KSSD_INFILE_FLAG_PIPECMD;
    if (path && strcmp(path, "-") == 0)
        flags |= KSSD_INFILE_FLAG_STDIN;
    if (path && isCompressfile((char *)path))
        flags |= KSSD_INFILE_FLAG_COMPRESSED;
    return flags;
}

static int8_t create_type_from_opt(const sketch_opt_t *opt)
{
    if (opt && opt->pipecmd && opt->pipecmd[0] != '\0')
        return KSSD_CREATE_PIPECMD;
    if (opt && opt->asone)
        return KSSD_CREATE_ASONE;
    if (opt && opt->split_mfa)
        return KSSD_CREATE_SPLITMFA;
    return KSSD_CREATE_NORMAL;
}

static inline bool should_collect_lengths_for_meta(const sketch_opt_t *opt, const char *path)
{
    return opt && opt->compute_meta && infile_fmt_from_path(path) == KSSD_INFILE_FMT_FASTA;
}

static void infile_meta_from_lengths(uint64_t *lengths, size_t n,
                                  int8_t infile_fmt, int8_t create_type,
                                  uint8_t infile_flags, infile_meta_t *out)
{
    memset(out, 0, sizeof(*out));
    out->meta_fmt_version = KSSD_INFILE_META_VERSION;
    out->infile_fmt = infile_fmt;
    out->create_type = create_type;
    out->infile_flags = infile_flags;
    out->asm_level = infile_fmt == KSSD_INFILE_FMT_FASTQ ? 0.0f : -1.0f;
    if (!lengths || n == 0)
        return;

    uint64_t total = 0;
    double mean = 0.0;
    double m2 = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        total += lengths[i];
        const double x = (double)lengths[i];
        const double delta = x - mean;
        mean += delta / (double)(i + 1);
        m2 += delta * (x - mean);
    }
    qsort(lengths, n, sizeof(lengths[0]), cmp_u64_asc_local);

    out->total_length_bp = total;
    out->record_count = n > UINT32_MAX ? UINT32_MAX : (uint32_t)n;
    if (n & 1)
        out->median_length_bp = lengths[n / 2] > UINT32_MAX ? UINT32_MAX : (uint32_t)lengths[n / 2];
    else
    {
        const uint64_t median = (lengths[n / 2 - 1] + lengths[n / 2]) / 2;
        out->median_length_bp = median > UINT32_MAX ? UINT32_MAX : (uint32_t)median;
    }
    if (mean > 0.0 && n > 1)
        out->length_cv = (float)(sqrt(m2 / (double)n) / mean);
    if (infile_fmt == KSSD_INFILE_FMT_FASTA && total > 0)
        out->asm_level = (float)((double)lengths[n - 1] / (double)total);
}

static void append_asm_length(uint64_t **lengths, size_t *n, size_t *cap, uint64_t len)
{
    if (len == 0)
        return;
    if (*n == *cap) {
        size_t new_cap = *cap ? (*cap * 2) : 128;
        uint64_t *new_lengths = (uint64_t *)realloc(*lengths, new_cap * sizeof((*lengths)[0]));
        if (!new_lengths)
            err(errno, "%s(): OOM assembly lengths", __func__);
        *lengths = new_lengths;
        *cap = new_cap;
    }
    (*lengths)[(*n)++] = len;
}

static void append_sketch_infile_meta_stat(infile_meta_t **stats, size_t *len, size_t *cap,
                                   infile_meta_t value)
{
    if (!stats || !len || !cap)
        return;
    if (*len == *cap) {
        size_t new_cap = *cap ? (*cap * 2) : 1024;
        infile_meta_t *new_stats =
            (infile_meta_t *)realloc(*stats, new_cap * sizeof((*stats)[0]));
        if (!new_stats)
            err(errno, "%s(): OOM assembly stats", __func__);
        *stats = new_stats;
        *cap = new_cap;
    }
    (*stats)[(*len)++] = value;
}

static void append_asm_lengths_for_path(const char *path, uint64_t **lengths,
                                        size_t *n, size_t *cap)
{
    if (!path || strcmp(path, "-") == 0)
        return;

    gzFile infile = gzopen(path, "r");
    if (!infile)
        return;
    kseq_t *seq = kseq_init(infile);
    if (!seq) {
        gzclose(infile);
        return;
    }

    while (kseq_read(seq) >= 0)
        append_asm_length(lengths, n, cap, (uint64_t)seq->seq.l);

    kseq_destroy(seq);
    gzclose(infile);
}

static infile_meta_t sketch_infile_meta_stat_for_path(const char *path)
{
    infile_meta_t stat = {0};
    uint64_t *lengths = NULL;
    size_t n = 0, cap = 0;
    if (infile_fmt_from_path(path) == KSSD_INFILE_FMT_FASTA)
        append_asm_lengths_for_path(path, &lengths, &n, &cap);
    infile_meta_from_lengths(lengths, n, infile_fmt_from_path(path),
                             KSSD_CREATE_NORMAL, infile_flags_from_path(path, NULL), &stat);
    free(lengths);
    return stat;
}

static void write_sketch_infile_meta_stats(const char *outdir, const infile_meta_t *stats,
                                   size_t sample_count)
{
    if (!stats || sample_count == 0)
        return;
    write_to_file(test_create_fullpath(outdir, sketch_infile_meta_stat), stats,
                  sample_count * sizeof(stats[0]));
}

static void write_sketch_input_infile_meta(const char *outdir, infile_tab_t *infile_stat)
{
    if (!infile_stat || infile_stat->infile_num <= 0)
        return;
    infile_meta_t *stats =
        (infile_meta_t *)calloc((size_t)infile_stat->infile_num, sizeof(stats[0]));
    if (!stats)
        err(errno, "%s(): OOM assembly stats", __func__);
    for (int i = 0; i < infile_stat->infile_num; ++i)
        stats[i] = sketch_infile_meta_stat_for_path(infile_stat->organized_infile_tab[i].fpath);
    write_sketch_infile_meta_stats(outdir, stats, (size_t)infile_stat->infile_num);
    free(stats);
}

static void write_sketch_asone_infile_meta(const char *outdir, infile_tab_t *infile_stat)
{
    if (!infile_stat || infile_stat->infile_num <= 0)
        return;
    uint64_t *lengths = NULL;
    size_t n = 0, cap = 0;
    bool all_fasta = true;
    for (int i = 0; i < infile_stat->infile_num; ++i) {
        const char *path = infile_stat->organized_infile_tab[i].fpath;
        if (infile_fmt_from_path(path) != KSSD_INFILE_FMT_FASTA) {
            all_fasta = false;
            break;
        }
    }
    if (all_fasta) {
        for (int i = 0; i < infile_stat->infile_num; ++i)
            append_asm_lengths_for_path(infile_stat->organized_infile_tab[i].fpath,
                                        &lengths, &n, &cap);
    }
    infile_meta_t stat = {0};
    infile_meta_from_lengths(lengths, n,
                             all_fasta ? KSSD_INFILE_FMT_FASTA : KSSD_INFILE_FMT_UNKNOWN,
                             KSSD_CREATE_ASONE,
                             all_fasta ? 0 : KSSD_INFILE_FLAG_MIXED_FORMAT,
                             &stat);
    write_sketch_infile_meta_stats(outdir, &stat, 1);
    free(lengths);
}

typedef struct
{
    gzFile gz;
    FILE *pipe_fp;
    char *cmd;
    bool is_pipe;
} sketch_stream_t;

static bool is_stdin_path(const char *path)
{
    return path && strcmp(path, "-") == 0;
}

static bool has_pipecmd(const char *pipecmd)
{
    return pipecmd && pipecmd[0] != '\0';
}

static bool uses_stream_input(const sketch_opt_t *opt, const char *path)
{
    return is_stdin_path(path) || has_pipecmd(opt ? opt->pipecmd : NULL);
}

static char *shell_quote_arg(const char *arg)
{
    size_t need = 3;
    for (const char *p = arg; *p; ++p)
        need += (*p == '\'') ? 4 : 1;

    char *out = malloc(need);
    if (!out)
        err(errno, "%s(): OOM shell quote", __func__);
    char *w = out;
    *w++ = '\'';
    for (const char *p = arg; *p; ++p)
    {
        if (*p == '\'')
        {
            memcpy(w, "'\\''", 4);
            w += 4;
        }
        else
        {
            *w++ = *p;
        }
    }
    *w++ = '\'';
    *w = '\0';
    return out;
}

static char *pipe_command_for_path(const char *pipecmd, const char *path)
{
    char *quoted = shell_quote_arg(path);
    const char *placeholder = strstr(pipecmd, "{}");
    if (!placeholder)
    {
        char *cmd = format_string("%s %s", pipecmd, quoted);
        free(quoted);
        return cmd;
    }

    size_t pipecmd_len = strlen(pipecmd);
    size_t quoted_len = strlen(quoted);
    size_t placeholders = 0;
    for (const char *p = pipecmd; (p = strstr(p, "{}")) != NULL; p += 2)
        placeholders++;

    size_t out_len = pipecmd_len - placeholders * 2 + placeholders * quoted_len;
    char *cmd = malloc(out_len + 1);
    if (!cmd)
        err(errno, "%s(): OOM pipe command", __func__);

    const char *src = pipecmd;
    char *dst = cmd;
    while ((placeholder = strstr(src, "{}")) != NULL)
    {
        size_t chunk = (size_t)(placeholder - src);
        memcpy(dst, src, chunk);
        dst += chunk;
        memcpy(dst, quoted, quoted_len);
        dst += quoted_len;
        src = placeholder + 2;
    }
    strcpy(dst, src);
    free(quoted);
    return cmd;
}

static sketch_stream_t open_sketch_stream(const char *path, const char *pipecmd)
{
    sketch_stream_t stream = {0};
    if (has_pipecmd(pipecmd))
    {
        stream.cmd = pipe_command_for_path(pipecmd, path);
        stream.pipe_fp = popen(stream.cmd, "r");
        if (!stream.pipe_fp)
            err(errno, "%s(): popen %s", __func__, stream.cmd);
        int fd = dup(fileno(stream.pipe_fp));
        if (fd < 0)
            err(errno, "%s(): dup pipe fd", __func__);
        stream.gz = gzdopen(fd, "rb");
        if (!stream.gz)
            err(errno, "%s(): gzdopen pipe %s", __func__, stream.cmd);
        stream.is_pipe = true;
        return stream;
    }

    if (is_stdin_path(path))
    {
        int fd = dup(STDIN_FILENO);
        if (fd < 0)
            err(errno, "%s(): dup stdin", __func__);
        stream.gz = gzdopen(fd, "rb");
        if (!stream.gz)
            err(errno, "%s(): gzdopen stdin", __func__);
        return stream;
    }

    stream.gz = gzopen(path, "r");
    if (!stream.gz)
        err(errno, "%s(): Cannot open %s", __func__, path);
    return stream;
}

static void close_sketch_stream(sketch_stream_t *stream)
{
    if (!stream)
        return;
    if (stream->gz)
    {
        int gz_rc = gzclose(stream->gz);
        stream->gz = NULL;
        if (gz_rc != Z_OK)
            errx(EXIT_FAILURE, "%s(): gzclose failed", __func__);
    }
    if (stream->pipe_fp)
    {
        int pipe_rc = pclose(stream->pipe_fp);
        stream->pipe_fp = NULL;
        if (pipe_rc == -1)
            err(errno, "%s(): pclose %s", __func__, stream->cmd ? stream->cmd : "pipe");
        if (pipe_rc != 0)
            errx(EXIT_FAILURE, "%s(): pipe command failed with status %d: %s",
                 __func__, pipe_rc, stream->cmd ? stream->cmd : "pipe");
    }
    free(stream->cmd);
    stream->cmd = NULL;
}

static bool sketch_has_stdin_input(const infile_tab_t *tab)
{
    if (!tab)
        return false;
    for (int i = 0; i < tab->infile_num; ++i)
        if (is_stdin_path(tab->organized_infile_tab[i].fpath))
            return true;
    return false;
}

static bool sketch_has_stream_input(const sketch_opt_t *opt, const infile_tab_t *tab)
{
    if (has_pipecmd(opt ? opt->pipecmd : NULL))
        return true;
    return sketch_has_stdin_input(tab);
}

static void write_sketch_input_annotations(const char *outdir, infile_tab_t *infile_stat)
{
    if (!infile_stat || infile_stat->infile_num <= 0)
        return;

    char (*annotations)[PATHLEN] = calloc((size_t)infile_stat->infile_num, PATHLEN);
    if (!annotations)
        err(errno, "%s(): OOM annotations", __func__);

    for (int i = 0; i < infile_stat->infile_num; ++i) {
        const char *path = infile_stat->organized_infile_tab[i].fpath;
        gzFile infile = gzopen(path, "r");
        if (!infile)
            err(errno, "%s(): Cannot open file %s", __func__, path);
        kseq_t *seq = kseq_init(infile);
        if (kseq_read(seq) >= 0) {
            sketch_annotation_copy(annotations[i], seq->name.s,
                                   seq->comment.l ? seq->comment.s : NULL);
        }
        kseq_destroy(seq);
        gzclose(infile);
    }

    write_sketch_annotations(outdir, annotations, (size_t)infile_stat->infile_num);
    free(annotations);
}

// keep k-mer using hashtable is faster(~16s) than using vector and sort(~24) per 1k genomes, probabaly due to cache competition
int seq2ht_sortedctxobj64(char *seqfname, char *outfname, bool abundance, int n)
{
    int TL = klen;
    gzFile infile = gzopen(seqfname, "r");
    if (!infile)
        err(errno, "reads2sketch64(): Cannot open file %s", seqfname);
    kseq_t *seq = kseq_init(infile);
    khash_t(sort64) *h = kh_init(sort64);
    uint64_t tuple, crvstuple, unituple, basenum, unictx;
    uint32_t len_mv = 2 * TL - 2;
    uint32_t sketch_size = 0;

    while (kseq_read(seq) >= 0)
    {
        const char *s = seq->seq.s;
        if (seq->seq.l < TL)
            continue;
        int base = 0;

        for (int pos = 0; pos < seq->seq.l; pos++)
        { // for(int pos = TL; pos < seq->seq.l ; pos++)
            if (Basemap[(unsigned short)s[pos]] == DEFAULT)
            {
                base = 0;
                continue;
            }
            basenum = Basemap[(unsigned short)s[pos]];
            tuple = ((tuple << 2) | basenum);
            crvstuple = ((crvstuple >> 2) | ((basenum ^ 3LLU) << len_mv));
            if (++base < TL)
                continue;
            // if base >=TL, namely, contiue ACGT TL-mer
            unituple = (tuple & ctxmask) < (crvstuple & ctxmask) ? tuple : crvstuple;
            unictx = unituple & ctxmask;
            if (SKETCH_HASH(unictx) > FILTER)
                continue;
            int ret;
            khint_t key = kh_put(sort64, h, uint64kmer2generic_ctxobj(unituple), &ret);

            if (ret)
            {
                kh_value(h, key) = 1;
                sketch_size++;
            }
            else
            {
                kh_value(h, key)++;
            }
        } // for line
    }; // while
    kseq_destroy(seq);
    gzclose(infile);
    // get sketch sorted
    SortedKV_Arrays_t lco_ab = sort_khash_u64(h);
    //	for(int i = 0; i<lco_ab.len;i++){		printf("%d\t%d\t%lx\n",i,lco_ab.len,lco_ab.keys[i]);}
    kh_destroy(sort64, h);
    if (n > 1)
        filter_n_SortedKV_Arrays(&lco_ab, n);
    write_to_file(outfname, lco_ab.keys, lco_ab.len * sizeof(lco_ab.keys[0]));
    if (abundance)
        write_to_file(format_string("%s.a", outfname), lco_ab.values, lco_ab.len * sizeof(lco_ab.values[0]));
    free_all(lco_ab.keys, lco_ab.values, NULL);
    return lco_ab.len; // sketch_size;
}
//
int seq2sortedsketch64(char *seqfname, char *outfname, bool abundance, int n)
{

    gzFile infile = gzopen(seqfname, "r");
    if (!infile)
        err(errno, "reads2sketch64(): Cannot open file %s", seqfname);
    kseq_t *seq = kseq_init(infile);
    Vector raw_sketch;
    vector_init(&raw_sketch, sizeof(uint64_t));

    uint64_t tuple, crvstuple, unituple, basenum, unictx;
    uint32_t len_mv = 2 * klen - 2;

    while (kseq_read(seq) >= 0)
    {
        const char *s = seq->seq.s;

        if (seq->seq.l < klen)
            continue;
        int base = 0;

        for (int pos = 0; pos < seq->seq.l; pos++)
        { // for(int pos = klen; pos < seq->seq.l ; pos++)
            if (Basemap[(unsigned short)s[pos]] == DEFAULT)
            {
                base = 0;
                continue;
            }
            basenum = Basemap[(unsigned short)s[pos]];
            tuple = ((tuple << 2) | basenum);
            crvstuple = ((crvstuple >> 2) | ((basenum ^ 3LLU) << len_mv));
            if (++base < klen)
                continue;
            // if base >=klen, namely, contiue ACGT klen-mer
            unituple = (tuple & ctxmask) < (crvstuple & ctxmask) ? tuple : crvstuple;
            unictx = unituple & ctxmask;
            if (SKETCH_HASH(unictx) > FILTER)
                continue;

            unituple = uint64kmer2generic_ctxobj(unituple);
            vector_push(&raw_sketch, &unituple);

        } // for line
    }; // while
    gzclose(infile);

    // write sketch and abundance
    uint64_t *mem_lco = raw_sketch.data;
    qsort(mem_lco, raw_sketch.size, sizeof(uint64_t), qsort_comparator_uint64);
    uint32_t *mem_ab; // if (abundance) mem_ab = malloc(sketch_size * sizeof(uint32_t));
    uint32_t sketch_size, kmer_ct;
    sketch_size = kmer_ct = dedup_with_counts(mem_lco, raw_sketch.size, &mem_ab);

    if (n > 1)
    {
        kmer_ct = 0;
        for (size_t i = 0; i < sketch_size; i++)
        {
            if (mem_ab[i] >= n)
            {
                mem_lco[kmer_ct] = mem_lco[i]; // Shift unique element forward
                mem_ab[kmer_ct] = mem_ab[i];
                kmer_ct++;
            }
        }
    }

    write_to_file(outfname, mem_lco, kmer_ct * sizeof(uint64_t));
    vector_free(&raw_sketch);

    if (abundance)
    {
        write_to_file(format_string("%s.a", outfname), mem_ab, kmer_ct * sizeof(uint32_t));
        free(mem_ab);
    }
    return kmer_ct; // sketch_size;
}


int merge_comblco(sketch_opt_t *sketch_opt_val)
{

    size_t first_stat_size = 0;
    void *mem_stat = read_from_file(test_get_fullpath(sketch_opt_val->remaining_args[0], sketch_stat), &first_stat_size);
    memcpy(&comblco_stat_one, mem_stat, sizeof(comblco_stat_one));
    comblco_stat_one.infile_num = 0;
    char (*tmpname)[PATHLEN] = malloc(PATHLEN);
    FILE *merge_out_fp = fopen(test_create_fullpath(sketch_opt_val->outdir, combined_sketch_suffix), "wb");
    if (merge_out_fp == NULL)
        err(errno, "%s():%s/%s", __func__, sketch_opt_val->outdir, combined_sketch_suffix);
    bool merge_all_have_positions = true;
    for (int i = 0; i < sketch_opt_val->num_remaining_args; ++i) {
        if (!file_exists_in_folder(sketch_opt_val->remaining_args[i], sketch_position_suffix)) {
            merge_all_have_positions = false;
            break;
        }
    }
    FILE *merge_pos_fp = NULL;
    if (merge_all_have_positions) {
        merge_pos_fp = fopen(test_create_fullpath(sketch_opt_val->outdir, sketch_position_suffix), "wb");
        if (merge_pos_fp == NULL)
            err(errno, "%s():%s/%s", __func__, sketch_opt_val->outdir, sketch_position_suffix);
    } else {
        remove_sketch_positions(sketch_opt_val->outdir);
    }

    uint64_t *index_arry = malloc(sizeof(uint64_t) * 100);
    index_arry[0] = 0;
    dim_sketch_qc_stat_t *merged_qc_stats = NULL;
    bool write_merged_qc_stats = false;
    infile_meta_t *merged_infile_meta = NULL;
    bool write_merged_infile_meta = false;
    bool merge_seen_meta = false;
    bool merge_seen_no_meta = false;
    char (*merged_annotations)[PATHLEN] = NULL;
    bool write_merged_annotations = false;
    for (int i = 0; i < sketch_opt_val->num_remaining_args; i++)
    {
        // read stat
        size_t stat_it_size = 0;
        void *mem_stat_it = read_from_file(test_get_fullpath(sketch_opt_val->remaining_args[i], sketch_stat), &stat_it_size);
        memcpy(&comblco_stat_it, mem_stat_it, sizeof(comblco_stat_it));
        if (comblco_stat_it.hash_id != comblco_stat_one.hash_id)
            err(EINVAL, "%uth %s hashid: %u != %u ", i, sketch_opt_val->remaining_args[i], comblco_stat_it.hash_id, comblco_stat_one.hash_id);
        if (comblco_stat_it.koc == 0)
            comblco_stat_one.koc = 0;
        const int old_infile_num = comblco_stat_one.infile_num;
        const int new_infile_num = old_infile_num + comblco_stat_it.infile_num;
        // collect filenames
        tmpname = realloc(tmpname, PATHLEN * new_infile_num);
        memcpy(tmpname + old_infile_num, mem_stat_it + sizeof(comblco_stat_it), PATHLEN * comblco_stat_it.infile_num);
        // set index
        size_t index_it_size = 0;
        uint64_t *mem_index_it = read_from_file(test_get_fullpath(sketch_opt_val->remaining_args[i], idx_sketch_suffix), &index_it_size);
        index_arry = (uint64_t *)realloc(index_arry, sizeof(uint64_t) * (new_infile_num + 1));
        for (int j = 1; j < comblco_stat_it.infile_num + 1; j++)
            index_arry[old_infile_num + j] = index_arry[old_infile_num] + mem_index_it[j];

        const bool input_has_qc_stats =
            file_exists_in_folder(sketch_opt_val->remaining_args[i], sketch_qc_stat);
        if (input_has_qc_stats || write_merged_qc_stats) {
            dim_sketch_qc_stat_t *new_qc_stats =
                (dim_sketch_qc_stat_t *)realloc(merged_qc_stats,
                                                sizeof(merged_qc_stats[0]) * (size_t)new_infile_num);
            if (!new_qc_stats)
                err(errno, "%s(): OOM merged sample QC stats", __func__);
            merged_qc_stats = new_qc_stats;

            if (!write_merged_qc_stats && old_infile_num > 0)
                memset(merged_qc_stats, 0, sizeof(merged_qc_stats[0]) * (size_t)old_infile_num);

            if (input_has_qc_stats) {
                size_t qc_file_size = 0;
                dim_sketch_qc_stat_t *mem_qc_it =
                    read_from_file(test_get_fullpath(sketch_opt_val->remaining_args[i], sketch_qc_stat),
                                   &qc_file_size);
                const size_t expected_size =
                    sizeof(mem_qc_it[0]) * (size_t)comblco_stat_it.infile_num;
                if (qc_file_size != expected_size)
                    err(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
                        __func__, sketch_opt_val->remaining_args[i], sketch_qc_stat,
                        qc_file_size, expected_size);
                memcpy(merged_qc_stats + old_infile_num, mem_qc_it, expected_size);
                free_read_from_file(mem_qc_it, qc_file_size);
            } else {
                memset(merged_qc_stats + old_infile_num, 0,
                       sizeof(merged_qc_stats[0]) * (size_t)comblco_stat_it.infile_num);
            }
            write_merged_qc_stats = true;
        }

        const bool input_has_infile_meta =
            file_exists_in_folder(sketch_opt_val->remaining_args[i], sketch_infile_meta_stat);
        merge_seen_meta = merge_seen_meta || input_has_infile_meta;
        merge_seen_no_meta = merge_seen_no_meta || !input_has_infile_meta;
        if (input_has_infile_meta || write_merged_infile_meta) {
            infile_meta_t *new_infile_meta =
                (infile_meta_t *)realloc(merged_infile_meta,
                                                 sizeof(merged_infile_meta[0]) * (size_t)new_infile_num);
            if (!new_infile_meta)
                err(errno, "%s(): OOM merged infile metadata", __func__);
            merged_infile_meta = new_infile_meta;

            if (!write_merged_infile_meta && old_infile_num > 0)
                memset(merged_infile_meta, 0, sizeof(merged_infile_meta[0]) * (size_t)old_infile_num);

            if (input_has_infile_meta) {
                size_t meta_file_size = 0;
                infile_meta_t *mem_meta_it =
                    read_from_file(test_get_fullpath(sketch_opt_val->remaining_args[i], sketch_infile_meta_stat),
                                   &meta_file_size);
                const size_t expected_size =
                    sizeof(mem_meta_it[0]) * (size_t)comblco_stat_it.infile_num;
                if (meta_file_size != expected_size)
                    err(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
                        __func__, sketch_opt_val->remaining_args[i], sketch_infile_meta_stat,
                        meta_file_size, expected_size);
                memcpy(merged_infile_meta + old_infile_num, mem_meta_it, expected_size);
                free_read_from_file(mem_meta_it, meta_file_size);
            } else {
                memset(merged_infile_meta + old_infile_num, 0,
                       sizeof(merged_infile_meta[0]) * (size_t)comblco_stat_it.infile_num);
            }
            write_merged_infile_meta = true;
        }

        const bool input_has_annotations =
            file_exists_in_folder(sketch_opt_val->remaining_args[i], sketch_anno_stat);
        if (input_has_annotations || write_merged_annotations) {
            char (*new_annotations)[PATHLEN] =
                (char (*)[PATHLEN])realloc(merged_annotations,
                                           (size_t)new_infile_num * PATHLEN);
            if (!new_annotations)
                err(errno, "%s(): OOM merged annotations", __func__);
            merged_annotations = new_annotations;

            if (!write_merged_annotations && old_infile_num > 0)
                memset(merged_annotations, 0, (size_t)old_infile_num * PATHLEN);

            if (input_has_annotations) {
                size_t anno_file_size = 0;
                char (*mem_anno_it)[PATHLEN] =
                    read_from_file(test_get_fullpath(sketch_opt_val->remaining_args[i],
                                                     sketch_anno_stat),
                                   &anno_file_size);
                const size_t expected_size = (size_t)comblco_stat_it.infile_num * PATHLEN;
                if (anno_file_size != expected_size)
                    err(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
                        __func__, sketch_opt_val->remaining_args[i], sketch_anno_stat,
                        anno_file_size, expected_size);
                memcpy(merged_annotations + old_infile_num, mem_anno_it, expected_size);
                free_read_from_file(mem_anno_it, anno_file_size);
            } else {
                memset(merged_annotations + old_infile_num, 0,
                       (size_t)comblco_stat_it.infile_num * PATHLEN);
            }
            write_merged_annotations = true;
        }

        // add file num
        comblco_stat_one.infile_num = new_infile_num;
        // write combined_sketch_suffix
        size_t comblco_it_size = 0;
        uint64_t *mem_comblco = read_from_file(test_get_fullpath(sketch_opt_val->remaining_args[i], combined_sketch_suffix), &comblco_it_size);
        fwrite(mem_comblco, comblco_it_size, 1, merge_out_fp);
        if (merge_pos_fp) {
            size_t pos_it_size = 0;
            uint64_t *mem_pos = read_from_file(test_get_fullpath(sketch_opt_val->remaining_args[i], sketch_position_suffix), &pos_it_size);
            if (pos_it_size != comblco_it_size)
                err(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
                    __func__, sketch_opt_val->remaining_args[i], sketch_position_suffix,
                    pos_it_size, comblco_it_size);
            fwrite(mem_pos, pos_it_size, 1, merge_pos_fp);
            free_read_from_file(mem_pos, pos_it_size);
        }
        free_read_from_file(mem_stat_it, stat_it_size);
        free_read_from_file(mem_index_it, index_it_size);
        free_read_from_file(mem_comblco, comblco_it_size);
    }
    fclose(merge_out_fp); // write combined_sketch_suffix complete
    if (merge_pos_fp)
        fclose(merge_pos_fp);

    if (comblco_stat_one.koc)
    {
        merge_out_fp = fopen(test_create_fullpath(sketch_opt_val->outdir, combined_ab_suffix), "wb");
        if (merge_out_fp == NULL)
            err(errno, "%s():%s/%s", __func__, sketch_opt_val->outdir, combined_ab_suffix);
        for (int i = 0; i < sketch_opt_val->num_remaining_args; i++)
        {
            size_t ab_it_size = 0;
            uint64_t *mem_ab = read_from_file(test_get_fullpath(sketch_opt_val->remaining_args[i], combined_ab_suffix), &ab_it_size);
            fwrite(mem_ab, ab_it_size, 1, merge_out_fp);
            free_read_from_file(mem_ab, ab_it_size);
        }
        fclose(merge_out_fp);
    }
    // write index and stat file
    write_to_file(test_create_fullpath(sketch_opt_val->outdir, idx_sketch_suffix), index_arry, sizeof(index_arry[0]) * (comblco_stat_one.infile_num + 1));
    concat_and_write_to_file(test_create_fullpath(sketch_opt_val->outdir, sketch_stat), &comblco_stat_one, sizeof(comblco_stat_one), tmpname, PATHLEN * (comblco_stat_one.infile_num));
    if (write_merged_qc_stats)
        write_sketch_qc_stats(sketch_opt_val->outdir, merged_qc_stats, (size_t)comblco_stat_one.infile_num);
    if (write_merged_infile_meta) {
        if (merge_seen_meta && merge_seen_no_meta)
            memset(merged_infile_meta, 0,
                   sizeof(merged_infile_meta[0]) * (size_t)comblco_stat_one.infile_num);
        write_sketch_infile_meta_stats(sketch_opt_val->outdir, merged_infile_meta, (size_t)comblco_stat_one.infile_num);
    }
    if (write_merged_annotations)
        write_sketch_annotations(sketch_opt_val->outdir, merged_annotations, (size_t)comblco_stat_one.infile_num);
    free_read_from_file(mem_stat, first_stat_size);
    free(index_arry);
    free(tmpname);
    free(merged_qc_stats);
    free(merged_infile_meta);
    free(merged_annotations);

    return comblco_stat_one.infile_num;
}

typedef struct append_sketch_part
{
    const char *path;
    dim_sketch_stat_t stat;
    void *stat_mem;
    size_t stat_size;
    uint64_t *index;
    size_t index_size;
    uint64_t entries;
    size_t comblco_size;
    bool has_positions;
    size_t position_size;
    size_t abundance_size;
    bool has_qc;
    dim_sketch_qc_stat_t *qc;
    size_t qc_size;
    bool has_meta;
    infile_meta_t *meta;
    size_t meta_size;
    bool has_anno;
    char (*anno)[PATHLEN];
    size_t anno_size;
} append_sketch_part_t;

typedef struct append_payload_rollback
{
    char *comblco_path;
    char *abundance_path;
    char *position_path;
    size_t comblco_size;
    size_t abundance_size;
    size_t position_size;
    bool has_abundance;
    bool has_positions;
} append_payload_rollback_t;

static char *append_join_path(const char *dir, const char *suffix)
{
    char *path = format_string("%s/%s", dir, suffix);
    if (!path)
        err(EXIT_FAILURE, "%s(): failed to allocate path", __func__);
    return path;
}

static bool append_file_exists(const char *dir, const char *suffix)
{
    char *path = append_join_path(dir, suffix);
    const bool exists = access(path, F_OK) == 0;
    free(path);
    return exists;
}

static size_t append_file_size_path(const char *path)
{
    struct stat st;
    if (stat(path, &st) != 0)
        err(errno, "%s(): cannot stat %s", __func__, path);
    if (!S_ISREG(st.st_mode))
        errx(EINVAL, "%s(): %s is not a regular file", __func__, path);
    if (st.st_size < 0)
        errx(EINVAL, "%s(): %s has invalid size", __func__, path);
    return (size_t)st.st_size;
}

static size_t append_file_size(const char *dir, const char *suffix)
{
    char *path = append_join_path(dir, suffix);
    const size_t size = append_file_size_path(path);
    free(path);
    return size;
}

static size_t append_entries_bytes(uint64_t entries, size_t elem_size, const char *path, const char *suffix)
{
    if (entries > SIZE_MAX / elem_size)
        errx(EINVAL, "%s(): %s/%s entry count is too large", __func__, path, suffix);
    return (size_t)entries * elem_size;
}

static void append_validate_regular_size(const char *dir, const char *suffix, size_t expected)
{
    const size_t actual = append_file_size(dir, suffix);
    if (actual != expected)
        errx(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
             __func__, dir, suffix, actual, expected);
}

static void append_load_optional_qc(append_sketch_part_t *part)
{
    part->has_qc = append_file_exists(part->path, sketch_qc_stat);
    if (!part->has_qc)
        return;
    char *path = append_join_path(part->path, sketch_qc_stat);
    part->qc = read_from_file(path, &part->qc_size);
    free(path);
    const size_t expected = (size_t)part->stat.infile_num * sizeof(part->qc[0]);
    if (part->qc_size != expected)
        errx(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
             __func__, part->path, sketch_qc_stat, part->qc_size, expected);
}

static void append_load_optional_meta(append_sketch_part_t *part)
{
    part->has_meta = append_file_exists(part->path, sketch_infile_meta_stat);
    if (!part->has_meta)
        return;
    char *path = append_join_path(part->path, sketch_infile_meta_stat);
    part->meta = read_from_file(path, &part->meta_size);
    free(path);
    const size_t expected = (size_t)part->stat.infile_num * sizeof(part->meta[0]);
    if (part->meta_size != expected)
        errx(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
             __func__, part->path, sketch_infile_meta_stat, part->meta_size, expected);
}

static void append_load_optional_anno(append_sketch_part_t *part)
{
    part->has_anno = append_file_exists(part->path, sketch_anno_stat);
    if (!part->has_anno)
        return;
    char *path = append_join_path(part->path, sketch_anno_stat);
    part->anno = read_from_file(path, &part->anno_size);
    free(path);
    const size_t expected = (size_t)part->stat.infile_num * PATHLEN;
    if (part->anno_size != expected)
        errx(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
             __func__, part->path, sketch_anno_stat, part->anno_size, expected);
}

static void append_load_part(append_sketch_part_t *part, const char *path)
{
    memset(part, 0, sizeof(*part));
    part->path = path;
    char *stat_path = append_join_path(path, sketch_stat);
    part->stat_mem = read_from_file(stat_path, &part->stat_size);
    free(stat_path);
    if (part->stat_size < sizeof(part->stat))
        errx(EINVAL, "%s(): %s/%s has %zu bytes, expected at least %zu",
             __func__, path, sketch_stat, part->stat_size, sizeof(part->stat));
    memcpy(&part->stat, part->stat_mem, sizeof(part->stat));
    if (part->stat.infile_num < 0)
        errx(EINVAL, "%s(): %s has negative sample count", __func__, path);
    const size_t expected_stat =
        sizeof(part->stat) + (size_t)part->stat.infile_num * PATHLEN;
    if (part->stat_size != expected_stat)
        errx(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
             __func__, path, sketch_stat, part->stat_size, expected_stat);

    char *index_path = append_join_path(path, idx_sketch_suffix);
    part->index = read_from_file(index_path, &part->index_size);
    free(index_path);
    const size_t expected_index =
        ((size_t)part->stat.infile_num + 1) * sizeof(part->index[0]);
    if (part->index_size != expected_index)
        errx(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
             __func__, path, idx_sketch_suffix, part->index_size, expected_index);
    if (part->index[0] != 0)
        errx(EINVAL, "%s(): %s/%s does not start at zero", __func__, path, idx_sketch_suffix);
    part->entries = part->index[part->stat.infile_num];
    for (int i = 0; i < part->stat.infile_num; ++i) {
        if (part->index[i + 1] < part->index[i])
            errx(EINVAL, "%s(): %s/%s is not monotonic at sample %d",
                 __func__, path, idx_sketch_suffix, i);
    }

    part->comblco_size =
        append_entries_bytes(part->entries, sizeof(uint64_t), path, combined_sketch_suffix);
    append_validate_regular_size(path, combined_sketch_suffix, part->comblco_size);
    if (part->stat.koc) {
        part->abundance_size =
            append_entries_bytes(part->entries, sizeof(uint32_t), path, combined_ab_suffix);
        append_validate_regular_size(path, combined_ab_suffix, part->abundance_size);
    }
    part->has_positions = append_file_exists(path, sketch_position_suffix);
    if (part->has_positions) {
        part->position_size =
            append_entries_bytes(part->entries, sizeof(uint64_t), path, sketch_position_suffix);
        append_validate_regular_size(path, sketch_position_suffix, part->position_size);
    }

    append_load_optional_qc(part);
    append_load_optional_meta(part);
    append_load_optional_anno(part);
}

static void append_free_part(append_sketch_part_t *part)
{
    if (part->stat_mem)
        free_read_from_file(part->stat_mem, part->stat_size);
    if (part->index)
        free_read_from_file(part->index, part->index_size);
    if (part->qc)
        free_read_from_file(part->qc, part->qc_size);
    if (part->meta)
        free_read_from_file(part->meta, part->meta_size);
    if (part->anno)
        free_read_from_file(part->anno, part->anno_size);
}

static bool append_same_sketch_dir(const char *a, const char *b)
{
    struct stat sa, sb;
    if (stat(a, &sa) != 0)
        err(errno, "%s(): cannot stat %s", __func__, a);
    if (stat(b, &sb) != 0)
        err(errno, "%s(): cannot stat %s", __func__, b);
    if (!S_ISDIR(sa.st_mode))
        errx(EINVAL, "%s(): %s is not a directory", __func__, a);
    if (!S_ISDIR(sb.st_mode))
        errx(EINVAL, "%s(): %s is not a directory", __func__, b);
    return sa.st_dev == sb.st_dev && sa.st_ino == sb.st_ino;
}

static bool append_same_existing_sketch_dir(const char *maybe_existing, const char *required)
{
    struct stat sa, sb;
    if (stat(maybe_existing, &sa) != 0) {
        if (errno == ENOENT)
            return false;
        err(errno, "%s(): cannot stat %s", __func__, maybe_existing);
    }
    if (stat(required, &sb) != 0)
        err(errno, "%s(): cannot stat %s", __func__, required);
    if (!S_ISDIR(sa.st_mode))
        errx(EINVAL, "%s(): %s exists but is not a directory", __func__, maybe_existing);
    if (!S_ISDIR(sb.st_mode))
        errx(EINVAL, "%s(): %s is not a directory", __func__, required);
    return sa.st_dev == sb.st_dev && sa.st_ino == sb.st_ino;
}

static void append_check_compatible_stat(const dim_sketch_stat_t *target,
                                         const dim_sketch_stat_t *source,
                                         const char *source_path)
{
    if (target->hash_id != source->hash_id || target->koc != source->koc ||
        target->conflict != source->conflict || target->coden_len != source->coden_len ||
        target->klen != source->klen || target->hclen != source->hclen ||
        target->holen != source->holen || target->drfold != source->drfold) {
        errx(EINVAL,
             "%s(): %s is not compatible with target "
             "(hash_id %u/%u, koc %d/%d, conflict %d/%d, coden_len %d/%d, "
             "klen %d/%d, hclen %d/%d, holen %d/%d, drfold %d/%d)",
             __func__, source_path,
             source->hash_id, target->hash_id,
             source->koc, target->koc,
             source->conflict, target->conflict,
             source->coden_len, target->coden_len,
             source->klen, target->klen,
             source->hclen, target->hclen,
             source->holen, target->holen,
             source->drfold, target->drfold);
    }
}

static char (*append_part_names(const append_sketch_part_t *part))[PATHLEN]
{
    return (char (*)[PATHLEN])((char *)part->stat_mem + sizeof(dim_sketch_stat_t));
}

static char *append_tmp_path(const char *target_dir, const char *suffix)
{
    char *path = format_string("%s/%s.append.%ld.tmp", target_dir, suffix, (long)getpid());
    if (!path)
        err(EXIT_FAILURE, "%s(): failed to allocate temporary path", __func__);
    return path;
}

static void append_rollback_payloads(const append_payload_rollback_t *rollback)
{
    if (rollback->comblco_path &&
        truncate(rollback->comblco_path, (off_t)rollback->comblco_size) != 0)
        warn("%s(): failed to truncate %s during rollback", __func__, rollback->comblco_path);
    if (rollback->has_abundance && rollback->abundance_path &&
        truncate(rollback->abundance_path, (off_t)rollback->abundance_size) != 0)
        warn("%s(): failed to truncate %s during rollback", __func__, rollback->abundance_path);
    if (rollback->has_positions && rollback->position_path &&
        truncate(rollback->position_path, (off_t)rollback->position_size) != 0)
        warn("%s(): failed to truncate %s during rollback", __func__, rollback->position_path);
}

static void append_copy_file_or_rollback(const char *src_path, const char *dst_path,
                                         const append_payload_rollback_t *rollback)
{
    FILE *src = fopen(src_path, "rb");
    if (!src) {
        append_rollback_payloads(rollback);
        err(errno, "%s(): cannot open %s", __func__, src_path);
    }
    FILE *dst = fopen(dst_path, "ab");
    if (!dst) {
        fclose(src);
        append_rollback_payloads(rollback);
        err(errno, "%s(): cannot open %s for append", __func__, dst_path);
    }
    const size_t buf_size = 8u << 20;
    char *buf = malloc(buf_size);
    if (!buf) {
        fclose(src);
        fclose(dst);
        append_rollback_payloads(rollback);
        err(EXIT_FAILURE, "%s(): OOM copy buffer", __func__);
    }
    while (true) {
        size_t n = fread(buf, 1, buf_size, src);
        if (n > 0 && fwrite(buf, 1, n, dst) != n) {
            free(buf);
            fclose(src);
            fclose(dst);
            append_rollback_payloads(rollback);
            err(errno, "%s(): failed writing %s", __func__, dst_path);
        }
        if (n < buf_size) {
            if (ferror(src)) {
                free(buf);
                fclose(src);
                fclose(dst);
                append_rollback_payloads(rollback);
                err(errno, "%s(): failed reading %s", __func__, src_path);
            }
            break;
        }
    }
    free(buf);
    if (fclose(src) != 0) {
        fclose(dst);
        append_rollback_payloads(rollback);
        err(errno, "%s(): failed closing %s", __func__, src_path);
    }
    if (fclose(dst) != 0) {
        append_rollback_payloads(rollback);
        err(errno, "%s(): failed closing %s", __func__, dst_path);
    }
}

static void append_replace_tmp_or_rollback(const char *tmp_path, const char *target_dir,
                                           const char *suffix,
                                           const append_payload_rollback_t *rollback)
{
    char *dst_path = append_join_path(target_dir, suffix);
    if (rename(tmp_path, dst_path) != 0) {
        append_rollback_payloads(rollback);
        err(errno, "%s(): failed to replace %s", __func__, dst_path);
    }
    free(dst_path);
}

typedef struct remove_name_entry
{
    char *name;
    bool found;
} remove_name_entry_t;

static int remove_name_entry_cmp(const void *a, const void *b)
{
    const remove_name_entry_t *ea = (const remove_name_entry_t *)a;
    const remove_name_entry_t *eb = (const remove_name_entry_t *)b;
    return strcmp(ea->name, eb->name);
}

static void remove_free_name_entries(remove_name_entry_t *entries, size_t n)
{
    if (!entries)
        return;
    for (size_t i = 0; i < n; ++i)
        free(entries[i].name);
    free(entries);
}

static remove_name_entry_t *remove_load_name_list(const char *path,
                                                  const char *operation,
                                                  size_t *n_out)
{
    FILE *fp = fopen(path, "r");
    if (!fp)
        err(errno, "%s(): cannot open %s list %s", __func__, operation, path);

    remove_name_entry_t *entries = NULL;
    size_t n = 0;
    size_t cap = 0;
    char *line = NULL;
    size_t line_cap = 0;
    ssize_t line_len;
    while ((line_len = getline(&line, &line_cap, fp)) >= 0) {
        while (line_len > 0 && (line[line_len - 1] == '\n' || line[line_len - 1] == '\r'))
            line[--line_len] = '\0';
        if (line_len == 0)
            continue;
        if ((size_t)line_len >= PATHLEN)
            errx(EINVAL, "%s(): %s-list sample name is too long for lcofiles.stat: %s",
                 __func__, operation, line);
        if (n == cap) {
            size_t new_cap = cap ? cap * 2 : 16;
            remove_name_entry_t *new_entries =
                (remove_name_entry_t *)realloc(entries, new_cap * sizeof(entries[0]));
            if (!new_entries)
                err(EXIT_FAILURE, "%s(): OOM %s-name list", __func__, operation);
            entries = new_entries;
            cap = new_cap;
        }
        entries[n].name = strdup(line);
        if (!entries[n].name)
            err(EXIT_FAILURE, "%s(): OOM %s-name entry", __func__, operation);
        entries[n].found = false;
        ++n;
    }
    if (ferror(fp))
        err(errno, "%s(): failed reading %s list %s", __func__, operation, path);
    free(line);
    fclose(fp);

    if (n == 0)
        errx(EINVAL, "%s(): %s list %s contains no sample names", __func__, operation, path);

    qsort(entries, n, sizeof(entries[0]), remove_name_entry_cmp);
    size_t w = 0;
    for (size_t r = 0; r < n; ++r) {
        if (w > 0 && strcmp(entries[w - 1].name, entries[r].name) == 0) {
            free(entries[r].name);
            continue;
        }
        if (w != r)
            entries[w] = entries[r];
        ++w;
    }

    *n_out = w;
    return entries;
}

static bool remove_find_name(remove_name_entry_t *entries, size_t n,
                             const char *name, size_t *idx_out)
{
    size_t lo = 0;
    size_t hi = n;
    while (lo < hi) {
        const size_t mid = lo + (hi - lo) / 2;
        const int cmp = strcmp(name, entries[mid].name);
        if (cmp == 0) {
            if (idx_out)
                *idx_out = mid;
            return true;
        }
        if (cmp < 0)
            hi = mid;
        else
            lo = mid + 1;
    }
    return false;
}

static uint64_t remove_payload_bytes(uint64_t entries, size_t elem_size,
                                     const char *suffix)
{
    if (elem_size == 0 || entries > UINT64_MAX / elem_size)
        errx(EINVAL, "%s(): %s byte count overflows uint64_t", __func__, suffix);
    return entries * (uint64_t)elem_size;
}

static void remove_copy_exact(FILE *src, FILE *dst, uint64_t bytes,
                              char *buf, size_t buf_size,
                              const char *src_path, const char *dst_path)
{
    uint64_t remaining = bytes;
    while (remaining > 0) {
        const size_t chunk = remaining > (uint64_t)buf_size ? buf_size : (size_t)remaining;
        const size_t got = fread(buf, 1, chunk, src);
        if (got != chunk) {
            if (ferror(src))
                err(errno, "%s(): failed reading %s", __func__, src_path);
            errx(EINVAL, "%s(): short read from %s", __func__, src_path);
        }
        if (fwrite(buf, 1, chunk, dst) != chunk)
            err(errno, "%s(): failed writing %s", __func__, dst_path);
        remaining -= chunk;
    }
}

static void remove_skip_exact(FILE *src, uint64_t bytes,
                              const char *src_path, const char *suffix)
{
    while (bytes > 0) {
        const uint64_t step = bytes > (uint64_t)LLONG_MAX ? (uint64_t)LLONG_MAX : bytes;
        if (fseeko(src, (off_t)step, SEEK_CUR) != 0)
            err(errno, "%s(): failed seeking over removed %s bytes in %s",
                __func__, suffix, src_path);
        bytes -= step;
    }
}

static void remove_copy_filtered_payload(const char *target_dir, const char *suffix,
                                         const char *tmp_path,
                                         const uint64_t *index, int nfiles,
                                         const bool *remove_sample,
                                         size_t elem_size)
{
    char *src_path = append_join_path(target_dir, suffix);
    FILE *src = fopen(src_path, "rb");
    if (!src)
        err(errno, "%s(): cannot open %s", __func__, src_path);
    FILE *dst = fopen(tmp_path, "wb");
    if (!dst) {
        fclose(src);
        err(errno, "%s(): cannot create %s", __func__, tmp_path);
    }
    setvbuf(src, NULL, _IOFBF, 8u << 20);
    setvbuf(dst, NULL, _IOFBF, 8u << 20);

    const size_t buf_size = 8u << 20;
    char *buf = (char *)malloc(buf_size);
    if (!buf) {
        fclose(src);
        fclose(dst);
        err(EXIT_FAILURE, "%s(): OOM payload copy buffer", __func__);
    }

    for (int i = 0; i < nfiles; ++i) {
        const uint64_t entries = index[i + 1] - index[i];
        const uint64_t bytes = remove_payload_bytes(entries, elem_size, suffix);
        if (remove_sample[i])
            remove_skip_exact(src, bytes, src_path, suffix);
        else
            remove_copy_exact(src, dst, bytes, buf, buf_size, src_path, tmp_path);
    }

    free(buf);
    if (fclose(src) != 0)
        err(errno, "%s(): failed closing %s", __func__, src_path);
    if (fclose(dst) != 0)
        err(errno, "%s(): failed closing %s", __func__, tmp_path);
    free(src_path);
}

static void remove_replace_tmp(const char *tmp_path, const char *target_dir,
                               const char *suffix)
{
    char *dst_path = append_join_path(target_dir, suffix);
    if (rename(tmp_path, dst_path) != 0)
        err(errno, "%s(): failed to replace %s", __func__, dst_path);
    free(dst_path);
}

static void remove_stale_sorted_index(const char *target_dir)
{
    char *path = append_join_path(target_dir, sorted_comb_ctxgid64obj32);
    if (unlink(path) != 0 && errno != ENOENT)
        err(errno, "%s(): cannot remove stale %s", __func__, path);
    free(path);
}

static void remove_file_if_exists(const char *target_dir, const char *suffix)
{
    char *path = append_join_path(target_dir, suffix);
    if (unlink(path) != 0 && errno != ENOENT)
        err(errno, "%s(): cannot remove stale %s", __func__, path);
    free(path);
}

static void remove_stale_optional_outputs(const char *target_dir,
                                          bool has_abundance, bool has_positions,
                                          bool has_qc, bool has_meta, bool has_anno)
{
    if (!has_abundance)
        remove_file_if_exists(target_dir, combined_ab_suffix);
    if (!has_positions)
        remove_file_if_exists(target_dir, sketch_position_suffix);
    if (!has_qc)
        remove_file_if_exists(target_dir, sketch_qc_stat);
    if (!has_meta)
        remove_file_if_exists(target_dir, sketch_infile_meta_stat);
    if (!has_anno)
        remove_file_if_exists(target_dir, sketch_anno_stat);
}

static void append_truncate_path(const char *path)
{
    FILE *fp = fopen(path, "wb");
    if (!fp)
        err(errno, "%s(): cannot create %s", __func__, path);
    if (fclose(fp) != 0)
        err(errno, "%s(): failed closing %s", __func__, path);
}

static void append_copy_parts_payload_to_tmp(const append_sketch_part_t *parts,
                                             int part_count,
                                             const char *suffix,
                                             const char *tmp_path)
{
    append_payload_rollback_t rollback = {0};
    append_truncate_path(tmp_path);
    for (int i = 0; i < part_count; ++i) {
        char *src = append_join_path(parts[i].path, suffix);
        append_copy_file_or_rollback(src, tmp_path, &rollback);
        free(src);
    }
}

static int apply_comblco_sample_filter(const char *input_dir, const char *output_dir,
                                       const append_sketch_part_t *target,
                                       const bool *remove_sample,
                                       bool drop_position,
                                       int *removed_samples_out,
                                       uint64_t *removed_entries_out)
{
    const bool keep_positions = target->has_positions && !drop_position;
    int removed_samples = 0;
    uint64_t removed_entries = 0;
    for (int i = 0; i < target->stat.infile_num; ++i) {
        if (remove_sample[i]) {
            ++removed_samples;
            removed_entries += target->index[i + 1] - target->index[i];
        }
    }
    if (removed_samples == target->stat.infile_num)
        errx(EINVAL, "%s(): sample filter would leave %s with zero samples; remove the sketch directory instead",
             __func__, input_dir);

    const int kept_samples = target->stat.infile_num - removed_samples;
    if (kept_samples < 0)
        errx(EINVAL, "%s(): sample filter produced a negative kept-sample count", __func__);
    const size_t kept_count = (size_t)kept_samples;
    uint64_t kept_entries = 0;
    uint64_t *new_index = (uint64_t *)calloc(kept_count + 1, sizeof(new_index[0]));
    char (*new_names)[PATHLEN] = kept_samples > 0 ? (char (*)[PATHLEN])calloc(kept_count, PATHLEN) : NULL;
    dim_sketch_qc_stat_t *new_qc =
        target->has_qc && kept_samples > 0 ? (dim_sketch_qc_stat_t *)calloc(kept_count, sizeof(new_qc[0])) : NULL;
    infile_meta_t *new_meta =
        target->has_meta && kept_samples > 0 ? (infile_meta_t *)calloc(kept_count, sizeof(new_meta[0])) : NULL;
    char (*new_anno)[PATHLEN] =
        target->has_anno && kept_samples > 0 ? (char (*)[PATHLEN])calloc(kept_count, PATHLEN) : NULL;
    if (!new_index || (kept_samples > 0 && !new_names) ||
        (target->has_qc && kept_samples > 0 && !new_qc) ||
        (target->has_meta && kept_samples > 0 && !new_meta) ||
        (target->has_anno && kept_samples > 0 && !new_anno))
        err(EXIT_FAILURE, "%s(): OOM filtered sketch auxiliary data", __func__);

    char (*target_names)[PATHLEN] = append_part_names(target);
    int out_i = 0;
    new_index[0] = 0;
    for (int i = 0; i < target->stat.infile_num; ++i) {
        if (remove_sample[i])
            continue;
        const uint64_t sample_entries = target->index[i + 1] - target->index[i];
        kept_entries += sample_entries;
        if (new_names)
            memcpy(new_names[out_i], target_names[i], PATHLEN);
        if (new_qc)
            new_qc[out_i] = target->qc[i];
        if (new_meta)
            new_meta[out_i] = target->meta[i];
        if (new_anno)
            memcpy(new_anno[out_i], target->anno[i], PATHLEN);
        new_index[out_i + 1] = kept_entries;
        ++out_i;
    }

    dim_sketch_stat_t new_stat = target->stat;
    new_stat.infile_num = kept_samples;

    char *tmp_comb = append_tmp_path(output_dir, combined_sketch_suffix);
    char *tmp_ab = target->stat.koc ? append_tmp_path(output_dir, combined_ab_suffix) : NULL;
    char *tmp_pos = keep_positions ? append_tmp_path(output_dir, sketch_position_suffix) : NULL;
    char *tmp_index = append_tmp_path(output_dir, idx_sketch_suffix);
    char *tmp_stat = append_tmp_path(output_dir, sketch_stat);
    char *tmp_qc = target->has_qc ? append_tmp_path(output_dir, sketch_qc_stat) : NULL;
    char *tmp_meta = target->has_meta ? append_tmp_path(output_dir, sketch_infile_meta_stat) : NULL;
    char *tmp_anno = target->has_anno ? append_tmp_path(output_dir, sketch_anno_stat) : NULL;

    remove_copy_filtered_payload(input_dir, combined_sketch_suffix, tmp_comb,
                                 target->index, target->stat.infile_num, remove_sample,
                                 sizeof(uint64_t));
    if (target->stat.koc)
        remove_copy_filtered_payload(input_dir, combined_ab_suffix, tmp_ab,
                                     target->index, target->stat.infile_num, remove_sample,
                                     sizeof(uint32_t));
    if (keep_positions)
        remove_copy_filtered_payload(input_dir, sketch_position_suffix, tmp_pos,
                                     target->index, target->stat.infile_num, remove_sample,
                                     sizeof(uint64_t));

    write_to_file(tmp_index, new_index, (kept_count + 1) * sizeof(new_index[0]));
    if (kept_samples > 0)
        concat_and_write_to_file(tmp_stat, &new_stat, sizeof(new_stat),
                                 new_names, kept_count * PATHLEN);
    else
        write_to_file(tmp_stat, &new_stat, sizeof(new_stat));
    if (target->has_qc)
        write_to_file(tmp_qc, new_qc ? (const void *)new_qc : "", kept_count * sizeof(new_qc[0]));
    if (target->has_meta)
        write_to_file(tmp_meta, new_meta ? (const void *)new_meta : "", kept_count * sizeof(new_meta[0]));
    if (target->has_anno)
        write_to_file(tmp_anno, new_anno ? (const void *)new_anno : "", kept_count * PATHLEN);

    remove_stale_sorted_index(output_dir);
    remove_replace_tmp(tmp_comb, output_dir, combined_sketch_suffix);
    if (target->stat.koc)
        remove_replace_tmp(tmp_ab, output_dir, combined_ab_suffix);
    if (keep_positions)
        remove_replace_tmp(tmp_pos, output_dir, sketch_position_suffix);
    remove_replace_tmp(tmp_index, output_dir, idx_sketch_suffix);
    remove_replace_tmp(tmp_stat, output_dir, sketch_stat);
    if (target->has_qc)
        remove_replace_tmp(tmp_qc, output_dir, sketch_qc_stat);
    if (target->has_meta)
        remove_replace_tmp(tmp_meta, output_dir, sketch_infile_meta_stat);
    if (target->has_anno)
        remove_replace_tmp(tmp_anno, output_dir, sketch_anno_stat);
    remove_stale_optional_outputs(output_dir, target->stat.koc, keep_positions,
                                  target->has_qc, target->has_meta, target->has_anno);

    free(tmp_comb);
    free(tmp_ab);
    free(tmp_pos);
    free(tmp_index);
    free(tmp_stat);
    free(tmp_qc);
    free(tmp_meta);
    free(tmp_anno);
    free(new_index);
    free(new_names);
    free(new_qc);
    free(new_meta);
    free(new_anno);

    if (removed_samples_out)
        *removed_samples_out = removed_samples;
    if (removed_entries_out)
        *removed_entries_out = removed_entries;
    return kept_samples;
}

static int filter_comblco_samples(sketch_opt_t *sketch_opt_val, bool keep_mode)
{
    const char *op = keep_mode ? "--keep" : "--remove";
    const char *op_name = keep_mode ? "keep" : "remove";
    const char *list_path = keep_mode ? sketch_opt_val->keep_list : sketch_opt_val->remove_list;
    const char *source_path = keep_mode ? sketch_opt_val->keep_source : sketch_opt_val->remove_source;
    const bool copy_mode = keep_mode ? sketch_opt_val->keep_copy_mode : sketch_opt_val->remove_copy_mode;

    if (!list_path || list_path[0] == '\0')
        errx(EINVAL, "%s(): %s requires a sample-name list file", __func__, op);

    const char *input_dir = source_path
                                ? source_path
                                : sketch_opt_val->outdir;
    const char *output_dir = sketch_opt_val->outdir;
    if (copy_mode) {
        if (append_same_existing_sketch_dir(output_dir, input_dir))
            errx(EINVAL, "%s(): %s copy-mode output %s is the same directory as input %s; omit -o for in-place filtering",
                 __func__, op, output_dir, input_dir);
        mkdir_p(output_dir);
    }

    append_sketch_part_t target = {0};
    append_load_part(&target, input_dir);

    size_t remove_name_count = 0;
    remove_name_entry_t *remove_names =
        remove_load_name_list(list_path, op_name, &remove_name_count);
    bool *remove_sample = (bool *)calloc((size_t)target.stat.infile_num, sizeof(remove_sample[0]));
    if (target.stat.infile_num > 0 && !remove_sample)
        err(EXIT_FAILURE, "%s(): OOM sample filter flags", __func__);

    char (*target_names)[PATHLEN] = append_part_names(&target);
    int removed_samples = 0;
    uint64_t removed_entries = 0;
    int listed_samples = 0;
    for (int i = 0; i < target.stat.infile_num; ++i) {
        char sample_name[PATHLEN];
        memcpy(sample_name, target_names[i], PATHLEN);
        sample_name[PATHLEN - 1] = '\0';
        size_t list_idx = 0;
        const bool listed = remove_find_name(remove_names, remove_name_count, sample_name, &list_idx);
        if (listed) {
            remove_names[list_idx].found = true;
            ++listed_samples;
        }
        const bool should_remove = keep_mode ? !listed : listed;
        remove_sample[i] = should_remove;
        if (should_remove) {
            ++removed_samples;
            removed_entries += target.index[i + 1] - target.index[i];
        }
    }

    for (size_t i = 0; i < remove_name_count; ++i) {
        if (!remove_names[i].found)
            errx(EINVAL, "%s(): sample listed for %s was not found in %s: %s",
                 __func__, op_name, input_dir, remove_names[i].name);
    }
    if (keep_mode && listed_samples == 0)
        errx(EINVAL, "%s(): no samples selected to keep", __func__);
    if (!keep_mode && removed_samples == 0)
        errx(EINVAL, "%s(): no samples selected for removal", __func__);
    if (removed_samples == target.stat.infile_num)
        errx(EINVAL, "%s(): %s would leave %s with zero samples; remove the sketch directory instead",
             __func__, op, input_dir);

    const int kept_samples = apply_comblco_sample_filter(input_dir, output_dir,
                                                         &target, remove_sample,
                                                         sketch_opt_val->drop_position,
                                                         &removed_samples,
                                                         &removed_entries);

    if (keep_mode) {
        if (copy_mode)
            printf("Kept %d samples from %s into %s; removed %d samples (%" PRIu64 " sketch entries); total samples=%d\n",
                   kept_samples, input_dir, output_dir, removed_samples, removed_entries, kept_samples);
        else
            printf("Kept %d samples in %s; removed %d samples (%" PRIu64 " sketch entries); total samples=%d\n",
                   kept_samples, output_dir, removed_samples, removed_entries, kept_samples);
    }
    else {
        if (copy_mode)
            printf("Removed %d samples (%" PRIu64 " sketch entries) from %s into %s; total samples=%d\n",
                   removed_samples, removed_entries, input_dir, output_dir, kept_samples);
        else
            printf("Removed %d samples (%" PRIu64 " sketch entries) from %s; total samples=%d\n",
                   removed_samples, removed_entries, output_dir, kept_samples);
    }

    free(remove_sample);
    remove_free_name_entries(remove_names, remove_name_count);
    append_free_part(&target);

    return kept_samples;
}

static const char *dedup_metric_name(sketch_dedup_metric_t metric)
{
    return pairwise_metric_expr_name(&metric);
}

typedef pairwise_eval_t dedup_pair_eval_t;

static uint32_t dedup_count_ctx_runs_sorted_ctxobj64(const uint64_t *a, size_t n)
{
    return pairwise_count_ctx_runs_sorted_ctxobj64(a, n);
}

static bool dedup_meta_rankable(const append_sketch_part_t *target, int idx)
{
    if (!target->has_meta || idx < 0 || idx >= target->stat.infile_num)
        return false;
    const infile_meta_t *meta = &target->meta[idx];
    return meta->meta_fmt_version == KSSD_INFILE_META_VERSION &&
           meta->total_length_bp > 0;
}

static double dedup_finite_or(double value, double fallback)
{
    return isfinite(value) ? value : fallback;
}

static int dedup_positive_lower_u32_cmp(uint32_t a, uint32_t b)
{
    if (a == b)
        return 0;
    if (a == 0)
        return -1;
    if (b == 0)
        return 1;
    return a < b ? 1 : -1;
}

static int dedup_quality_compare(const append_sketch_part_t *target, int a, int b)
{
    if (a == b)
        return 0;

    const bool a_rankable = dedup_meta_rankable(target, a);
    const bool b_rankable = dedup_meta_rankable(target, b);
    if (a_rankable != b_rankable)
        return a_rankable ? 1 : -1;

    if (a_rankable && b_rankable) {
        const infile_meta_t *ma = &target->meta[a];
        const infile_meta_t *mb = &target->meta[b];
        const bool a_fasta = ma->infile_fmt == KSSD_INFILE_FMT_FASTA;
        const bool b_fasta = mb->infile_fmt == KSSD_INFILE_FMT_FASTA;
        if (a_fasta != b_fasta)
            return a_fasta ? 1 : -1;

        const double a_asm = dedup_finite_or((double)ma->asm_level, -1.0);
        const double b_asm = dedup_finite_or((double)mb->asm_level, -1.0);
        if (a_asm != b_asm)
            return a_asm > b_asm ? 1 : -1;

        if (ma->total_length_bp != mb->total_length_bp)
            return ma->total_length_bp > mb->total_length_bp ? 1 : -1;

        const int record_cmp = dedup_positive_lower_u32_cmp(ma->record_count, mb->record_count);
        if (record_cmp != 0)
            return record_cmp;

        if (ma->median_length_bp != mb->median_length_bp)
            return ma->median_length_bp > mb->median_length_bp ? 1 : -1;

        const double a_cv = dedup_finite_or((double)ma->length_cv, INFINITY);
        const double b_cv = dedup_finite_or((double)mb->length_cv, INFINITY);
        if (a_cv != b_cv)
            return a_cv < b_cv ? 1 : -1;
    }

    const uint64_t a_entries = target->index[a + 1] - target->index[a];
    const uint64_t b_entries = target->index[b + 1] - target->index[b];
    if (a_entries != b_entries)
        return a_entries > b_entries ? 1 : -1;

    return a < b ? 1 : -1;
}

typedef struct dedup_component_ctx
{
    sketch_dedup_metric_t metric;
    const dim_sketch_stat_t *stat;
    const uint64_t *comb;
    const uint64_t *index;
    const uint32_t *ctx_counts;
    const append_sketch_part_t *target;
} dedup_component_ctx_t;

static dedup_pair_eval_t dedup_component_eval_pair(void *ctx, int a, int b)
{
    const dedup_component_ctx_t *dctx = ctx;
    const uint64_t a_begin = dctx->index[a];
    const uint64_t a_end = dctx->index[a + 1];
    const uint64_t b_begin = dctx->index[b];
    const uint64_t b_end = dctx->index[b + 1];
    return pairwise_eval_expr_arrays(&dctx->metric, dctx->stat,
                                     dctx->comb + a_begin, (size_t)(a_end - a_begin),
                                     dctx->comb + b_begin, (size_t)(b_end - b_begin),
                                     dctx->ctx_counts[a], dctx->ctx_counts[b],
                                     false);
}

static int dedup_component_quality_compare(void *ctx, int a, int b)
{
    const dedup_component_ctx_t *dctx = ctx;
    return dedup_quality_compare(dctx->target, a, b);
}

static void dedup_index_edge_observer(void *ctx, int qry, int ref,
                                      const pairwise_eval_t *eval)
{
    (void)eval;
    pairwise_edge_list_add((pairwise_edge_list_t *)ctx, qry, ref);
}

#define SKETCH_DEDUP_PROGRESS_AUTO_MIN 1000ULL

int dedup_comblco_samples(sketch_opt_t *sketch_opt_val)
{
    const char *input_dir = sketch_opt_val->dedup_source
                                ? sketch_opt_val->dedup_source
                                : sketch_opt_val->outdir;
    const char *output_dir = sketch_opt_val->outdir;
    const bool copy_mode = sketch_opt_val->dedup_copy_mode;

    if (copy_mode) {
        if (append_same_existing_sketch_dir(output_dir, input_dir))
            errx(EINVAL, "%s(): --dedup copy-mode output %s is the same directory as input %s; omit -o for in-place dedup",
                 __func__, output_dir, input_dir);
        mkdir_p(output_dir);
    }

    append_sketch_part_t target = {0};
    append_load_part(&target, input_dir);
    if (target.stat.infile_num <= 0)
        errx(EINVAL, "%s(): %s contains no samples", __func__, input_dir);

    char *comb_path = append_join_path(input_dir, combined_sketch_suffix);
    size_t comb_size = 0;
    uint64_t *comb = read_from_file(comb_path, &comb_size);
    free(comb_path);
    if (comb_size != target.comblco_size)
        errx(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
             __func__, input_dir, combined_sketch_suffix, comb_size, target.comblco_size);

    const int nfiles = target.stat.infile_num;
    const_comask_init(&target.stat);
    ani_model_drfold = target.stat.drfold;

    uint32_t *ctx_counts = NULL;
    if (!sketch_opt_val->dedup_index) {
        ctx_counts = (uint32_t *)calloc((size_t)nfiles, sizeof(ctx_counts[0]));
        if (!ctx_counts)
            err(EXIT_FAILURE, "%s(): OOM dedup bookkeeping", __func__);
        for (int i = 0; i < nfiles; ++i) {
            const uint64_t begin = target.index[i];
            const uint64_t end = target.index[i + 1];
            ctx_counts[i] = dedup_count_ctx_runs_sorted_ctxobj64(comb + begin, (size_t)(end - begin));
        }
    }

    const double guard_max_afcut = sketch_opt_val->dedup_max_afcut;
    const uint32_t guard_ctxcut = sketch_opt_val->dedup_ctxcut;
    dedup_component_ctx_t ctx = {
        .metric = sketch_opt_val->dedup_metric,
        .stat = &target.stat,
        .comb = comb,
        .index = target.index,
        .ctx_counts = ctx_counts,
        .target = &target,
    };
    pairwise_component_result_t comp = {0};
    pairwise_index_scan_stats_t index_stats = {0};
    if (sketch_opt_val->dedup_index) {
        if (target.stat.conflict)
            errx(EINVAL, "%s(): --dedup-index requires a non-conflict lco sketch", __func__);
        if (!append_file_exists(input_dir, sorted_comb_ctxgid64obj32))
            errx(EINVAL, "%s(): --dedup-index requires %s/%s; build it with 'kssd3a sketch -i %s'",
                 __func__, input_dir, sorted_comb_ctxgid64obj32, input_dir);

        unify_sketch_t sketch_view = {0};
        sketch_view.stat_type = 2;
        sketch_view.hash_id = target.stat.hash_id;
        sketch_view.kmerlen = target.stat.klen;
        sketch_view.infile_num = nfiles;
        sketch_view.conflict = target.stat.conflict;
        sketch_view.comb_sketch = comb;
        sketch_view.sketch_index = target.index;
        sketch_view.stats.lco_stat_val = target.stat;

        pairwise_edge_list_t dedup_edges = {0};
        pairwise_index_scan_options_t scan_opt = {
            .metric = &sketch_opt_val->dedup_metric,
            .cut = sketch_opt_val->dedup_cutoff,
            .ctxcut = guard_ctxcut,
            .max_afcut = guard_max_afcut,
            .index_max_ctx_freq = sketch_opt_val->dedup_index_max_ctx_freq,
            .index_min_votes = sketch_opt_val->dedup_index_min_votes,
            .index_sample_step = sketch_opt_val->dedup_index_sample_step,
            .threads = sketch_opt_val->p,
        };
        char *sorted_index_path = append_join_path(input_dir, sorted_comb_ctxgid64obj32);
        kssd_progress_t progress =
            kssd_progress_start((uint64_t)nfiles >= SKETCH_DEDUP_PROGRESS_AUTO_MIN,
                                "sketch", "indexed dedup rows", "rows",
                                (uint64_t)nfiles);
        pairwise_indexed_self_scan(&sketch_view, sorted_index_path, &scan_opt,
                                   dedup_index_edge_observer, &dedup_edges,
                                   kssd_progress_cb, &progress,
                                   &index_stats);
        kssd_progress_done(&progress);
        free(sorted_index_path);
        if (sketch_opt_val->dedup_strategy == PAIRWISE_DEDUP_COMPLETE_LINKAGE) {
            pairwise_complete_linkage_dedup_from_edges(nfiles, dedup_edges.edges, dedup_edges.n,
                                                       &ctx, dedup_component_quality_compare,
                                                       &comp);
        } else {
            pairwise_greedy_dedup_from_edges(nfiles, dedup_edges.edges, dedup_edges.n,
                                             &ctx, dedup_component_quality_compare,
                                             &comp);
        }
        pairwise_edge_list_free(&dedup_edges);
        comp.distance_edges = index_stats.distance_edges;
        comp.duplicate_edges = index_stats.accepted_edges;
        comp.ctx_rejects = index_stats.ctx_rejects;
        comp.max_af_rejects = index_stats.max_af_rejects;
    } else {
        if (sketch_opt_val->dedup_strategy == PAIRWISE_DEDUP_COMPLETE_LINKAGE) {
            pairwise_build_complete_linkage_dedup(nfiles, sketch_opt_val->dedup_cutoff,
                                                  guard_ctxcut, guard_max_afcut, &ctx,
                                                  dedup_component_eval_pair,
                                                  dedup_component_quality_compare,
                                                  NULL, NULL, NULL, NULL, &comp);
        } else {
            pairwise_build_greedy_dedup(nfiles, sketch_opt_val->dedup_cutoff,
                                        guard_ctxcut, guard_max_afcut, &ctx,
                                        dedup_component_eval_pair,
                                        dedup_component_quality_compare,
                                        NULL, NULL, NULL, NULL, &comp);
        }
    }

    uint64_t removed_entries = 0;
    int kept_samples = nfiles;
    int removed_samples = comp.removed_samples;
    if (removed_samples > 0 || copy_mode || sketch_opt_val->drop_position) {
        kept_samples = apply_comblco_sample_filter(input_dir, output_dir,
                                                   &target, comp.remove_sample,
                                                   sketch_opt_val->drop_position,
                                                   &removed_samples,
                                                   &removed_entries);
    }

    if (copy_mode)
        printf("Deduplicated %s into %s; metric=%s strategy=%s cutoff=%.12g dedup_max_afcut=%.12g dedup_ctxcut=%u dedup_index=%d index_max_ctx_freq=%u index_min_votes=%u index_sample_step=%u index_candidate_pairs=%" PRIu64 " index_exact_pairs=%" PRIu64 " kept=%d removed=%d duplicate_clusters=%d duplicate_edges=%" PRIu64 " distance_edges=%" PRIu64 " dedup_ctx_rejects=%" PRIu64 " dedup_max_af_rejects=%" PRIu64 " removed_entries=%" PRIu64 "\n",
               input_dir, output_dir, dedup_metric_name(sketch_opt_val->dedup_metric),
               pairwise_dedup_strategy_name(sketch_opt_val->dedup_strategy),
               sketch_opt_val->dedup_cutoff, guard_max_afcut, guard_ctxcut,
               sketch_opt_val->dedup_index ? 1 : 0,
               sketch_opt_val->dedup_index_max_ctx_freq,
               sketch_opt_val->dedup_index_min_votes,
               sketch_opt_val->dedup_index_sample_step,
               index_stats.candidate_pairs, index_stats.exact_pairs,
               kept_samples, removed_samples,
               comp.duplicate_clusters, comp.duplicate_edges, comp.distance_edges,
               comp.ctx_rejects, comp.max_af_rejects, removed_entries);
    else
        printf("Deduplicated %s; metric=%s strategy=%s cutoff=%.12g dedup_max_afcut=%.12g dedup_ctxcut=%u dedup_index=%d index_max_ctx_freq=%u index_min_votes=%u index_sample_step=%u index_candidate_pairs=%" PRIu64 " index_exact_pairs=%" PRIu64 " kept=%d removed=%d duplicate_clusters=%d duplicate_edges=%" PRIu64 " distance_edges=%" PRIu64 " dedup_ctx_rejects=%" PRIu64 " dedup_max_af_rejects=%" PRIu64 " removed_entries=%" PRIu64 "\n",
               output_dir, dedup_metric_name(sketch_opt_val->dedup_metric),
               pairwise_dedup_strategy_name(sketch_opt_val->dedup_strategy),
               sketch_opt_val->dedup_cutoff, guard_max_afcut, guard_ctxcut,
               sketch_opt_val->dedup_index ? 1 : 0,
               sketch_opt_val->dedup_index_max_ctx_freq,
               sketch_opt_val->dedup_index_min_votes,
               sketch_opt_val->dedup_index_sample_step,
               index_stats.candidate_pairs, index_stats.exact_pairs,
               kept_samples, removed_samples,
               comp.duplicate_clusters, comp.duplicate_edges, comp.distance_edges,
               comp.ctx_rejects, comp.max_af_rejects, removed_entries);

    pairwise_component_result_free(&comp);
    free(ctx_counts);
    free_read_from_file(comb, comb_size);
    append_free_part(&target);

    return kept_samples;
}

int remove_comblco_samples(sketch_opt_t *sketch_opt_val)
{
    return filter_comblco_samples(sketch_opt_val, false);
}

int keep_comblco_samples(sketch_opt_t *sketch_opt_val)
{
    return filter_comblco_samples(sketch_opt_val, true);
}

int append_comblco(sketch_opt_t *sketch_opt_val)
{
    const bool copy_mode = sketch_opt_val->append_copy_mode;
    if (copy_mode) {
        if (sketch_opt_val->num_remaining_args < 2)
            errx(EINVAL, "%s(): --append copy mode requires a base sketch and at least one source sketch", __func__);
        for (int i = 0; i < sketch_opt_val->num_remaining_args; ++i) {
            const char *src = sketch_opt_val->remaining_args[i];
            if (append_same_existing_sketch_dir(sketch_opt_val->outdir, src))
                errx(EINVAL, "%s(): --append copy-mode output %s is the same directory as input %s; omit -o for in-place append",
                     __func__, sketch_opt_val->outdir, src);
        }
        mkdir_p(sketch_opt_val->outdir);
    } else if (sketch_opt_val->num_remaining_args <= 0) {
        errx(EINVAL, "%s(): --append requires at least one source sketch directory", __func__);
    }

    const int part_count = copy_mode ? sketch_opt_val->num_remaining_args
                                     : sketch_opt_val->num_remaining_args + 1;
    append_sketch_part_t *parts = calloc((size_t)part_count, sizeof(parts[0]));
    if (!parts)
        err(EXIT_FAILURE, "%s(): OOM append parts", __func__);

    if (copy_mode) {
        for (int i = 0; i < part_count; ++i)
            append_load_part(&parts[i], sketch_opt_val->remaining_args[i]);
    } else {
        append_load_part(&parts[0], sketch_opt_val->outdir);
        for (int i = 1; i < part_count; ++i) {
            const char *src = sketch_opt_val->remaining_args[i - 1];
            if (append_same_sketch_dir(sketch_opt_val->outdir, src))
                errx(EINVAL, "%s(): source sketch %s is the same directory as target %s",
                     __func__, src, sketch_opt_val->outdir);
            append_load_part(&parts[i], src);
        }
    }

    for (int i = 1; i < part_count; ++i) {
        append_check_compatible_stat(&parts[0].stat, &parts[i].stat, parts[i].path);
        if (parts[i].has_positions != parts[0].has_positions)
            errx(EINVAL, "%s(): mixed %s sidecars are not supported by --append "
                 "(base has %s, source %s has %s)",
                 __func__, sketch_position_suffix,
                 parts[0].has_positions ? "positions" : "no positions",
                 parts[i].path,
                 parts[i].has_positions ? "positions" : "no positions");
    }

    int total_samples = 0;
    uint64_t total_entries = 0;
    uint64_t appended_entries = 0;
    int appended_samples = 0;
    bool any_qc = false;
    bool any_meta = false;
    bool any_anno = false;
    bool seen_meta = false;
    bool seen_no_meta = false;
    for (int i = 0; i < part_count; ++i) {
        if (parts[i].stat.infile_num > INT_MAX - total_samples)
            errx(EINVAL, "%s(): appended sample count exceeds INT_MAX", __func__);
        total_samples += parts[i].stat.infile_num;
        if (UINT64_MAX - total_entries < parts[i].entries)
            errx(EINVAL, "%s(): appended entry count overflows uint64_t", __func__);
        total_entries += parts[i].entries;
        any_qc = any_qc || parts[i].has_qc;
        any_meta = any_meta || parts[i].has_meta;
        any_anno = any_anno || parts[i].has_anno;
        seen_meta = seen_meta || parts[i].has_meta;
        seen_no_meta = seen_no_meta || !parts[i].has_meta;
        if (i > 0) {
            appended_samples += parts[i].stat.infile_num;
            appended_entries += parts[i].entries;
        }
    }

    uint64_t *merged_index = calloc((size_t)total_samples + 1, sizeof(merged_index[0]));
    char (*merged_names)[PATHLEN] = calloc((size_t)total_samples, PATHLEN);
    dim_sketch_qc_stat_t *merged_qc =
        any_qc ? calloc((size_t)total_samples, sizeof(merged_qc[0])) : NULL;
    infile_meta_t *merged_meta =
        any_meta ? calloc((size_t)total_samples, sizeof(merged_meta[0])) : NULL;
    char (*merged_anno)[PATHLEN] = any_anno ? calloc((size_t)total_samples, PATHLEN) : NULL;
    if (!merged_index || !merged_names || (any_qc && !merged_qc) ||
        (any_meta && !merged_meta) || (any_anno && !merged_anno))
        err(EXIT_FAILURE, "%s(): OOM merged append auxiliary data", __func__);

    int sample_offset = 0;
    uint64_t entry_offset = 0;
    merged_index[0] = 0;
    for (int i = 0; i < part_count; ++i) {
        const append_sketch_part_t *part = &parts[i];
        memcpy(merged_names + sample_offset, append_part_names(part),
               (size_t)part->stat.infile_num * PATHLEN);
        for (int j = 1; j <= part->stat.infile_num; ++j)
            merged_index[sample_offset + j] = entry_offset + part->index[j];
        if (part->has_qc)
            memcpy(merged_qc + sample_offset, part->qc,
                   (size_t)part->stat.infile_num * sizeof(merged_qc[0]));
        if (part->has_meta)
            memcpy(merged_meta + sample_offset, part->meta,
                   (size_t)part->stat.infile_num * sizeof(merged_meta[0]));
        if (part->has_anno)
            memcpy(merged_anno + sample_offset, part->anno,
                   (size_t)part->stat.infile_num * PATHLEN);
        sample_offset += part->stat.infile_num;
        entry_offset += part->entries;
    }
    if (any_meta && seen_meta && seen_no_meta)
        memset(merged_meta, 0, (size_t)total_samples * sizeof(merged_meta[0]));

    dim_sketch_stat_t merged_stat = parts[0].stat;
    merged_stat.infile_num = total_samples;

    char *tmp_index = append_tmp_path(sketch_opt_val->outdir, idx_sketch_suffix);
    char *tmp_stat = append_tmp_path(sketch_opt_val->outdir, sketch_stat);
    char *tmp_qc = any_qc ? append_tmp_path(sketch_opt_val->outdir, sketch_qc_stat) : NULL;
    char *tmp_meta = any_meta ? append_tmp_path(sketch_opt_val->outdir, sketch_infile_meta_stat) : NULL;
    char *tmp_anno = any_anno ? append_tmp_path(sketch_opt_val->outdir, sketch_anno_stat) : NULL;
    char *tmp_comb = copy_mode ? append_tmp_path(sketch_opt_val->outdir, combined_sketch_suffix) : NULL;
    char *tmp_ab = copy_mode && parts[0].stat.koc ? append_tmp_path(sketch_opt_val->outdir, combined_ab_suffix) : NULL;
    char *tmp_pos = copy_mode && parts[0].has_positions ? append_tmp_path(sketch_opt_val->outdir, sketch_position_suffix) : NULL;

    write_to_file(tmp_index, merged_index,
                  ((size_t)total_samples + 1) * sizeof(merged_index[0]));
    concat_and_write_to_file(tmp_stat, &merged_stat, sizeof(merged_stat),
                             merged_names, (size_t)total_samples * PATHLEN);
    if (any_qc)
        write_to_file(tmp_qc, merged_qc, (size_t)total_samples * sizeof(merged_qc[0]));
    if (any_meta)
        write_to_file(tmp_meta, merged_meta, (size_t)total_samples * sizeof(merged_meta[0]));
    if (any_anno)
        write_to_file(tmp_anno, merged_anno, (size_t)total_samples * PATHLEN);

    append_payload_rollback_t rollback = {0};
    if (copy_mode) {
        append_copy_parts_payload_to_tmp(parts, part_count, combined_sketch_suffix, tmp_comb);
        if (parts[0].stat.koc)
            append_copy_parts_payload_to_tmp(parts, part_count, combined_ab_suffix, tmp_ab);
        if (parts[0].has_positions)
            append_copy_parts_payload_to_tmp(parts, part_count, sketch_position_suffix, tmp_pos);

        remove_stale_sorted_index(sketch_opt_val->outdir);
        remove_replace_tmp(tmp_comb, sketch_opt_val->outdir, combined_sketch_suffix);
        if (parts[0].stat.koc)
            remove_replace_tmp(tmp_ab, sketch_opt_val->outdir, combined_ab_suffix);
        if (parts[0].has_positions)
            remove_replace_tmp(tmp_pos, sketch_opt_val->outdir, sketch_position_suffix);
        remove_replace_tmp(tmp_index, sketch_opt_val->outdir, idx_sketch_suffix);
        remove_replace_tmp(tmp_stat, sketch_opt_val->outdir, sketch_stat);
        if (any_qc)
            remove_replace_tmp(tmp_qc, sketch_opt_val->outdir, sketch_qc_stat);
        if (any_meta)
            remove_replace_tmp(tmp_meta, sketch_opt_val->outdir, sketch_infile_meta_stat);
        if (any_anno)
            remove_replace_tmp(tmp_anno, sketch_opt_val->outdir, sketch_anno_stat);
        remove_stale_optional_outputs(sketch_opt_val->outdir, parts[0].stat.koc,
                                      parts[0].has_positions, any_qc, any_meta, any_anno);
    } else {
        rollback.comblco_path = append_join_path(sketch_opt_val->outdir, combined_sketch_suffix);
        rollback.comblco_size = parts[0].comblco_size;
        rollback.has_abundance = parts[0].stat.koc;
        rollback.has_positions = parts[0].has_positions;
        if (rollback.has_abundance) {
            rollback.abundance_path = append_join_path(sketch_opt_val->outdir, combined_ab_suffix);
            rollback.abundance_size = parts[0].abundance_size;
        }
        if (rollback.has_positions) {
            rollback.position_path = append_join_path(sketch_opt_val->outdir, sketch_position_suffix);
            rollback.position_size = parts[0].position_size;
        }

        for (int i = 1; i < part_count; ++i) {
            char *src_comb = append_join_path(parts[i].path, combined_sketch_suffix);
            append_copy_file_or_rollback(src_comb, rollback.comblco_path, &rollback);
            free(src_comb);
            if (rollback.has_abundance) {
                char *src_ab = append_join_path(parts[i].path, combined_ab_suffix);
                append_copy_file_or_rollback(src_ab, rollback.abundance_path, &rollback);
                free(src_ab);
            }
            if (rollback.has_positions) {
                char *src_pos = append_join_path(parts[i].path, sketch_position_suffix);
                append_copy_file_or_rollback(src_pos, rollback.position_path, &rollback);
                free(src_pos);
            }
        }

        remove_stale_sorted_index(sketch_opt_val->outdir);
        append_replace_tmp_or_rollback(tmp_index, sketch_opt_val->outdir, idx_sketch_suffix, &rollback);
        append_replace_tmp_or_rollback(tmp_stat, sketch_opt_val->outdir, sketch_stat, &rollback);
        if (any_qc)
            append_replace_tmp_or_rollback(tmp_qc, sketch_opt_val->outdir, sketch_qc_stat, &rollback);
        if (any_meta)
            append_replace_tmp_or_rollback(tmp_meta, sketch_opt_val->outdir, sketch_infile_meta_stat, &rollback);
        if (any_anno)
            append_replace_tmp_or_rollback(tmp_anno, sketch_opt_val->outdir, sketch_anno_stat, &rollback);
    }

    if (copy_mode)
        printf("Appended %d samples (%" PRIu64 " sketch entries) from %s into %s; total samples=%d\n",
               appended_samples, appended_entries, parts[0].path, sketch_opt_val->outdir, total_samples);
    else
        printf("Appended %d samples (%" PRIu64 " sketch entries) into %s; total samples=%d\n",
               appended_samples, appended_entries, sketch_opt_val->outdir, total_samples);

    free(tmp_index);
    free(tmp_stat);
    free(tmp_qc);
    free(tmp_meta);
    free(tmp_anno);
    free(tmp_comb);
    free(tmp_ab);
    free(tmp_pos);
    free(rollback.comblco_path);
    free(rollback.abundance_path);
    free(rollback.position_path);
    free(merged_index);
    free(merged_names);
    free(merged_qc);
    free(merged_meta);
    free(merged_anno);
    for (int i = 0; i < part_count; ++i)
        append_free_part(&parts[i]);
    free(parts);

    return total_samples;
}

int sketch_qc_comblco(sketch_opt_t *sketch_opt_val)
{
    if (sketch_opt_val->num_remaining_args != 1)
        err(EINVAL, "%s(): --sketchQC requires exactly one input sketch directory", __func__);

    const char *indir = sketch_opt_val->remaining_args[0];
    unify_sketch_t *in_sketch = generic_sketch_parse(
        indir, SKETCH_PARSE_ABUNDANCE | SKETCH_PARSE_SAMPLE_QC);
    if (in_sketch->stat_type != 2)
        err(EINVAL, "%s(): --sketchQC only supports 64-bit comblco sketches", __func__);
    if (!in_sketch->stats.lco_stat_val.koc || !in_sketch->abundance)
        err(EINVAL, "%s(): --sketchQC requires an abundance sketch with %s",
            __func__, combined_ab_suffix);
    const bool infer_missing_qc = in_sketch->sample_qc == NULL;

    mkdir_p(sketch_opt_val->outdir);

    FILE *comb = fopen(test_create_fullpath(sketch_opt_val->outdir, combined_sketch_suffix), "wb");
    if (!comb)
        err(errno, "%s(): cannot create %s/%s", __func__, sketch_opt_val->outdir,
            combined_sketch_suffix);
    setvbuf(comb, NULL, _IOFBF, 8u << 20);

    FILE *comb_ab = fopen(test_create_fullpath(sketch_opt_val->outdir, combined_ab_suffix), "wb");
    if (!comb_ab)
        err(errno, "%s(): cannot create %s/%s", __func__, sketch_opt_val->outdir,
            combined_ab_suffix);
    setvbuf(comb_ab, NULL, _IOFBF, 8u << 20);

    FILE *comb_pos = NULL;
    if (in_sketch->positions) {
        comb_pos = fopen(test_create_fullpath(sketch_opt_val->outdir, sketch_position_suffix), "wb");
        if (!comb_pos)
            err(errno, "%s(): cannot create %s/%s", __func__, sketch_opt_val->outdir,
                sketch_position_suffix);
        setvbuf(comb_pos, NULL, _IOFBF, 8u << 20);
    } else {
        remove_sketch_positions(sketch_opt_val->outdir);
    }

    const int nfiles = in_sketch->infile_num;
    infile_meta_t *in_infile_meta = NULL;
    if (file_exists_in_folder(indir, sketch_infile_meta_stat)) {
        size_t meta_file_size = 0;
        char *meta_path = test_get_fullpath(indir, sketch_infile_meta_stat);
        in_infile_meta = read_from_file(meta_path, &meta_file_size);
        free(meta_path);
        const size_t expected_size = sizeof(in_infile_meta[0]) * (size_t)nfiles;
        if (meta_file_size != expected_size)
            err(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
                __func__, indir, sketch_infile_meta_stat, meta_file_size, expected_size);
    }
    uint64_t *out_index = (uint64_t *)calloc((size_t)nfiles + 1, sizeof(out_index[0]));
    dim_sketch_qc_stat_t *out_qc =
        (dim_sketch_qc_stat_t *)calloc((size_t)nfiles, sizeof(out_qc[0]));
    if (!out_index || !out_qc)
        err(errno, "%s(): OOM output index/QC buffers", __func__);

    uint64_t total_in = 0;
    uint64_t total_kept = 0;
    int invalid_ranges = 0;
    for (int sample = 0; sample < nfiles; ++sample) {
        const uint64_t start = in_sketch->sketch_index[sample];
        const uint64_t end = in_sketch->sketch_index[sample + 1];
        dim_sketch_qc_stat_t qc = {0};
        if (infer_missing_qc) {
            SortedKV_Arrays_t sample_kv = {
                .keys = in_sketch->comb_sketch + start,
                .values = in_sketch->abundance + start,
                .len = end - start,
            };
            const reads_qc_range_t range = infer_reads_qc_count_range(&sample_kv, 2u);
            qc = sample_qc_from_range(&range, false);
        } else {
            qc = in_sketch->sample_qc[sample];
        }
        const bool valid_range = (qc.flags & DIM_SKETCH_QC_RANGE_VALID) != 0;
        out_qc[sample] = qc;
        if (valid_range)
            out_qc[sample].flags |= DIM_SKETCH_QC_RANGE_APPLIED;
        else
            ++invalid_ranges;

        for (uint64_t idx = start; idx < end; ++idx) {
            const uint32_t count = in_sketch->abundance[idx];
            const bool keep = !valid_range
                              || (count >= qc.reads_qc_lower
                                  && (!qc.reads_qc_upper || count <= qc.reads_qc_upper));
            if (!keep)
                continue;

            if (fwrite(&in_sketch->comb_sketch[idx], sizeof(in_sketch->comb_sketch[idx]), 1, comb) != 1)
                err(errno, "%s(): write %s failed", __func__, combined_sketch_suffix);
            if (fwrite(&in_sketch->abundance[idx], sizeof(in_sketch->abundance[idx]), 1, comb_ab) != 1)
                err(errno, "%s(): write %s failed", __func__, combined_ab_suffix);
            if (comb_pos &&
                fwrite(&in_sketch->positions[idx], sizeof(in_sketch->positions[idx]), 1, comb_pos) != 1)
                err(errno, "%s(): write %s failed", __func__, sketch_position_suffix);
            ++total_kept;
        }
        total_in += end - start;
        out_index[sample + 1] = total_kept;
    }

    fclose(comb);
    fclose(comb_ab);
    if (comb_pos)
        fclose(comb_pos);

    write_to_file(test_create_fullpath(sketch_opt_val->outdir, idx_sketch_suffix),
                  out_index, ((size_t)nfiles + 1) * sizeof(out_index[0]));
    const size_t stat_size = sizeof(dim_sketch_stat_t) + (size_t)nfiles * PATHLEN;
    write_to_file(test_create_fullpath(sketch_opt_val->outdir, sketch_stat),
                  in_sketch->mem_stat, stat_size);
    write_sketch_qc_stats(sketch_opt_val->outdir, out_qc, (size_t)nfiles);
    if (in_infile_meta)
        write_sketch_infile_meta_stats(sketch_opt_val->outdir, in_infile_meta, (size_t)nfiles);
    if (in_sketch->annotation)
        write_sketch_annotations(sketch_opt_val->outdir, in_sketch->annotation, (size_t)nfiles);

    if (!in_sketch->conflict) {
        warnx("%s(): input sketch has already removed conflicting context-objects; "
              "post-QC may not exactly reproduce sketch-time --readsQC conflict filtering.",
              __func__);
    }
    if (invalid_ranges) {
        warnx("%s(): %d/%d samples had no valid QC range and were copied unchanged",
              __func__, invalid_ranges, nfiles);
    }
    fprintf(stderr, "--sketchQC applied %s to %d samples: kept %lu/%lu sketch entries\n",
            sketch_qc_stat, nfiles, (unsigned long)total_kept, (unsigned long)total_in);

    free(out_index);
    free(out_qc);
    if (in_infile_meta)
        free_read_from_file(in_infile_meta, sizeof(in_infile_meta[0]) * (size_t)nfiles);
    free_unify_sketch(in_sketch);
    return 0;
}



KHASH_MAP_INIT_INT64(kmer_hash, int)
void mfa2sortedctxobj64(sketch_opt_t *sketch_opt_val, infile_tab_t *infile_stat)
{

    uint64_t tuple, crvstuple, unituple, basenum, unictx, totle_sketch_size = 0;
    uint32_t len_mv = 2 * klen - 2; // uint32_t saved_genome_num = infile_stat->infile_num, genome_num = 0;
    Vector sketch_index;
    vector_init(&sketch_index, sizeof(uint64_t));
    vector_push(&sketch_index, &totle_sketch_size);
    int tmpname_size_alloc = sketch_index.capacity;
    char (*tmpname)[PATHLEN] = malloc(tmpname_size_alloc * PATHLEN);
    if (!tmpname)
        err(errno, "%s():Memory allocation failed for tmpname", __func__);
    FILE *comb_sketch_fp;
    if ((comb_sketch_fp = fopen(format_string("%s/%s", sketch_opt_val->outdir, combined_sketch_suffix), "wb")) == NULL)
        err(errno, "%s() open file error: %s/%s", __func__, sketch_opt_val->outdir, combined_sketch_suffix);

    for (int i = 0; i < infile_stat->infile_num; i++)
    {
        char *seqfname = (infile_stat->organized_infile_tab)[i].fpath;
        gzFile infile = gzopen(seqfname, "r");
        if (!infile)
            err(errno, "reads2sketch64(): Cannot open file %s", seqfname);
        kseq_t *seq = kseq_init(infile);
        while (kseq_read(seq) >= 0 && seq->seq.l > klen)
        {
            // seq name handle
            if (sketch_index.size >= tmpname_size_alloc)
            {
                tmpname_size_alloc += 1000;
                tmpname = realloc(tmpname, tmpname_size_alloc * PATHLEN);
            }
            replace_special_chars_with_underscore(seq->name.s);
            strncpy(tmpname[sketch_index.size - 1], seq->name.s, PATHLEN);
            // sketching each fasta sequence
            khash_t(sort64) *h = kh_init(sort64);
            const char *s = seq->seq.s;
            int base = 0;
            uint32_t sketch_size = 0;

            for (int pos = 0; pos < seq->seq.l; pos++)
            {
                if (Basemap[(unsigned short)s[pos]] == DEFAULT)
                {
                    base = 0;
                    continue;
                }
                basenum = Basemap[(unsigned short)s[pos]];
                tuple = ((tuple << 2) | basenum);
                crvstuple = ((crvstuple >> 2) | ((basenum ^ 3LLU) << len_mv));
                if (++base < klen)
                    continue;

                unituple = (tuple & ctxmask) < (crvstuple & ctxmask) ? tuple : crvstuple;
                unictx = unituple & ctxmask;
                if (SKETCH_HASH(unictx) > FILTER)
                    continue;
                int ret;
                khint_t key = kh_put(sort64, h, uint64kmer2generic_ctxobj(unituple), &ret);

                if (ret)
                    kh_value(h, key) = 1;

            } // for line
            SortedKV_Arrays_t lco_ab = sort_khash_u64(h);
            totle_sketch_size += lco_ab.len;
            vector_push(&sketch_index, &totle_sketch_size);
            kh_destroy(sort64, h);
            fwrite(lco_ab.keys, sizeof(lco_ab.keys[0]), lco_ab.len, comb_sketch_fp);
            free_all(lco_ab.keys, lco_ab.values, NULL);
        } // while
        kseq_destroy(seq);
        gzclose(infile);
        printf("\r%dth/%d multifasta sketching %s completed!\t #genomes=%lu", i + 1,
               infile_stat->infile_num, (infile_stat->organized_infile_tab)[i].fpath, sketch_index.size - 1);
        if (i == infile_stat->infile_num - 1)
            printf("\n");
    } // for infile
    fclose(comb_sketch_fp);
    // write index	and stat file
    write_to_file(test_create_fullpath(sketch_opt_val->outdir, idx_sketch_suffix), sketch_index.data, sizeof(uint64_t) * sketch_index.size);
    comblco_stat_one.infile_num = sketch_index.size - 1;
    concat_and_write_to_file(test_create_fullpath(sketch_opt_val->outdir, sketch_stat),
                             &comblco_stat_one, sizeof(comblco_stat_one), tmpname, comblco_stat_one.infile_num * PATHLEN);

    vector_free(&sketch_index);
    free(tmpname);
} // mfa2sortedctxobj64()

void read_genomes2mem2sortedctxobj64(sketch_opt_t *sketch_opt_val, infile_tab_t *infile_stat, int batch_size)
{
    //   printf("%lx\n", ctxmask);
    uint32_t len_mv = 2 * klen - 2;
    uint64_t *sketch_index = calloc(infile_stat->infile_num + 1, sizeof(uint64_t));
    int *gseq_nums = calloc(batch_size + 1, sizeof(int));
    uint64_t **batch_sketches = malloc(batch_size * sizeof(uint64_t *));
    Vector all_reads;
    vector_init(&all_reads, sizeof(char *));

    FILE *comb_sketch_fp;
    if ((comb_sketch_fp = fopen(format_string("%s/%s", sketch_opt_val->outdir, combined_sketch_suffix), "wb")) == NULL)
        err(errno, "%s() open file error: %s/%s", __func__, sketch_opt_val->outdir, combined_sketch_suffix);
    for (int infile_num_p = 0; infile_num_p < infile_stat->infile_num; infile_num_p++)
    {
        gzFile infile = gzopen((infile_stat->organized_infile_tab)[infile_num_p].fpath, "r");
        if (!infile)
            err(errno, "%s(): Cannot open file %s", __func__, (infile_stat->organized_infile_tab)[infile_num_p].fpath);
        kseq_t *seq = kseq_init(infile);

        while (kseq_read(seq) >= 0)
        {
            char *read = malloc(seq->seq.l + 1);
            if (!read)
                err(errno, "%dth genome malloc failed", infile_num_p);
            memcpy(read, seq->seq.s, seq->seq.l);
            read[seq->seq.l] = '\0';
            char **add = &read;
            vector_push(&all_reads, add);
            gseq_nums[(infile_num_p % batch_size) + 1]++; // gseq_nums[infile_num_p+1]++;
        }
        kseq_destroy(seq);
        gzclose(infile);
        if ((infile_num_p < infile_stat->infile_num - 1) && infile_num_p % batch_size < batch_size - 1)
            continue;

        int num = infile_num_p % batch_size + 1;
        for (int i = 0; i < num; i++)
            gseq_nums[i + 1] += gseq_nums[i];
#pragma omp parallel for num_threads(sketch_opt_val->p)
        for (uint32_t i = 0; i < num; i++)
        {
            khash_t(sort64) *h = kh_init(sort64);
            for (uint32_t j = gseq_nums[i]; j < gseq_nums[i + 1]; j++)
            {
                char *s = *(char **)vector_get(&all_reads, j);
                int len = strlen(s);
                if (len < klen)
                    continue;
                int base = 0;
                uint64_t tuple, crvstuple, unituple, basenum, unictx;
                for (int pos = 0; pos < len; pos++)
                {
                    if (Basemap[(unsigned short)s[pos]] == DEFAULT)
                    {
                        base = 0;
                        continue;
                    }
                    basenum = Basemap[(unsigned short)s[pos]];
                    tuple = ((tuple << 2) | basenum);
                    crvstuple = ((crvstuple >> 2) | ((basenum ^ 3LLU) << len_mv));
                    if (++base < klen)
                        continue;

                    unituple = (tuple & ctxmask) < (crvstuple & ctxmask) ? tuple : crvstuple;
                    unictx = unituple & ctxmask;
                    // printf("%lx\t%lx\t%lx\t%lx\t%lx\t%lx\n", unituple, unictx, tuple,tuple & ctxmask, crvstuple, crvstuple & ctxmask);
                    if (SKETCH_HASH(unictx) > FILTER)
                        continue;

                    int ret;
                    khint_t key = kh_put(sort64, h, uint64kmer2generic_ctxobj(unituple & tupmask), &ret);
                    if (ret)
                        kh_value(h, key) = 1;
                } // nt pos loop
                free(s);
            } // seq j loop
            SortedKV_Arrays_t lco_ab = gpt_sort_khash_u64(h);
            // remove context with conflict object
            if (!sketch_opt_val->conflict)
                remove_ctx_with_conflict_obj(&lco_ab, Bitslen.obj);
            // may filter n for lco_ab
            batch_sketches[i] = lco_ab.keys;
            sketch_index[infile_num_p - num + i + 2] = lco_ab.len;
            kh_destroy(sort64, h);
            free(lco_ab.values);
        } // genome i loop

        for (uint32_t i = 0; i < num; i++)
        {
            fwrite(batch_sketches[i], sizeof(uint64_t), sketch_index[infile_num_p - num + i + 2], comb_sketch_fp);
            free(batch_sketches[i]);
        }

        memset(gseq_nums, 0, (batch_size + 1) * sizeof(int));
        vector_free(&all_reads);
        vector_init(&all_reads, sizeof(char *));

        printf("\r%d/%d genomes sketched\n", infile_num_p + 1, infile_stat->infile_num);

    } // infile_num_p loop

    for (int i = 0; i < infile_stat->infile_num; i++)
        sketch_index[i + 1] += sketch_index[i];
    write_to_file(format_string("%s/%s", sketch_opt_val->outdir, idx_sketch_suffix), sketch_index, (infile_stat->infile_num + 1) * sizeof(sketch_index[0]));

    fclose(comb_sketch_fp);
    vector_free(&all_reads);
    free_all(sketch_index, gseq_nums, batch_sketches, NULL);
    write_sketch_stat_ex(sketch_opt_val->outdir, infile_stat, sketch_opt_val->anno,
                         sketch_opt_val->compute_meta);
}

// for kssd3 ani use directly
simple_sketch_t *old_simple_genomes2mem2sortedctxobj64_mem(infile_tab_t *infile_stat, int drfold)
{ // only for coden pattern sketching*
    FILTER = UINT32_MAX >> drfold;
    uint32_t len_mv = 2 * klen - 2;

    uint64_t *sketch_index = calloc(infile_stat->infile_num + 1, sizeof(uint64_t));

    int *gseq_nums = calloc(infile_stat->infile_num + 1, sizeof(int));

    Vector all_reads;
    vector_init(&all_reads, sizeof(char *));

    int num = infile_stat->infile_num;
    uint64_t **batch_sketches = malloc(num * sizeof(uint64_t *));

    for (int infile_num_p = 0; infile_num_p < num; infile_num_p++)
    {
        gzFile infile = gzopen((infile_stat->organized_infile_tab)[infile_num_p].fpath, "r");
        if (!infile)
            err(errno, "%s(): Cannot open file %s", __func__, (infile_stat->organized_infile_tab)[infile_num_p].fpath);
        kseq_t *seq = kseq_init(infile);

        while (kseq_read(seq) >= 0)
        {
            char *read = malloc(seq->seq.l + 1);
            if (!read)
                err(errno, "%dth genome malloc failed", infile_num_p);
            memcpy(read, seq->seq.s, seq->seq.l);
            read[seq->seq.l] = '\0';
            char **add = &read;
            vector_push(&all_reads, add);
            gseq_nums[(infile_num_p) + 1]++; // gseq_nums[infile_num_p+1]++;
        }
        kseq_destroy(seq);
        gzclose(infile);

        if (infile_num_p < num - 1)
            continue;

        for (int i = 0; i < num; i++)
            gseq_nums[i + 1] += gseq_nums[i];

#pragma omp parallel for num_threads(2)
        for (uint32_t i = 0; i < num; i++)
        {

            khash_t(sort64) *h = kh_init(sort64);
            for (uint32_t j = gseq_nums[i]; j < gseq_nums[i + 1]; j++)
            {
                char *s = *(char **)vector_get(&all_reads, j);
                int len = strlen(s);
                if (len < klen)
                    continue;
                int base = 0;
                uint64_t tuple, crvstuple, unituple, basenum, unictx;
                for (int pos = 0; pos < len; pos++)
                {
                    if (Basemap[(unsigned short)s[pos]] == DEFAULT)
                    {
                        base = 0;
                        continue;
                    }
                    basenum = Basemap[(unsigned short)s[pos]];
                    tuple = ((tuple << 2) | basenum);
                    crvstuple = ((crvstuple >> 2) | ((basenum ^ 3LLU) << len_mv));
                    if (++base < klen)
                        continue;

                    unituple = (tuple & ctxmask) < (crvstuple & ctxmask) ? tuple : crvstuple;
                    unictx = unituple & ctxmask;
                    if (SKETCH_HASH(unictx) > FILTER)
                        continue;

                    int ret;
                    khint_t key = kh_put(sort64, h, reorder_unituple_by_coden_pattern64(unituple & tupmask), &ret);
                    if (ret)
                        kh_value(h, key) = 1;
                } // nt pos loop
                free(s);
            } // seq j loop

            SortedKV_Arrays_t lco_ab = sort_khash_u64(h);

            // remove context with conflict object
            remove_ctx_with_conflict_obj(&lco_ab, Bitslen.obj);

            // may filter n for lco_ab
            batch_sketches[i] = lco_ab.keys;
            sketch_index[i + 1] = lco_ab.len;

            kh_destroy(sort64, h);
            free(lco_ab.values);
        } // genome i loop

    } // infile_num_p loop
    for (int i = 0; i < num; i++)
        sketch_index[i + 1] += sketch_index[i];

    uint64_t *combined_sketch = malloc(sizeof(uint64_t) * sketch_index[num]);
    if (!combined_sketch)
        err(EXIT_FAILURE, "combined_sketch allocation failed");

    for (int i = 0; i < num; i++)
    {
        memcpy(combined_sketch + sketch_index[i], batch_sketches[i], sizeof(uint64_t) * (sketch_index[i + 1] - sketch_index[i]));
        free(batch_sketches[i]);
    }

    vector_free(&all_reads);
    free_all(gseq_nums, batch_sketches, NULL);

    simple_sketch_t *return_sketch = malloc(sizeof(simple_sketch_t));
    if (!return_sketch)
        err(EXIT_FAILURE, "return_sketch allocation failed");

    return_sketch->comb_sketch = combined_sketch;
    return_sketch->sketch_index = sketch_index;
    return_sketch->infile_num = num;

    return return_sketch;
}

// drop-in replacement for read_genomes2mem2sortedctxobj64()
// keeps identical interface and output semantics
void test_read_genomes2mem2sortedctxobj64(sketch_opt_t *sketch_opt_val, infile_tab_t *infile_stat, int batch_size)
{
    uint8_t nobjbits = Bitslen.obj;
    const uint32_t len_mv = (uint32_t)(2 * klen - 2); // shift for reverse rolling
    uint64_t *sketch_index = (uint64_t *)calloc((size_t)infile_stat->infile_num + 1, sizeof(uint64_t));
    if (!sketch_index)
        err(errno, "%s(): OOM sketch_index", __func__);

    // reuse across batches; only batch-local slots are touched
    uint64_t **batch_sketches = (uint64_t **)malloc((size_t)batch_size * sizeof(uint64_t *));
    if (!batch_sketches)
        err(errno, "%s(): OOM batch_sketches", __func__);

    // single combined output file, fully buffered
    FILE *comb_sketch_fp = fopen(format_string("%s/%s", sketch_opt_val->outdir, combined_sketch_suffix), "wb");
    if (!comb_sketch_fp)
        err(errno, "%s() open file error: %s/%s", __func__, sketch_opt_val->outdir, combined_sketch_suffix);
    // 8 MB buffered writes improve throughput on large outputs
    (void)setvbuf(comb_sketch_fp, NULL, _IOFBF, 8u << 20);
    //
    const int nfiles = infile_stat->infile_num;
    for (int batch_start = 0; batch_start < nfiles; batch_start += batch_size)
    {
        const int batch_end = (batch_start + batch_size <= nfiles) ? (batch_start + batch_size) : nfiles;
        const int this_batch = batch_end - batch_start;

        // lengths per genome in this batch
        uint64_t *lens = (uint64_t *)calloc((size_t)this_batch, sizeof(uint64_t));
        if (!lens)
            err(errno, "%s(): OOM lens", __func__);

// Process one genome per thread; no staging of reads in RAM.
#pragma omp parallel for num_threads(sketch_opt_val->p) schedule(dynamic, 1)
        for (int bi = 0; bi < this_batch; ++bi)
        {
            const int file_idx = batch_start + bi;
            const char *fpath = infile_stat->organized_infile_tab[file_idx].fpath;

            gzFile infile = gzopen(fpath, "r");
            if (!infile)
                err(errno, "%s(): Cannot open file %s", __func__, fpath);

            // enlarge zlib internal buffer to reduce syscalls
            (void)gzbuffer(infile, 1u << 20); // 1 MB

            kseq_t *seq = kseq_init(infile);
            if (!seq)
                err(errno, "%s(): kseq_init failed on %s", __func__, fpath);

            khash_t(sort64) *h = kh_init(sort64);
            if (!h)
                err(errno, "%s(): kh_init sort64 OOM", __func__);

            // Heuristic reserve: assume ~1/8 of bases produce valid kept k-mers after filtering.
            // Tweak if you know your sampling rate (c) to reduce rehash.
            kh_resize(sort64, h, 1u << 15); // start with 32K; grows as needed

            while (kseq_read(seq) >= 0)
            {
                const char *s = seq->seq.s;
                const int len = (int)seq->seq.l;
                if (len < klen)
                    continue;
                // rolling 2-bit canonical
                uint64_t tuple, crv; // forward and reverse-complement
                int base = 0;

                for (int pos = 0; pos < len; ++pos)
                {

                    const int bmap = Basemap[(unsigned char)s[pos]];
                    if (bmap == DEFAULT)
                    { // non-ACGT → reset window
                        base = 0;
                        continue;
                    }

                    const uint64_t b2 = (uint64_t)bmap;
                    tuple = (tuple << 2) | b2;
                    crv = (crv >> 2) | ((b2 ^ 3ull) << len_mv);
                    if (++base < klen)
                        continue;

                    // canonical by ctx, same as your original:
                    // compare (tuple & ctxmask) vs (crv & ctxmask)
                    const uint64_t t_ctx = (tuple & ctxmask);
                    const uint64_t r_ctx = (crv & ctxmask);
                    uint64_t unituple = (t_ctx < r_ctx) ? tuple : crv;
                    const uint64_t unictx = unituple & ctxmask;

                    // sketching decision based on context
                    if (SKETCH_HASH(unictx) > FILTER)
                        continue;
                    unituple &= tupmask;
                    // Rearrange to context-object (same length) and store encoded k-mer (not a hash)
                    // Your helper packs it; unchanged API:
                    // const uint64_t ctxobj = uint64kmer2generic_ctxobj(unituple & tupmask);
                    const uint64_t ctxobj = make_ctxobj(unituple, tupmask, ctxmask, nobjbits);

                    int ret;
                    khint_t key = kh_put(sort64, h, ctxobj, &ret);
                    if (ret)
                        kh_value(h, key) = 1; // presence; keep as in original
                } // end bases of this read
            } // end reads

            // finalize this genome: sort + (optionally) remove conflicts
            SortedKV_Arrays_t lco_ab = gpt_sort_khash_u64(h); // sort_khash_u64(h);
            if (!sketch_opt_val->conflict)
                remove_ctx_with_conflict_obj(&lco_ab, nobjbits);

            batch_sketches[bi] = lco_ab.keys; // will be written by main thread
            lens[bi] = lco_ab.len;

            kh_destroy(sort64, h);
            free(lco_ab.values);

            kseq_destroy(seq);
            gzclose(infile);
        } // end omp per-genome

        // Write this batch in file order, build sketch_index
        for (int bi = 0; bi < this_batch; ++bi)
        {
            const int global_i = batch_start + bi;
            const uint64_t n = lens[bi];

            // index prefix sums are built after all batches (like your original),
            // here we record per-genome lengths at sketch_index[i+1]
            sketch_index[global_i + 1] = n;

            if (n)
            {
                fwrite(batch_sketches[bi], sizeof(uint64_t), n, comb_sketch_fp);
                free(batch_sketches[bi]);
                batch_sketches[bi] = NULL;
            }
        }

        free(lens);

        printf("\r%d/%d genomes sketched\n", batch_end, nfiles);
        fflush(stdout);
    } // end batches

    // prefix sum → absolute offsets (same as your original)
    for (int i = 0; i < nfiles; ++i)
        sketch_index[i + 1] += sketch_index[i];

    write_to_file(format_string("%s/%s", sketch_opt_val->outdir, idx_sketch_suffix),
                  sketch_index,
                  (size_t)(nfiles + 1) * sizeof(sketch_index[0]));

    fclose(comb_sketch_fp);
    free(batch_sketches);
    free(sketch_index);

    write_sketch_stat_ex(sketch_opt_val->outdir, infile_stat, sketch_opt_val->anno,
                         sketch_opt_val->compute_meta);
}


// === unified helpers + both modes (A and B) ===
// ---- conflict filters provided by your code:
// void remove_ctx_with_conflict_obj(SortedKV_Arrays_t *kv, uint32_t n_obj_bits);
// void remove_ctx_with_conflict_obj_noabund(uint64_t *keys, size_t *len_io, uint32_t n_obj_bits);
#include <omp.h>
typedef struct{char *s; int l; uint64_t pos0;} read_span_t;
typedef struct{uint64_t *keys; uint64_t *positions; size_t n, cap;} u64posvec;

static inline void pv_init(u64posvec *v, size_t cap)
{
    v->keys = cap ? (uint64_t *)malloc(cap * sizeof(uint64_t)) : NULL;
    v->positions = cap ? (uint64_t *)malloc(cap * sizeof(uint64_t)) : NULL;
    if (cap && (!v->keys || !v->positions))
        err(errno, "%s(): OOM position vector", __func__);
    v->n = 0;
    v->cap = cap;
}

static inline void pv_free(u64posvec *v)
{
    free(v->keys);
    free(v->positions);
    v->keys = NULL;
    v->positions = NULL;
    v->n = v->cap = 0;
}

static inline void pv_reserve(u64posvec *v, size_t need)
{
    if (need <= v->cap)
        return;
    size_t nc = v->cap ? v->cap : 8192;
    while (nc < need)
        nc <<= 1;
    uint64_t *new_keys = (uint64_t *)realloc(v->keys, nc * sizeof(uint64_t));
    if (!new_keys)
        err(errno, "%s(): OOM position vector keys", __func__);
    v->keys = new_keys;
    uint64_t *new_positions = (uint64_t *)realloc(v->positions, nc * sizeof(uint64_t));
    if (!new_positions)
        err(errno, "%s(): OOM position vector positions", __func__);
    v->positions = new_positions;
    v->cap = nc;
}

static inline void pv_push(u64posvec *v, uint64_t key, uint64_t position)
{
    if (v->n == v->cap)
        pv_reserve(v, v->cap ? (v->cap << 1) : 8192);
    v->keys[v->n] = key;
    v->positions[v->n] = position;
    ++v->n;
}

static inline void u64vec_clear(u64vec *v) { v->n = 0; }

static inline void u64vec_destroy(u64vec *v) {
    free(v->a);
    v->a = NULL;
    v->n = v->cap = 0;
}

// ---- In-file parallel helpers (used only when #files < threads) ----

// Hot loop: push kept ctxobj into a vector (no hashing)
static inline void sketch_read_into_vec(const char *restrict s, int len, u64vec *restrict vec,
                                        uint64_t ctxmask, uint64_t tupmask, uint8_t nobjbits, uint32_t klen)
{
    if (len < (int)klen) return;
    const uint32_t len_mv = (uint32_t)(2 * klen - 2);

    uint64_t tuple = 0, crv = 0;
    int base = 0;

    for (int pos = 0; pos < len; ++pos)
    {
        const int bmap = Basemap[(unsigned char)s[pos]];
        if (unlikely(bmap == DEFAULT))
        {
            base = 0;
            tuple = 0;
            crv = 0;
            continue;
        }

        const uint64_t b2 = (uint64_t)bmap;
        tuple = (tuple << 2) | b2;
        crv = (crv >> 2) | ((b2 ^ 3ull) << len_mv);
        if (unlikely(++base < (int)klen)) continue;

        const uint64_t t_ctx = (tuple & ctxmask);
        const uint64_t r_ctx = (crv & ctxmask);
        const int use_fwd = (t_ctx < r_ctx);
        const uint64_t unictx = use_fwd ? t_ctx : r_ctx;

        if (unlikely(SKETCH_HASH(unictx) > FILTER)) continue;

        const uint64_t unituple = (use_fwd ? tuple : crv) & tupmask;
        const uint64_t ctxobj = make_ctxobj(unituple, tupmask, ctxmask, nobjbits);
        v_push(vec, ctxobj);
    }
}

static inline void sketch_read_into_posvec(const char *restrict s, int len, u64posvec *restrict vec,
                                           uint64_t seq_offset,
                                           uint64_t ctxmask, uint64_t tupmask,
                                           uint8_t nobjbits, uint32_t klen)
{
    if (len < (int)klen) return;
    const uint32_t len_mv = (uint32_t)(2 * klen - 2);

    uint64_t tuple = 0, crv = 0;
    int base = 0;

    for (int pos = 0; pos < len; ++pos)
    {
        const int bmap = Basemap[(unsigned char)s[pos]];
        if (unlikely(bmap == DEFAULT))
        {
            base = 0;
            tuple = 0;
            crv = 0;
            continue;
        }

        const uint64_t b2 = (uint64_t)bmap;
        tuple = (tuple << 2) | b2;
        crv = (crv >> 2) | ((b2 ^ 3ull) << len_mv);
        if (unlikely(++base < (int)klen)) continue;

        const uint64_t t_ctx = (tuple & ctxmask);
        const uint64_t r_ctx = (crv & ctxmask);
        const int use_fwd = (t_ctx < r_ctx);
        const uint64_t unictx = use_fwd ? t_ctx : r_ctx;

        if (unlikely(SKETCH_HASH(unictx) > FILTER)) continue;

        const uint64_t unituple = (use_fwd ? tuple : crv) & tupmask;
        const uint64_t ctxobj = make_ctxobj(unituple, tupmask, ctxmask, nobjbits);
        const uint64_t kmer_pos = seq_offset + (uint64_t)pos + 1u - (uint64_t)klen;
        pv_push(vec, ctxobj, kmer_pos);
    }
}

static inline void radix_sort_kp_u64(uint64_t *k, uint64_t *p, size_t n)
{
    if (n < 2) return;
    uint64_t *tk = (uint64_t *)malloc(n * sizeof(uint64_t));
    uint64_t *tp = (uint64_t *)malloc(n * sizeof(uint64_t));
    if (!tk || !tp)
        err(errno, "%s(): OOM radix key/position buffers", __func__);

    size_t cnt[256];
    for (unsigned pass = 0; pass < 8; ++pass) {
        for (int i = 0; i < 256; ++i) cnt[i] = 0;
        const unsigned sh = pass * 8;
        for (size_t i = 0; i < n; ++i) ++cnt[(k[i] >> sh) & 0xFFu];
        size_t sum = 0;
        for (int i = 0; i < 256; ++i) {
            const size_t c = cnt[i];
            cnt[i] = sum;
            sum += c;
        }
        for (size_t i = 0; i < n; ++i) {
            const unsigned b = (unsigned)((k[i] >> sh) & 0xFFu);
            const size_t dst = cnt[b]++;
            tk[dst] = k[i];
            tp[dst] = p[i];
        }
        uint64_t *swk = k; k = tk; tk = swk;
        uint64_t *swp = p; p = tp; tp = swp;
    }
    free(tk);
    free(tp);
}

static inline void radix_sort_kvp_u64(uint64_t *k, uint32_t *v, uint64_t *p, size_t n)
{
    if (n < 2) return;
    uint64_t *tk = (uint64_t *)malloc(n * sizeof(uint64_t));
    uint32_t *tv = (uint32_t *)malloc(n * sizeof(uint32_t));
    uint64_t *tp = (uint64_t *)malloc(n * sizeof(uint64_t));
    if (!tk || !tv || !tp)
        err(errno, "%s(): OOM radix key/count/position buffers", __func__);

    size_t cnt[256];
    for (unsigned pass = 0; pass < 8; ++pass) {
        for (int i = 0; i < 256; ++i) cnt[i] = 0;
        const unsigned sh = pass * 8;
        for (size_t i = 0; i < n; ++i) ++cnt[(k[i] >> sh) & 0xFFu];
        size_t sum = 0;
        for (int i = 0; i < 256; ++i) {
            const size_t c = cnt[i];
            cnt[i] = sum;
            sum += c;
        }
        for (size_t i = 0; i < n; ++i) {
            const unsigned b = (unsigned)((k[i] >> sh) & 0xFFu);
            const size_t dst = cnt[b]++;
            tk[dst] = k[i];
            tv[dst] = v[i];
            tp[dst] = p[i];
        }
        uint64_t *swk = k; k = tk; tk = swk;
        uint32_t *swv = v; v = tv; tv = swv;
        uint64_t *swp = p; p = tp; tp = swp;
    }
    free(tk);
    free(tv);
    free(tp);
}

static inline size_t dedup_with_positions(uint64_t *keys, uint64_t *positions, size_t n)
{
    if (n <= 1) return n;
    size_t j = 0;
    for (size_t i = 1; i < n; ++i) {
        if (keys[i] == keys[j]) {
            if (positions[i] < positions[j])
                positions[j] = positions[i];
        } else {
            ++j;
            keys[j] = keys[i];
            positions[j] = positions[i];
        }
    }
    return j + 1;
}

static inline size_t dedup_with_counts_and_positions(uint64_t *keys, uint64_t *positions,
                                                     size_t n, uint32_t **counts)
{
    if (n == 0) {
        *counts = NULL;
        return 0;
    }
    uint32_t *cnt = (uint32_t *)malloc(n * sizeof(uint32_t));
    if (!cnt)
        err(errno, "%s(): OOM count buffer", __func__);

    size_t j = 0;
    cnt[j] = 1;
    for (size_t i = 1; i < n; ++i) {
        if (keys[i] == keys[j]) {
            ++cnt[j];
            if (positions[i] < positions[j])
                positions[j] = positions[i];
        } else {
            ++j;
            keys[j] = keys[i];
            positions[j] = positions[i];
            cnt[j] = 1;
        }
    }
    uint32_t *tmp = (uint32_t *)realloc(cnt, (j + 1) * sizeof(uint32_t));
    *counts = tmp ? tmp : cnt;
    return j + 1;
}

static inline size_t shrink_kvp_inplace_u64_u32_pos(uint64_t *k, uint32_t *v, uint64_t *p, size_t n)
{
    if (!n) return 0;
    size_t w = 0;
    for (size_t r = 1; r < n; ++r) {
        if (k[r] == k[w]) {
            v[w] += v[r];
            if (p[r] < p[w])
                p[w] = p[r];
        } else {
            ++w;
            k[w] = k[r];
            v[w] = v[r];
            p[w] = p[r];
        }
    }
    return w + 1;
}

static inline void shrink_thread_vec(u64vec *vec, uint32_t **counts_out)
{
    // vec->a contains keys; produce counts_out aligned to deduped keys
    if (vec->n == 0) { *counts_out = NULL; return;}
    if (vec->n == 1){ // fast path
        *counts_out = (uint32_t *)malloc(sizeof(uint32_t));
        if (*counts_out) (*counts_out)[0] = 1;
        return;
    }
    radix_sort_u64(vec->a, vec->n);
    vec->n = dedup_with_counts(vec->a, vec->n, counts_out);
    if (!*counts_out)
    { // fallback if alloc failed inside helper
        *counts_out = (uint32_t *)malloc(vec->n * sizeof(uint32_t));
        if (*counts_out)
            for (size_t i = 0; i < vec->n; ++i)
                (*counts_out)[i] = 1u;
    }
}

static inline void shrink_thread_posvec(u64posvec *vec, uint32_t **counts_out)
{
    if (vec->n == 0) {
        *counts_out = NULL;
        return;
    }
    radix_sort_kp_u64(vec->keys, vec->positions, vec->n);
    vec->n = dedup_with_counts_and_positions(vec->keys, vec->positions, vec->n, counts_out);
}

static inline SortedKV_Arrays_t build_kv_from_vec(u64vec *vec, bool has_abundance)
{
    SortedKV_Arrays_t kv = (SortedKV_Arrays_t){0};
    if (vec->n == 0) return kv;

    radix_sort_u64(vec->a, vec->n);

    if (has_abundance) {
        uint32_t *counts = NULL;
        vec->n = dedup_with_counts(vec->a, vec->n, &counts);
        if (!counts) {
            counts = (uint32_t*)malloc(vec->n * sizeof(uint32_t));
            if (!counts) return kv; // OOM; kv stays empty, caller handles
            for (size_t i=0;i<vec->n;++i) counts[i] = 1u;
        }
        kv.keys   = (uint64_t*)malloc(vec->n * sizeof(uint64_t));
        kv.values = (uint32_t*)malloc(vec->n * sizeof(uint32_t));
        if (!kv.keys || !kv.values) { free(kv.keys); free(kv.values); free(counts); return (SortedKV_Arrays_t){0}; }
        memcpy(kv.keys, vec->a, vec->n * sizeof(uint64_t));
        memcpy(kv.values, counts, vec->n * sizeof(uint32_t));
        kv.len = vec->n;
        free(counts);
    } else {
        vec->n = dedup_sorted_uint64(vec->a, vec->n);
        kv.keys = (uint64_t*)malloc(vec->n * sizeof(uint64_t));
        if (!kv.keys) return (SortedKV_Arrays_t){0};
        memcpy(kv.keys, vec->a, vec->n * sizeof(uint64_t));
        kv.values = NULL;
        kv.len    = vec->n;
    }
    return kv;
}

static inline SortedKV_Arrays_t build_kv_from_posvec(u64posvec *vec, bool need_counts)
{
    SortedKV_Arrays_t kv = (SortedKV_Arrays_t){0};
    if (vec->n == 0) return kv;

    radix_sort_kp_u64(vec->keys, vec->positions, vec->n);

    uint32_t *counts = NULL;
    if (need_counts)
        vec->n = dedup_with_counts_and_positions(vec->keys, vec->positions, vec->n, &counts);
    else
        vec->n = dedup_with_positions(vec->keys, vec->positions, vec->n);

    kv.keys = (uint64_t *)malloc(vec->n * sizeof(uint64_t));
    kv.positions = (uint64_t *)malloc(vec->n * sizeof(uint64_t));
    if (need_counts)
        kv.values = (uint32_t *)malloc(vec->n * sizeof(uint32_t));
    if (!kv.keys || !kv.positions || (need_counts && !kv.values))
        err(errno, "%s(): OOM positional KV arrays", __func__);

    memcpy(kv.keys, vec->keys, vec->n * sizeof(uint64_t));
    memcpy(kv.positions, vec->positions, vec->n * sizeof(uint64_t));
    if (need_counts) {
        memcpy(kv.values, counts, vec->n * sizeof(uint32_t));
        free(counts);
    }
    kv.len = vec->n;
    return kv;
}

static inline void remove_ctx_with_conflict_obj_noabund_pos(SortedKV_Arrays_t *kv, uint32_t n_obj_bits)
{
    if (!kv || !kv->len)
        return;
    size_t write_idx = 0, i = 0;
    while (i < kv->len) {
        size_t count = 1;
        while (i + count < kv->len &&
               ((kv->keys[i] ^ kv->keys[i + count]) >> n_obj_bits) == 0) {
            ++count;
        }
        if (count == 1) {
            if (write_idx != i) {
                kv->keys[write_idx] = kv->keys[i];
                if (kv->positions)
                    kv->positions[write_idx] = kv->positions[i];
            }
            ++write_idx;
        }
        i += count;
    }
    kv->len = write_idx;
}

// ---------- unified FASTQ loader (gz or plain) → u64vec ----------
// Requires:
//   - sketch_read_into_vec(const char*, int, u64vec*, uint64_t, uint64_t, uint8_t, uint32_t)
//   - globals/params: ctxmask, tupmask, Bitslen.obj, klen
//   - zlib/kseq headers linked (for gz path)
// Build: add -fopenmp if you parallelize elsewhere; this function itself is single-threaded.

// small helper: check suffix (case-sensitive)
static inline int has_suffix(const char *s, const char *suf) {
    size_t ns=0, nf=0; while (s[ns]) ++ns; while (suf[nf]) ++nf;
    return (nf<=ns) && (memcmp(s+ns-nf, suf, nf)==0);
}

// =====================================================
// Mode-A helpers: single-vector collectors (gz / plain)
// =====================================================

// Single-vector collector for gz FASTQ
static inline void load_fastx_into_single_vec(
    char *path, const char *pipecmd, u64vec *vec,
    uint64_t ctxmask, uint64_t tupmask, uint8_t nobjbits, uint32_t klen,
    uint64_t **asm_lengths, size_t *asm_n, size_t *asm_cap)
{
    sketch_stream_t stream = open_sketch_stream(path, pipecmd);
    gzFile in = stream.gz;
    (void)gzbuffer(in, 4u << 20);

    kseq_t *seq = kseq_init(in);
    if (!seq) err(errno, "%s(): kseq_init %s", __func__, path);

    while (kseq_read(seq) >= 0) {
        if (asm_lengths && asm_n && asm_cap)
            append_asm_length(asm_lengths, asm_n, asm_cap, (uint64_t)seq->seq.l);
        if (seq->seq.l) {
            sketch_read_into_vec(seq->seq.s, (int)seq->seq.l, vec,
                                 ctxmask, tupmask, nobjbits, klen);
        }
    }
    kseq_destroy(seq);
    close_sketch_stream(&stream);
}

static inline uint64_t load_fastx_into_single_posvec(
    char *path, const char *pipecmd, u64posvec *vec, uint64_t start_offset,
    uint64_t ctxmask, uint64_t tupmask, uint8_t nobjbits, uint32_t klen,
    uint64_t **asm_lengths, size_t *asm_n, size_t *asm_cap)
{
    sketch_stream_t stream = open_sketch_stream(path, pipecmd);
    gzFile in = stream.gz;
    (void)gzbuffer(in, 4u << 20);

    kseq_t *seq = kseq_init(in);
    if (!seq) err(errno, "%s(): kseq_init %s", __func__, path);

    uint64_t seq_offset = start_offset;
    while (kseq_read(seq) >= 0) {
        if (asm_lengths && asm_n && asm_cap)
            append_asm_length(asm_lengths, asm_n, asm_cap, (uint64_t)seq->seq.l);
        if (seq->seq.l) {
            sketch_read_into_posvec(seq->seq.s, (int)seq->seq.l, vec, seq_offset,
                                    ctxmask, tupmask, nobjbits, klen);
        }
        seq_offset += (uint64_t)seq->seq.l;
    }
    kseq_destroy(seq);
    close_sketch_stream(&stream);
    return seq_offset;
}

// Single-vector collector for plain FASTQ via mmap (no copies)
static inline void load_fastq_plain_mmap_into_single_vec(char *path, u64vec *vec, uint64_t ctxmask, uint64_t tupmask, uint8_t nobjbits, uint32_t klen,
                                                         uint64_t **asm_lengths, size_t *asm_n, size_t *asm_cap)
{
    int fd = open(path, O_RDONLY);
    if (fd < 0) err(errno, "%s(): open %s", __func__, path);
#ifdef POSIX_FADV_SEQUENTIAL
    posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
    struct stat st;
    if (fstat(fd, &st) != 0) err(errno, "%s(): fstat %s", __func__, path);
    const size_t fsz = (size_t)st.st_size;
    if (fsz == 0) { close(fd); return; }

    char *base = (char*)mmap(NULL, fsz, PROT_READ, MAP_PRIVATE, fd, 0);
    if (base == MAP_FAILED) err(errno, "%s(): mmap %s", __func__, path);
#ifdef MADV_SEQUENTIAL
    madvise(base, fsz, MADV_SEQUENTIAL);
#endif
#ifdef MADV_WILLNEED
    madvise(base, fsz, MADV_WILLNEED);
#endif

    // Walk lines; every 2nd of each 4 is the sequence
    size_t line_start = 0;
    unsigned line_mod = 0; // 0=@hdr,1=seq,2=+,3=qual
    for (size_t i = 0; i < fsz; ++i) {
        if (base[i] != '\n') continue;
        size_t L = i - line_start;
        if (L && base[i-1] == '\r') --L; // trim CR
        if (line_mod == 1 && L) {
            if (asm_lengths && asm_n && asm_cap)
                append_asm_length(asm_lengths, asm_n, asm_cap, (uint64_t)L);
            sketch_read_into_vec(base + line_start, (int)L, vec,
                                 ctxmask, tupmask, nobjbits, klen);
        }
        line_start = i + 1;
        line_mod = (line_mod + 1) & 3;
    }
    // (Optional) handle last unterminated line — FASTQ typically ends with '\n'

    munmap(base, fsz);
    close(fd);
}

// Tiny wrapper to choose gz vs plain
static inline void load_genome_into_single_vec(
    char *path, const sketch_opt_t *opt, u64vec *vec,
    uint64_t ctxmask, uint64_t tupmask, uint8_t nobjbits, uint32_t klen,
    uint64_t **asm_lengths, size_t *asm_n, size_t *asm_cap)
{
    if (uses_stream_input(opt, path) || isCompressfile(path) || isOK_fmt_infile(path, fasta_fmt, FAS_FMT_SZ))
        load_fastx_into_single_vec(path, opt ? opt->pipecmd : NULL, vec, ctxmask, tupmask, nobjbits, klen,
                                   asm_lengths, asm_n, asm_cap);
    else if (isOK_fmt_infile(path, fastq_fmt, FQ_FMT_SZ))
        load_fastq_plain_mmap_into_single_vec(path, vec, ctxmask, tupmask, nobjbits, klen,
                                              asm_lengths, asm_n, asm_cap);
    else  err(errno, "%s(): %s is not accept format(.fasta,.fastq)",__func__, path);    
}

static inline uint64_t load_genome_into_single_posvec(
    char *path, const sketch_opt_t *opt, u64posvec *vec, uint64_t start_offset,
    uint64_t ctxmask, uint64_t tupmask, uint8_t nobjbits, uint32_t klen,
    uint64_t **asm_lengths, size_t *asm_n, size_t *asm_cap)
{
    return load_fastx_into_single_posvec(path, opt ? opt->pipecmd : NULL, vec, start_offset,
                                         ctxmask, tupmask, nobjbits, klen,
                                         asm_lengths, asm_n, asm_cap);
}

// =====================================================
// ---- Mode A: per-genome parallel (vector + counts) --
// =====================================================
static void sketch_many_files_in_parallel(sketch_opt_t *opt, infile_tab_t *tab, int batch_size)
{
    const int nfiles = tab->infile_num;
    const bool has_abundance = opt->abundance ;
    dim_sketch_qc_stat_t *sample_qc_stats = NULL;
    if (sketch_records_sample_qc(opt)) {
        sample_qc_stats = (dim_sketch_qc_stat_t *)calloc((size_t)nfiles, sizeof(sample_qc_stats[0]));
        if (!sample_qc_stats) err(errno, "%s(): OOM sample_qc_stats", __func__);
    }

    FILE *comb = fopen(format_string("%s/%s", opt->outdir, combined_sketch_suffix), "wb");
    if (!comb) err(errno, "%s() open file error: %s/%s", __func__, opt->outdir, combined_sketch_suffix);
    setvbuf(comb, NULL, _IOFBF, 8u << 20);

    FILE *comb_ab = NULL;
    if (has_abundance) {
        comb_ab = fopen(format_string("%s/%s", opt->outdir, combined_ab_suffix), "wb");
        if (!comb_ab) err(errno, "%s() open file error: %s/%s", __func__, opt->outdir, combined_ab_suffix);
        setvbuf(comb_ab, NULL, _IOFBF, 8u << 20);
    }
    FILE *comb_pos = NULL;
    if (opt->position) {
        comb_pos = fopen(format_string("%s/%s", opt->outdir, sketch_position_suffix), "wb");
        if (!comb_pos) err(errno, "%s() open file error: %s/%s", __func__, opt->outdir, sketch_position_suffix);
        setvbuf(comb_pos, NULL, _IOFBF, 8u << 20);
    }

    uint64_t *sketch_index = (uint64_t *)calloc((size_t)nfiles + 1, sizeof(uint64_t));
    if (!sketch_index) err(errno, "%s(): OOM index", __func__);
    infile_meta_t *sample_infile_meta = NULL;
    if (opt->compute_meta) {
        sample_infile_meta =
            (infile_meta_t *)calloc((size_t)nfiles, sizeof(sample_infile_meta[0]));
        if (!sample_infile_meta) err(errno, "%s(): OOM sample_infile_meta", __func__);
    }

    // Process in batches of files to cap memory if needed
    for (int batch_start = 0; batch_start < nfiles; batch_start += batch_size) {
        const int batch_end  = (batch_start + batch_size <= nfiles) ? (batch_start + batch_size) : nfiles;
        const int this_batch = batch_end - batch_start;

        uint64_t *lens = (uint64_t *)calloc((size_t)this_batch, sizeof(uint64_t));
        if (!lens) err(errno, "%s(): OOM lens", __func__);

        uint64_t **batch_keys = (uint64_t **)calloc((size_t)this_batch, sizeof(uint64_t *));
        if (!batch_keys) err(errno, "%s(): OOM batch_keys", __func__);

        uint64_t **batch_pos = NULL;
        if (opt->position) {
            batch_pos = (uint64_t **)calloc((size_t)this_batch, sizeof(uint64_t *));
            if (!batch_pos) err(errno, "%s(): OOM batch_pos", __func__);
        }

        kmer_count_filter_result_t *batch_count_cutoffs = NULL;
        if (sketch_reports_count_filter(opt) || sketch_records_sample_qc(opt)) {
            batch_count_cutoffs = (kmer_count_filter_result_t *)calloc((size_t)this_batch, sizeof(kmer_count_filter_result_t));
            if (!batch_count_cutoffs) err(errno, "%s(): OOM batch_count_cutoffs", __func__);
        }

        uint32_t **batch_vals = NULL;
        if (has_abundance) {
            batch_vals = (uint32_t **)calloc((size_t)this_batch, sizeof(uint32_t *));
            if (!batch_vals) err(errno, "%s(): OOM batch_vals", __func__);
        }

        // One thread per file (parallel across files)
        #pragma omp parallel for num_threads(opt->p) schedule(dynamic, 1)
        for (int bi = 0; bi < this_batch; ++bi) {
            const int file_idx = batch_start + bi;
            char *fpath  = tab->organized_infile_tab[file_idx].fpath;

            uint64_t *asm_lengths = NULL;
            size_t asm_n = 0, asm_cap = 0;
            const bool collect_meta_lengths = should_collect_lengths_for_meta(opt, fpath);
            SortedKV_Arrays_t kv = (SortedKV_Arrays_t){0};
            if (opt->position) {
                u64posvec pvec;
                pv_init(&pvec, 1u << 15);
                load_genome_into_single_posvec(fpath, opt, &pvec, 0, ctxmask, tupmask, Bitslen.obj, klen,
                                               collect_meta_lengths ? &asm_lengths : NULL,
                                               collect_meta_lengths ? &asm_n : NULL,
                                               collect_meta_lengths ? &asm_cap : NULL);
                if (pvec.n)
                    kv = build_kv_from_posvec(&pvec, sketch_needs_kmer_counts(opt));
                pv_free(&pvec);
            } else {
                u64vec vec; v_init(&vec, 1u << 15);
                load_genome_into_single_vec(fpath, opt, &vec, ctxmask, tupmask, Bitslen.obj, klen,
                                            collect_meta_lengths ? &asm_lengths : NULL,
                                            collect_meta_lengths ? &asm_n : NULL,
                                            collect_meta_lengths ? &asm_cap : NULL);
                if (vec.n)
                    kv = build_kv_from_vec(&vec, sketch_needs_kmer_counts(opt));
                v_free(&vec);
            }
            if (opt->compute_meta)
                infile_meta_from_lengths(asm_lengths, asm_n, infile_fmt_from_path(fpath),
                                         create_type_from_opt(opt),
                                         infile_flags_from_path(fpath, opt),
                                         &sample_infile_meta[file_idx]);
            free(asm_lengths);

            if (kv.len) {
                const kmer_count_filter_result_t count_cutoff = apply_kmer_count_filters(&kv, opt);
                if (batch_count_cutoffs) batch_count_cutoffs[bi] = count_cutoff;
                if (!opt->conflict){
                    if(has_abundance || kv.values) remove_ctx_with_conflict_obj(&kv, Bitslen.obj);
                    else if (kv.positions) remove_ctx_with_conflict_obj_noabund_pos(&kv, Bitslen.obj);
                    else remove_ctx_with_conflict_obj_noabund(kv.keys, &(kv.len), Bitslen.obj);
                }
                if (!kv.len) {
                    free(kv.keys);
                    free(kv.values);
                    free(kv.positions);
                    lens[bi] = 0;
                    continue;
                }
                batch_keys[bi] = kv.keys;
                lens[bi]       = kv.len;
                if (opt->position) batch_pos[bi] = kv.positions;
                if (has_abundance) batch_vals[bi] = kv.values;
                else               free(kv.values);
            } else {
                lens[bi] = 0;
                free(kv.keys);
                free(kv.values);
                free(kv.positions);
            }
        }

        // write batch in order
        for (int bi = 0; bi < this_batch; ++bi) {
            const int gi = batch_start + bi;
            sketch_index[gi + 1] = lens[bi];
            if (sample_qc_stats && batch_count_cutoffs)
                sample_qc_stats[gi] = batch_count_cutoffs[bi].sample_qc;
            if (batch_count_cutoffs && sketch_reports_count_filter(opt)) {
                print_kmer_count_filter_result(opt, tab->organized_infile_tab[gi].fpath,
                                               &batch_count_cutoffs[bi], true);
            }
            if (lens[bi]) {
                fwrite(batch_keys[bi], sizeof(uint64_t), lens[bi], comb);
                free(batch_keys[bi]);
                if (opt->position) {
                    fwrite(batch_pos[bi], sizeof(uint64_t), lens[bi], comb_pos);
                    free(batch_pos[bi]);
                }
                if (has_abundance) {
                    fwrite(batch_vals[bi], sizeof(uint32_t), lens[bi], comb_ab);
                    free(batch_vals[bi]);
                }
            }
        }

        free(batch_keys);
        free(batch_pos);
        free(lens);
        free(batch_count_cutoffs);
        if (has_abundance) free(batch_vals);

        fprintf(stderr, "\r%d/%d genomes sketched", batch_end, nfiles);
        fflush(stderr);
    }
    fprintf(stderr, "\n");

    // prefix-sum index and finish
    for (int i = 0; i < nfiles; ++i) sketch_index[i + 1] += sketch_index[i];
    write_to_file(format_string("%s/%s", opt->outdir, idx_sketch_suffix),
                  sketch_index, (size_t)(nfiles + 1) * sizeof(uint64_t));
    fclose(comb);
    if (has_abundance) fclose(comb_ab);
    if (opt->position) fclose(comb_pos);
    free(sketch_index);
    write_sketch_stat_ex(opt->outdir, tab, opt->anno && !sketch_has_stream_input(opt, tab), false);
    if (opt->compute_meta)
        write_sketch_infile_meta_stats(opt->outdir, sample_infile_meta, (size_t)nfiles);
    write_sketch_qc_stats(opt->outdir, sample_qc_stats, (size_t)nfiles);
    free(sample_infile_meta);
    free(sample_qc_stats);
}

// ============================
// Unified helpers + Mode-B
// ============================


// C helper for finding '\n' between [off, lim)
static inline char *find_nl(const char *base, size_t off, size_t lim) {
    if (off >= lim) return NULL;
    return (char*)memchr(base + off, '\n', lim - off);
}

// ---------- shrink helpers ----------
static inline unsigned clamp_bucket_bits_from_total(size_t total_distinct){
    // target ~4k entries per bucket; clamp to [8..14]
    size_t target_buckets = (total_distinct / 4096) ? (total_distinct / 4096) : 256;
    unsigned BITS = 8;
    while (((size_t)1 << BITS) < target_buckets && BITS < 14) ++BITS;
    return BITS;
}

static inline size_t shrink_kv_inplace_u64_u32(uint64_t *k, uint32_t *v, size_t n){
    if (!n) return 0;
    size_t w=0;
    for (size_t r=1; r<n; ++r){
        if (k[r]==k[w]) v[w] += v[r];
        else { ++w; k[w]=k[r]; v[w]=v[r]; }
    }
    return w+1;
}

// ---------- collectors (no finalization here) ----------
//nth = how many independent k-mer buffers you have
//nthreads_for_omp = how many OpenMP threads may run concurrently
static void collect_gz_into_vectors(const char *path, const char *pipecmd, int nth, int nthreads_for_omp,
                                    u64vec *V_thr,
                                    uint64_t ctxmask, uint64_t tupmask, uint8_t nobjbits, uint32_t klen,
                                    int BATCH_READS,
                                    uint64_t **asm_lengths, size_t *asm_n, size_t *asm_cap)
{
    sketch_stream_t stream = open_sketch_stream(path, pipecmd);
    gzFile in = stream.gz;
    (void)gzbuffer(in, 4u<<20);
    kseq_t *seq = kseq_init(in);
    if (!seq) err(errno, "%s(): kseq_init %s", __func__, path);

    read_span_t *R = (read_span_t*)malloc((size_t)BATCH_READS*sizeof(*R));
    if (!R) err(errno, "%s(): OOM batch reads", __func__);
    int rcnt=0;

    int ret;
    do{
        ret = kseq_read(seq);
        if (ret >= 0){
            if (asm_lengths && asm_n && asm_cap)
                append_asm_length(asm_lengths, asm_n, asm_cap, (uint64_t)seq->seq.l);
            char *buf = (char*)malloc(seq->seq.l + 1);
            if (!buf) err(errno, "%s(): OOM read buf", __func__);
            memcpy(buf, seq->seq.s, seq->seq.l);
            buf[seq->seq.l] = '\0';
            R[rcnt].s = buf; R[rcnt].l = (int)seq->seq.l;
            ++rcnt;
        }
        if (rcnt >= BATCH_READS || ret < 0){
            #pragma omp parallel for if (nth>1) num_threads(nthreads_for_omp) schedule(static)
            for (int ri=0; ri<rcnt; ++ri){
                int tid = 0;
                #ifdef _OPENMP
                tid = omp_get_thread_num();
                #endif
                sketch_read_into_vec(R[ri].s, R[ri].l, &V_thr[tid],
                                     ctxmask, tupmask, nobjbits, klen);
            }
            for (int i=0;i<rcnt;++i) free(R[i].s);
            rcnt=0;
        }
    } while (ret >= 0);
    free(R);
    kseq_destroy(seq);
    close_sketch_stream(&stream);
}

static uint64_t collect_gz_into_posvectors(const char *path, const char *pipecmd,
                                           int nth, int nthreads_for_omp,
                                           u64posvec *V_thr, uint64_t start_offset,
                                           uint64_t ctxmask, uint64_t tupmask,
                                           uint8_t nobjbits, uint32_t klen,
                                           int BATCH_READS,
                                           uint64_t **asm_lengths, size_t *asm_n, size_t *asm_cap)
{
    sketch_stream_t stream = open_sketch_stream(path, pipecmd);
    gzFile in = stream.gz;
    (void)gzbuffer(in, 4u<<20);
    kseq_t *seq = kseq_init(in);
    if (!seq) err(errno, "%s(): kseq_init %s", __func__, path);

    read_span_t *R = (read_span_t*)malloc((size_t)BATCH_READS*sizeof(*R));
    if (!R) err(errno, "%s(): OOM batch reads", __func__);
    int rcnt = 0;
    uint64_t seq_offset = start_offset;

    int ret;
    do {
        ret = kseq_read(seq);
        if (ret >= 0) {
            if (asm_lengths && asm_n && asm_cap)
                append_asm_length(asm_lengths, asm_n, asm_cap, (uint64_t)seq->seq.l);
            char *buf = (char*)malloc(seq->seq.l + 1);
            if (!buf) err(errno, "%s(): OOM read buf", __func__);
            memcpy(buf, seq->seq.s, seq->seq.l);
            buf[seq->seq.l] = '\0';
            R[rcnt].s = buf;
            R[rcnt].l = (int)seq->seq.l;
            R[rcnt].pos0 = seq_offset;
            seq_offset += (uint64_t)seq->seq.l;
            ++rcnt;
        }
        if (rcnt >= BATCH_READS || ret < 0) {
            #pragma omp parallel for if (nth>1) num_threads(nthreads_for_omp) schedule(static)
            for (int ri=0; ri<rcnt; ++ri) {
                int tid = 0;
                #ifdef _OPENMP
                tid = omp_get_thread_num();
                #endif
                sketch_read_into_posvec(R[ri].s, R[ri].l, &V_thr[tid], R[ri].pos0,
                                        ctxmask, tupmask, nobjbits, klen);
            }
            for (int i=0; i<rcnt; ++i) free(R[i].s);
            rcnt = 0;
        }
    } while (ret >= 0);
    free(R);
    kseq_destroy(seq);
    close_sketch_stream(&stream);
    return seq_offset;
}

static void collect_plain_mmap_into_vectors(const char *path, int nth, int nthreads_for_omp,
                                            u64vec *V_thr,
                                            uint64_t ctxmask, uint64_t tupmask, uint8_t nobjbits, uint32_t klen,
                                            uint64_t **asm_lengths, size_t *asm_n, size_t *asm_cap)
{
    int fd = open(path, O_RDONLY);
    if (fd < 0) err(errno, "%s(): open %s", __func__, path);
#ifdef POSIX_FADV_SEQUENTIAL
    posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
    struct stat st;
    if (fstat(fd, &st) != 0) err(errno, "%s(): fstat %s", __func__, path);
    const size_t fsz = (size_t)st.st_size;
    if (fsz == 0) { close(fd); return; }

    char *base = (char*)mmap(NULL, fsz, PROT_READ, MAP_PRIVATE, fd, 0);
    if (base == MAP_FAILED) err(errno, "%s(): mmap %s", __func__, path);
#ifdef MADV_SEQUENTIAL
    madvise(base, fsz, MADV_SEQUENTIAL);
#endif
#ifdef MADV_WILLNEED
    madvise(base, fsz, MADV_WILLNEED);
#endif

    const size_t chunk = (fsz + (size_t)nth - 1) / (size_t)nth;

    #pragma omp parallel for if (nth>1) num_threads(nthreads_for_omp) schedule(static)
    for (int tid=0; tid<nth; ++tid){
        size_t start = (size_t)tid * chunk;
        size_t end   = (tid == nth-1) ? fsz : start + chunk;
        if (start >= fsz) { 
            // V_thr[tid].n = 0;  caller controls clearing; do NOT reset in collector -- for asone mode considerarion
            continue;
        }
        if (start > 0) { while (start < end && base[start-1] != '\n') ++start; }

        size_t cursor = start;
        while (cursor < end) {
            if (base[cursor] == '@') {
                char *nl1 = find_nl(base, cursor, fsz);                 if (!nl1) break;
                size_t seq_beg = (size_t)(nl1 - base) + 1;
                char *nl2 = find_nl(base, seq_beg, fsz);                if (!nl2) break;
                size_t plus_beg = (size_t)(nl2 - base) + 1;
                char *nl3 = find_nl(base, plus_beg, fsz);               if (!nl3) break;
                size_t qual_beg = (size_t)(nl3 - base) + 1;
                char *nl4 = find_nl(base, qual_beg, fsz);               if (!nl4) break;

                if (cursor >= start) {
                    size_t L = (size_t)(nl2 - (base + seq_beg));
                    if (L && base[seq_beg + L - 1] == '\r') --L;
                    if (L) {
                        if (asm_lengths && asm_n && asm_cap) {
                            #pragma omp critical (append_mmap_meta_length)
                            append_asm_length(asm_lengths, asm_n, asm_cap, (uint64_t)L);
                        }
                        sketch_read_into_vec(base + seq_beg, (int)L, &V_thr[tid],
                                             ctxmask, tupmask, nobjbits, klen);
                    }
                }
                cursor = (size_t)(nl4 - base) + 1;
                if (cursor >= end && tid != nth-1) break;
            } else {
                char *nl = find_nl(base, cursor, end);
                if (!nl) break;
                cursor = (size_t)(nl - base) + 1;
            }
        }
    }

    munmap(base, fsz);
    close(fd);
}

// ---------- shared finalizer (shrink + bucketed merge) ----------
static SortedKV_Arrays_t finalize_vectors_bucketed(u64vec *V_thr, int nth, int nthreads_for_omp)
{
    SortedKV_Arrays_t out = (SortedKV_Arrays_t){0};

    uint32_t **C_thr = (uint32_t**)calloc((size_t)nth, sizeof(uint32_t*));
    if (!C_thr) return out;

    size_t pre_elems = 0;
    for (int t = 0; t < nth; ++t) pre_elems += V_thr[t].n;

    #pragma omp parallel for if (nth>1 && pre_elems >= (1u<<20)) num_threads(nthreads_for_omp) schedule(dynamic)
    for (int t = 0; t < nth; ++t)
        shrink_thread_vec(&V_thr[t], &C_thr[t]);

    size_t total_distinct = 0;
    for (int t = 0; t < nth; ++t) total_distinct += V_thr[t].n;

    if (total_distinct == 0) {
        for (int t = 0; t < nth; ++t) {
            v_free(&V_thr[t]);
            free(C_thr[t]);
        }
        free(C_thr);
        return out;
    }

    const unsigned BITS = clamp_bucket_bits_from_total(total_distinct);
    const size_t   NB   = (size_t)1u << BITS;

    size_t *bucket_totals = (size_t*)calloc(NB, sizeof(size_t));
    size_t **thr_counts   = (size_t**)malloc((size_t)nth * sizeof(size_t*));
    if (!bucket_totals || !thr_counts) goto CLEAN_EMPTY;

    for (int t = 0; t < nth; ++t) {
        thr_counts[t] = (size_t*)calloc(NB, sizeof(size_t));
        if (!thr_counts[t]) goto CLEAN_EMPTY;
    }

    #pragma omp parallel for num_threads(nthreads_for_omp) schedule(static)
    for (int t = 0; t < nth; ++t) {
        uint64_t *A = V_thr[t].a;
        size_t n = V_thr[t].n;
        for (size_t i = 0; i < n; ++i) {
            unsigned b = (unsigned)(A[i] >> (64 - BITS));
            ++thr_counts[t][b];
        }
    }

    for (size_t b = 0; b < NB; ++b) {
        size_t sum = 0;
        for (int t = 0; t < nth; ++t) sum += thr_counts[t][b];
        bucket_totals[b] = sum;
    }

    size_t *bucket_off = (size_t*)malloc((NB + 1) * sizeof(size_t));
    if (!bucket_off) goto CLEAN_EMPTY;

    bucket_off[0] = 0;
    for (size_t b = 0; b < NB; ++b)
        bucket_off[b + 1] = bucket_off[b] + bucket_totals[b];

    const size_t TOTAL = bucket_off[NB];

    uint64_t *B_keys = (uint64_t*)malloc(TOTAL * sizeof(uint64_t));
    uint32_t *B_vals = (uint32_t*)malloc(TOTAL * sizeof(uint32_t));
    if (!B_keys || !B_vals) goto CLEAN_EMPTY;

    size_t **thr_write = (size_t**)malloc((size_t)nth * sizeof(size_t*));
    if (!thr_write) goto CLEAN_EMPTY;

    for (int t = 0; t < nth; ++t) {
        thr_write[t] = (size_t*)malloc(NB * sizeof(size_t));
        if (!thr_write[t]) goto CLEAN_EMPTY;
    }

    for (size_t b = 0; b < NB; ++b) {
        size_t off = bucket_off[b];
        for (int t = 0; t < nth; ++t) {
            thr_write[t][b] = off;
            off += thr_counts[t][b];
        }
    }

    #pragma omp parallel for num_threads(nthreads_for_omp) schedule(static)
    for (int t = 0; t < nth; ++t) {
        uint64_t *Ak = V_thr[t].a;
        uint32_t *Av = C_thr[t];
        size_t n = V_thr[t].n;
        size_t *pos = thr_write[t];
        for (size_t i = 0; i < n; ++i) {
            uint64_t k = Ak[i];
            unsigned b = (unsigned)(k >> (64 - BITS));
            size_t p = pos[b]++;
            B_keys[p] = k;
            B_vals[p] = Av ? Av[i] : 1u;
        }
    }

    for (int t = 0; t < nth; ++t) {
        v_free(&V_thr[t]);          // free per-thread buffer only
        free(C_thr[t]);
        free(thr_counts[t]);
        free(thr_write[t]);
    }
    free(C_thr);
    free(thr_counts);
    free(thr_write);

    uint64_t **bk_keys = (uint64_t**)malloc(NB * sizeof(uint64_t*));
    uint32_t **bk_vals = (uint32_t**)malloc(NB * sizeof(uint32_t*));
    size_t    *bk_len  = (size_t*)   malloc(NB * sizeof(size_t));
    if (!bk_keys || !bk_vals || !bk_len) goto CLEAN_BUCKET;

    #pragma omp parallel for num_threads(nthreads_for_omp) schedule(dynamic)
    for (size_t b = 0; b < NB; ++b) {
        size_t begin = bucket_off[b];
        size_t end   = bucket_off[b + 1];
        size_t n = end - begin;
        if (!n) {
            bk_keys[b] = NULL;
            bk_vals[b] = NULL;
            bk_len[b]  = 0;
            continue;
        }
        uint64_t *kseg = B_keys + begin;
        uint32_t *vseg = B_vals + begin;
        radix_sort_kv_u64(kseg, vseg, n);
        size_t m = shrink_kv_inplace_u64_u32(kseg, vseg, n);
        bk_keys[b] = (uint64_t*)malloc(m * sizeof(uint64_t));
        bk_vals[b] = (uint32_t*)malloc(m * sizeof(uint32_t));
        memcpy(bk_keys[b], kseg, m * sizeof(uint64_t));
        memcpy(bk_vals[b], vseg, m * sizeof(uint32_t));
        bk_len[b] = m;
    }

    free(B_keys);
    free(B_vals);
    free(bucket_totals);

    size_t *out_off = (size_t*)malloc((NB + 1) * sizeof(size_t));
    out_off[0] = 0;
    for (size_t b = 0; b < NB; ++b)
        out_off[b + 1] = out_off[b] + bk_len[b];

    const size_t M = out_off[NB];

    out.keys   = (uint64_t*)malloc(M * sizeof(uint64_t));
    out.values = (uint32_t*)malloc(M * sizeof(uint32_t));

    #pragma omp parallel for num_threads(nthreads_for_omp) schedule(static)
    for (size_t b = 0; b < NB; ++b) {
        if (!bk_len[b]) continue;
        memcpy(out.keys   + out_off[b], bk_keys[b], bk_len[b] * sizeof(uint64_t));
        memcpy(out.values + out_off[b], bk_vals[b], bk_len[b] * sizeof(uint32_t));
        free(bk_keys[b]);
        free(bk_vals[b]);
    }

    free(bk_keys);
    free(bk_vals);
    free(bk_len);
    free(bucket_off);
    free(out_off);

    out.len = M;
    return out;

CLEAN_BUCKET:
    free(B_keys);
    free(B_vals);

CLEAN_EMPTY:
    for (int t = 0; t < nth; ++t) {
        v_free(&V_thr[t]);
        free(C_thr[t]);
    }
    free(C_thr);
    return (SortedKV_Arrays_t){0};
}

static SortedKV_Arrays_t finalize_pos_vectors_bucketed(u64posvec *V_thr, int nth, int nthreads_for_omp)
{
    SortedKV_Arrays_t out = (SortedKV_Arrays_t){0};

    uint32_t **C_thr = (uint32_t**)calloc((size_t)nth, sizeof(uint32_t*));
    if (!C_thr) return out;

    size_t pre_elems = 0;
    for (int t = 0; t < nth; ++t) pre_elems += V_thr[t].n;

    #pragma omp parallel for if (nth>1 && pre_elems >= (1u<<20)) num_threads(nthreads_for_omp) schedule(dynamic)
    for (int t = 0; t < nth; ++t)
        shrink_thread_posvec(&V_thr[t], &C_thr[t]);

    size_t total_distinct = 0;
    for (int t = 0; t < nth; ++t) total_distinct += V_thr[t].n;

    if (total_distinct == 0) {
        for (int t = 0; t < nth; ++t) {
            pv_free(&V_thr[t]);
            free(C_thr[t]);
        }
        free(C_thr);
        return out;
    }

    const unsigned BITS = clamp_bucket_bits_from_total(total_distinct);
    const size_t NB = (size_t)1u << BITS;

    size_t *bucket_totals = (size_t*)calloc(NB, sizeof(size_t));
    size_t **thr_counts = (size_t**)malloc((size_t)nth * sizeof(size_t*));
    if (!bucket_totals || !thr_counts) goto CLEAN_EMPTY;

    for (int t = 0; t < nth; ++t) {
        thr_counts[t] = (size_t*)calloc(NB, sizeof(size_t));
        if (!thr_counts[t]) goto CLEAN_EMPTY;
    }

    #pragma omp parallel for num_threads(nthreads_for_omp) schedule(static)
    for (int t = 0; t < nth; ++t) {
        uint64_t *A = V_thr[t].keys;
        size_t n = V_thr[t].n;
        for (size_t i = 0; i < n; ++i) {
            unsigned b = (unsigned)(A[i] >> (64 - BITS));
            ++thr_counts[t][b];
        }
    }

    for (size_t b = 0; b < NB; ++b) {
        size_t sum = 0;
        for (int t = 0; t < nth; ++t) sum += thr_counts[t][b];
        bucket_totals[b] = sum;
    }

    size_t *bucket_off = (size_t*)malloc((NB + 1) * sizeof(size_t));
    if (!bucket_off) goto CLEAN_EMPTY;
    bucket_off[0] = 0;
    for (size_t b = 0; b < NB; ++b)
        bucket_off[b + 1] = bucket_off[b] + bucket_totals[b];

    const size_t TOTAL = bucket_off[NB];
    uint64_t *B_keys = (uint64_t*)malloc(TOTAL * sizeof(uint64_t));
    uint32_t *B_vals = (uint32_t*)malloc(TOTAL * sizeof(uint32_t));
    uint64_t *B_pos = (uint64_t*)malloc(TOTAL * sizeof(uint64_t));
    if (!B_keys || !B_vals || !B_pos) goto CLEAN_EMPTY;

    size_t **thr_write = (size_t**)malloc((size_t)nth * sizeof(size_t*));
    if (!thr_write) goto CLEAN_EMPTY;
    for (int t = 0; t < nth; ++t) {
        thr_write[t] = (size_t*)malloc(NB * sizeof(size_t));
        if (!thr_write[t]) goto CLEAN_EMPTY;
    }

    for (size_t b = 0; b < NB; ++b) {
        size_t off = bucket_off[b];
        for (int t = 0; t < nth; ++t) {
            thr_write[t][b] = off;
            off += thr_counts[t][b];
        }
    }

    #pragma omp parallel for num_threads(nthreads_for_omp) schedule(static)
    for (int t = 0; t < nth; ++t) {
        uint64_t *Ak = V_thr[t].keys;
        uint64_t *Ap = V_thr[t].positions;
        uint32_t *Av = C_thr[t];
        size_t n = V_thr[t].n;
        size_t *pos = thr_write[t];
        for (size_t i = 0; i < n; ++i) {
            uint64_t k = Ak[i];
            unsigned b = (unsigned)(k >> (64 - BITS));
            size_t dst = pos[b]++;
            B_keys[dst] = k;
            B_vals[dst] = Av ? Av[i] : 1u;
            B_pos[dst] = Ap[i];
        }
    }

    for (int t = 0; t < nth; ++t) {
        pv_free(&V_thr[t]);
        free(C_thr[t]);
        free(thr_counts[t]);
        free(thr_write[t]);
    }
    free(C_thr);
    free(thr_counts);
    free(thr_write);

    uint64_t **bk_keys = (uint64_t**)malloc(NB * sizeof(uint64_t*));
    uint32_t **bk_vals = (uint32_t**)malloc(NB * sizeof(uint32_t*));
    uint64_t **bk_pos = (uint64_t**)malloc(NB * sizeof(uint64_t*));
    size_t *bk_len = (size_t*)malloc(NB * sizeof(size_t));
    if (!bk_keys || !bk_vals || !bk_pos || !bk_len) goto CLEAN_BUCKET;

    #pragma omp parallel for num_threads(nthreads_for_omp) schedule(dynamic)
    for (size_t b = 0; b < NB; ++b) {
        size_t begin = bucket_off[b];
        size_t end = bucket_off[b + 1];
        size_t n = end - begin;
        if (!n) {
            bk_keys[b] = NULL;
            bk_vals[b] = NULL;
            bk_pos[b] = NULL;
            bk_len[b] = 0;
            continue;
        }
        uint64_t *kseg = B_keys + begin;
        uint32_t *vseg = B_vals + begin;
        uint64_t *pseg = B_pos + begin;
        radix_sort_kvp_u64(kseg, vseg, pseg, n);
        size_t m = shrink_kvp_inplace_u64_u32_pos(kseg, vseg, pseg, n);
        bk_keys[b] = (uint64_t*)malloc(m * sizeof(uint64_t));
        bk_vals[b] = (uint32_t*)malloc(m * sizeof(uint32_t));
        bk_pos[b] = (uint64_t*)malloc(m * sizeof(uint64_t));
        if (!bk_keys[b] || !bk_vals[b] || !bk_pos[b])
            err(errno, "%s(): OOM position bucket", __func__);
        memcpy(bk_keys[b], kseg, m * sizeof(uint64_t));
        memcpy(bk_vals[b], vseg, m * sizeof(uint32_t));
        memcpy(bk_pos[b], pseg, m * sizeof(uint64_t));
        bk_len[b] = m;
    }

    free(B_keys);
    free(B_vals);
    free(B_pos);
    free(bucket_totals);

    size_t *out_off = (size_t*)malloc((NB + 1) * sizeof(size_t));
    if (!out_off) err(errno, "%s(): OOM output offsets", __func__);
    out_off[0] = 0;
    for (size_t b = 0; b < NB; ++b)
        out_off[b + 1] = out_off[b] + bk_len[b];

    const size_t M = out_off[NB];
    out.keys = (uint64_t*)malloc(M * sizeof(uint64_t));
    out.values = (uint32_t*)malloc(M * sizeof(uint32_t));
    out.positions = (uint64_t*)malloc(M * sizeof(uint64_t));
    if (!out.keys || !out.values || !out.positions)
        err(errno, "%s(): OOM positional output arrays", __func__);

    #pragma omp parallel for num_threads(nthreads_for_omp) schedule(static)
    for (size_t b = 0; b < NB; ++b) {
        if (!bk_len[b]) continue;
        memcpy(out.keys + out_off[b], bk_keys[b], bk_len[b] * sizeof(uint64_t));
        memcpy(out.values + out_off[b], bk_vals[b], bk_len[b] * sizeof(uint32_t));
        memcpy(out.positions + out_off[b], bk_pos[b], bk_len[b] * sizeof(uint64_t));
        free(bk_keys[b]);
        free(bk_vals[b]);
        free(bk_pos[b]);
    }

    free(bk_keys);
    free(bk_vals);
    free(bk_pos);
    free(bk_len);
    free(bucket_off);
    free(out_off);

    out.len = M;
    return out;

CLEAN_BUCKET:
    free(B_keys);
    free(B_vals);
    free(B_pos);

CLEAN_EMPTY:
    for (int t = 0; t < nth; ++t) {
        pv_free(&V_thr[t]);
        free(C_thr[t]);
    }
    free(C_thr);
    return (SortedKV_Arrays_t){0};
}


static inline void write_index_payload_and_free(FILE *comb, FILE *comb_ab, FILE *comb_pos,
                bool write_ab, bool write_pos, uint64_t *sketch_index_slot, SortedKV_Arrays_t *kv)
{
    *sketch_index_slot = kv->len;
    if (kv->len) {
        fwrite(kv->keys,   sizeof(uint64_t), kv->len, comb);
        if (write_ab) fwrite(kv->values, sizeof(uint32_t), kv->len, comb_ab);
        if (write_pos) fwrite(kv->positions, sizeof(uint64_t), kv->len, comb_pos);
    }
    free(kv->keys);   kv->keys   = NULL;
    free(kv->values); kv->values = NULL;
    free(kv->positions); kv->positions = NULL;
    kv->len = 0;
}

// new helps for asone mode
static inline void collect_one_file_into_vectors(
    const char *path, int nth, sketch_opt_t *opt, u64vec *V_thr, int BATCH_READS,
    uint64_t **asm_lengths, size_t *asm_n, size_t *asm_cap)
{
    if (uses_stream_input(opt, path) || isCompressfile((char*)path) || isOK_fmt_infile((char*)path, fasta_fmt, FAS_FMT_SZ)) {
        const int nth = (opt->p > 0 ? opt->p : 1);
        collect_gz_into_vectors(path, opt->pipecmd, nth, nth, V_thr,
                                ctxmask, tupmask, Bitslen.obj, klen, BATCH_READS,
                                asm_lengths, asm_n, asm_cap);
    } else {
        collect_plain_mmap_into_vectors(path, nth, opt->p, V_thr,
                                        ctxmask, tupmask, Bitslen.obj, klen,
                                        asm_lengths, asm_n, asm_cap);
    }
}

static inline uint64_t collect_one_file_into_posvectors(
    const char *path, int nth, sketch_opt_t *opt, u64posvec *V_thr, int BATCH_READS,
    uint64_t start_offset,
    uint64_t **asm_lengths, size_t *asm_n, size_t *asm_cap)
{
    return collect_gz_into_posvectors(path, opt->pipecmd, nth, nth, V_thr, start_offset,
                                      ctxmask, tupmask, Bitslen.obj, klen, BATCH_READS,
                                      asm_lengths, asm_n, asm_cap);
}

static inline uint64_t finalize_filter_conflict_and_write(
    FILE *comb, FILE *comb_ab, sketch_opt_t *opt, u64vec *V_thr, int nth,
    kmer_count_filter_result_t *count_cutoff_out)
{
    SortedKV_Arrays_t kv = finalize_vectors_bucketed(V_thr, nth, nth);

    const kmer_count_filter_result_t count_cutoff = apply_kmer_count_filters(&kv, opt);
    if (count_cutoff_out) *count_cutoff_out = count_cutoff;

    if (kv.len && !opt->conflict) {
        if (opt->abundance) remove_ctx_with_conflict_obj(&kv, Bitslen.obj);
        else remove_ctx_with_conflict_obj_noabund(kv.keys, &(kv.len), Bitslen.obj);
    }

    uint64_t bytes = 0;
    if (kv.len) {
        write_index_payload_and_free(comb, (opt->abundance ? comb_ab : NULL), NULL,
                                     opt->abundance, false, &bytes, &kv);
    }
    return bytes;
}

static inline uint64_t finalize_pos_filter_conflict_and_write(
    FILE *comb, FILE *comb_ab, FILE *comb_pos, sketch_opt_t *opt, u64posvec *V_thr, int nth,
    kmer_count_filter_result_t *count_cutoff_out)
{
    SortedKV_Arrays_t kv = finalize_pos_vectors_bucketed(V_thr, nth, nth);

    const kmer_count_filter_result_t count_cutoff = apply_kmer_count_filters(&kv, opt);
    if (count_cutoff_out) *count_cutoff_out = count_cutoff;

    if (kv.len && !opt->conflict) {
        if (opt->abundance || kv.values) remove_ctx_with_conflict_obj(&kv, Bitslen.obj);
        else remove_ctx_with_conflict_obj_noabund_pos(&kv, Bitslen.obj);
    }

    uint64_t bytes = 0;
    if (kv.len) {
        write_index_payload_and_free(comb, (opt->abundance ? comb_ab : NULL), comb_pos,
                                     opt->abundance, true, &bytes, &kv);
    } else {
        free(kv.keys);
        free(kv.values);
        free(kv.positions);
    }
    return bytes;
}

static void sketch_few_files_with_intrafile_parallel_pos(sketch_opt_t *opt, infile_tab_t *tab, int BATCH_READS)
{
    const int nfiles = tab->infile_num;
    const int nth = (opt->p > 0 ? opt->p : 1);
    const size_t sample_qc_count = opt->asone ? 1u : (size_t)nfiles;

    dim_sketch_qc_stat_t *sample_qc_stats = NULL;
    if (sketch_records_sample_qc(opt)) {
        sample_qc_stats = (dim_sketch_qc_stat_t *)calloc(sample_qc_count, sizeof(sample_qc_stats[0]));
        if (!sample_qc_stats) err(errno, "%s(): OOM sample_qc_stats", __func__);
    }
    infile_meta_t *sample_infile_meta = NULL;
    if (opt->compute_meta) {
        sample_infile_meta = (infile_meta_t *)calloc(sample_qc_count, sizeof(sample_infile_meta[0]));
        if (!sample_infile_meta) err(errno, "%s(): OOM sample_infile_meta", __func__);
    }

    FILE *comb = fopen(format_string("%s/%s", opt->outdir, combined_sketch_suffix), "wb");
    if (!comb) err(errno, "%s() open file error: %s/%s", __func__, opt->outdir, combined_sketch_suffix);
    setvbuf(comb, NULL, _IOFBF, 8u << 20);

    FILE *comb_ab = NULL;
    if (opt->abundance) {
        comb_ab = fopen(format_string("%s/%s", opt->outdir, combined_ab_suffix), "wb");
        if (!comb_ab) err(errno, "%s() open file error: %s/%s", __func__, opt->outdir, combined_ab_suffix);
        setvbuf(comb_ab, NULL, _IOFBF, 8u << 20);
    }

    FILE *comb_pos = fopen(format_string("%s/%s", opt->outdir, sketch_position_suffix), "wb");
    if (!comb_pos) err(errno, "%s() open file error: %s/%s", __func__, opt->outdir, sketch_position_suffix);
    setvbuf(comb_pos, NULL, _IOFBF, 8u << 20);

    u64posvec *V_thr = (u64posvec*)malloc((size_t)nth * sizeof(u64posvec));
    if (!V_thr) err(errno, "%s(): OOM V_thr", __func__);
    for (int t = 0; t < nth; ++t) pv_init(&V_thr[t], 1u << 15);

    if (opt->asone) {
        uint8_t flags = KSSD_INFILE_FLAG_MIXED_FORMAT;
        int8_t fmt = KSSD_INFILE_FMT_UNKNOWN;
        for (int f = 0; f < nfiles; ++f) {
            const char *path = tab->organized_infile_tab[f].fpath;
            const int8_t this_fmt = infile_fmt_from_path(path);
            flags |= infile_flags_from_path(path, opt);
            if (f == 0)
                fmt = this_fmt;
            else if (fmt != this_fmt)
                fmt = KSSD_INFILE_FMT_UNKNOWN;
        }
        if (fmt != KSSD_INFILE_FMT_UNKNOWN)
            flags &= (uint8_t)~KSSD_INFILE_FLAG_MIXED_FORMAT;
        const bool collect_asone_meta_lengths =
            opt->compute_meta && fmt == KSSD_INFILE_FMT_FASTA;
        uint64_t *asm_lengths = NULL;
        size_t asm_n = 0, asm_cap = 0;
        uint64_t stream_offset = 0;

        for (int f = 0; f < nfiles; ++f) {
            const char *path = tab->organized_infile_tab[f].fpath;
            stream_offset = collect_one_file_into_posvectors(
                path, nth, opt, V_thr, BATCH_READS, stream_offset,
                collect_asone_meta_lengths ? &asm_lengths : NULL,
                collect_asone_meta_lengths ? &asm_n : NULL,
                collect_asone_meta_lengths ? &asm_cap : NULL);
            fprintf(stderr, "\r%d/%d parts processed (as-one)", f + 1, nfiles);
            fflush(stderr);
        }
        fprintf(stderr, "\n");

        kmer_count_filter_result_t count_cutoff = {0};
        uint64_t bytes = finalize_pos_filter_conflict_and_write(comb, comb_ab, comb_pos, opt, V_thr, nth, &count_cutoff);
        if (sample_qc_stats)
            sample_qc_stats[0] = count_cutoff.sample_qc;
        if (sketch_reports_count_filter(opt))
            print_kmer_count_filter_result(opt, "as-one sketch", &count_cutoff, false);
        if (opt->compute_meta)
            infile_meta_from_lengths(asm_lengths, asm_n, fmt, create_type_from_opt(opt),
                                     flags, &sample_infile_meta[0]);
        free(asm_lengths);

        uint64_t idx2[2] = {0, bytes};
        write_to_file(format_string("%s/%s", opt->outdir, idx_sketch_suffix), idx2, sizeof(idx2));

        fclose(comb);
        if (opt->abundance) fclose(comb_ab);
        fclose(comb_pos);
        free(V_thr);

        infile_tab_t tab_asone = *tab;
        tab_asone.infile_num = 1;
        write_sketch_stat_ex(opt->outdir, &tab_asone, opt->anno && !sketch_has_stream_input(opt, tab), false);
        if (opt->compute_meta)
            write_sketch_infile_meta_stats(opt->outdir, sample_infile_meta, sample_qc_count);
        write_sketch_qc_stats(opt->outdir, sample_qc_stats, sample_qc_count);
        free(sample_infile_meta);
        free(sample_qc_stats);
        return;
    }

    uint64_t *sketch_index = (uint64_t*)calloc((size_t)nfiles + 1, sizeof(uint64_t));
    if (!sketch_index) err(errno, "%s(): OOM sketch_index", __func__);

    for (int f = 0; f < nfiles; ++f) {
        for (int t = 0; t < nth; ++t) {
            if (V_thr[t].keys == NULL || V_thr[t].cap == 0) pv_init(&V_thr[t], 1u << 15);
            else V_thr[t].n = 0;
        }

        const char *path = tab->organized_infile_tab[f].fpath;
        uint64_t *asm_lengths = NULL;
        size_t asm_n = 0, asm_cap = 0;
        const bool collect_meta_lengths = should_collect_lengths_for_meta(opt, path);
        collect_one_file_into_posvectors(path, nth, opt, V_thr, BATCH_READS, 0,
                                         collect_meta_lengths ? &asm_lengths : NULL,
                                         collect_meta_lengths ? &asm_n : NULL,
                                         collect_meta_lengths ? &asm_cap : NULL);
        if (opt->compute_meta)
            infile_meta_from_lengths(asm_lengths, asm_n, infile_fmt_from_path(path),
                                     create_type_from_opt(opt),
                                     infile_flags_from_path(path, opt),
                                     &sample_infile_meta[f]);
        free(asm_lengths);

        kmer_count_filter_result_t count_cutoff = {0};
        sketch_index[f + 1] = finalize_pos_filter_conflict_and_write(comb, comb_ab, comb_pos, opt, V_thr, nth, &count_cutoff);
        if (sample_qc_stats)
            sample_qc_stats[f] = count_cutoff.sample_qc;
        if (sketch_reports_count_filter(opt))
            print_kmer_count_filter_result(opt, path, &count_cutoff, true);

        fprintf(stderr, "\r%d/%d genomes sketched", f + 1, nfiles);
        fflush(stderr);
    }
    fprintf(stderr, "\n");

    for (int i = 0; i < nfiles; ++i) sketch_index[i + 1] += sketch_index[i];
    write_to_file(format_string("%s/%s", opt->outdir, idx_sketch_suffix),
                  sketch_index, (size_t)(nfiles + 1) * sizeof(uint64_t));

    fclose(comb);
    if (opt->abundance) fclose(comb_ab);
    fclose(comb_pos);
    free(sketch_index);

    for (int t = 0; t < nth; ++t) pv_free(&V_thr[t]);
    free(V_thr);

    write_sketch_stat_ex(opt->outdir, tab, opt->anno && !sketch_has_stream_input(opt, tab), false);
    if (opt->compute_meta)
        write_sketch_infile_meta_stats(opt->outdir, sample_infile_meta, sample_qc_count);
    write_sketch_qc_stats(opt->outdir, sample_qc_stats, sample_qc_count);
    free(sample_infile_meta);
    free(sample_qc_stats);
}
// =================== Mode B (in-file parallel; gz + mmap) ===================
static void sketch_few_files_with_intrafile_parallel(sketch_opt_t *opt, infile_tab_t *tab, int BATCH_READS)
{
    if (opt->position) {
        sketch_few_files_with_intrafile_parallel_pos(opt, tab, BATCH_READS);
        return;
    }
    const int nfiles   = tab->infile_num;
    const int nth      = (opt->p > 0 ? opt->p : 1);
    const size_t sample_qc_count = opt->asone ? 1u : (size_t)nfiles;
    dim_sketch_qc_stat_t *sample_qc_stats = NULL;
    if (sketch_records_sample_qc(opt)) {
        sample_qc_stats =
            (dim_sketch_qc_stat_t *)calloc(sample_qc_count, sizeof(sample_qc_stats[0]));
        if (!sample_qc_stats) err(errno, "%s(): OOM sample_qc_stats", __func__);
    }
    infile_meta_t *sample_infile_meta = NULL;
    if (opt->compute_meta) {
        sample_infile_meta =
            (infile_meta_t *)calloc(sample_qc_count, sizeof(sample_infile_meta[0]));
        if (!sample_infile_meta) err(errno, "%s(): OOM sample_infile_meta", __func__);
    }

    FILE *comb = fopen(format_string("%s/%s", opt->outdir, combined_sketch_suffix), "wb");
    if (!comb) err(errno, "%s() open file error: %s/%s", __func__, opt->outdir, combined_sketch_suffix);
    setvbuf(comb, NULL, _IOFBF, 8u << 20);

    FILE *comb_ab = NULL;
    if (opt->abundance) {
        comb_ab = fopen(format_string("%s/%s", opt->outdir, combined_ab_suffix), "wb");
        if (!comb_ab) err(errno, "%s() open file error: %s/%s", __func__, opt->outdir, combined_ab_suffix);
        setvbuf(comb_ab, NULL, _IOFBF, 8u << 20);
    }

    // allocate thread vectors once
    u64vec *V_thr = (u64vec*)malloc((size_t)nth * sizeof(u64vec));
    if (!V_thr) err(errno, "%s(): OOM V_thr", __func__);
    for (int t = 0; t < nth; ++t) v_init(&V_thr[t], 1u << 15);

    // -------------------- AS-ONE MODE --------------------
    if (opt->asone) {
        uint8_t flags = KSSD_INFILE_FLAG_MIXED_FORMAT;
        int8_t fmt = KSSD_INFILE_FMT_UNKNOWN;
        for (int f = 0; f < nfiles; ++f) {
            const char *path = tab->organized_infile_tab[f].fpath;
            const int8_t this_fmt = infile_fmt_from_path(path);
            flags |= infile_flags_from_path(path, opt);
            if (f == 0)
                fmt = this_fmt;
            else if (fmt != this_fmt)
                fmt = KSSD_INFILE_FMT_UNKNOWN;
        }
        if (fmt != KSSD_INFILE_FMT_UNKNOWN)
            flags &= (uint8_t)~KSSD_INFILE_FLAG_MIXED_FORMAT;
        const bool collect_asone_meta_lengths =
            opt->compute_meta && fmt == KSSD_INFILE_FMT_FASTA;
        uint64_t *asm_lengths = NULL;
        size_t asm_n = 0, asm_cap = 0;
        // do NOT clear between files: append all parts into one logical input
        for (int f = 0; f < nfiles; ++f) {
            const char *path = tab->organized_infile_tab[f].fpath;
            collect_one_file_into_vectors(path, nth, opt, V_thr, BATCH_READS,
                                          collect_asone_meta_lengths ? &asm_lengths : NULL,
                                          collect_asone_meta_lengths ? &asm_n : NULL,
                                          collect_asone_meta_lengths ? &asm_cap : NULL);
            fprintf(stderr, "\r%d/%d parts processed (as-one)", f+1, nfiles);
            fflush(stderr);
        }
        fprintf(stderr, "\n");

        kmer_count_filter_result_t count_cutoff = {0};
        uint64_t bytes = finalize_filter_conflict_and_write(comb, comb_ab, opt, V_thr, nth, &count_cutoff);
        if (sample_qc_stats)
            sample_qc_stats[0] = count_cutoff.sample_qc;
        if (sketch_reports_count_filter(opt)) {
            print_kmer_count_filter_result(opt, "as-one sketch", &count_cutoff, false);
        }
        if (opt->compute_meta) {
            infile_meta_from_lengths(asm_lengths, asm_n, fmt, create_type_from_opt(opt),
                                     flags, &sample_infile_meta[0]);
        }
        free(asm_lengths);

        // index for asone has only 2 entries: [0, total_bytes]
        uint64_t idx2[2] = {0, bytes};
        write_to_file(format_string("%s/%s", opt->outdir, idx_sketch_suffix), idx2, sizeof(idx2));

        fclose(comb);
        if (opt->abundance) fclose(comb_ab);

        // finalize_vectors_bucketed() already v_free()'d each V_thr[t], so only free the struct array
        free(V_thr);
        V_thr = NULL;

        // stats should see only 1 logical input, first filename is representative
        infile_tab_t tab_asone = *tab;
        tab_asone.infile_num = 1;
        write_sketch_stat_ex(opt->outdir, &tab_asone, opt->anno && !sketch_has_stream_input(opt, tab), false);
        if (opt->compute_meta)
            write_sketch_infile_meta_stats(opt->outdir, sample_infile_meta, sample_qc_count);
        write_sketch_qc_stats(opt->outdir, sample_qc_stats, sample_qc_count);
        free(sample_infile_meta);
        free(sample_qc_stats);
        return;
    }

    // -------------------- NORMAL MODE (per-file) --------------------
    uint64_t *sketch_index = (uint64_t*)calloc((size_t)nfiles + 1, sizeof(uint64_t));
    if (!sketch_index) err(errno, "%s(): OOM sketch_index", __func__);

    for (int f = 0; f < nfiles; ++f) {
        // After each finalize, V_thr[t].a might be freed (a=NULL, cap=0). Re-init if needed.
        for (int t = 0; t < nth; ++t) {
            if (V_thr[t].a == NULL || V_thr[t].cap == 0) v_init(&V_thr[t], 1u << 15);
            else V_thr[t].n = 0;
        }

        const char *path = tab->organized_infile_tab[f].fpath;
        uint64_t *asm_lengths = NULL;
        size_t asm_n = 0, asm_cap = 0;
        const bool collect_meta_lengths = should_collect_lengths_for_meta(opt, path);
        collect_one_file_into_vectors(path, nth, opt, V_thr, BATCH_READS,
                                      collect_meta_lengths ? &asm_lengths : NULL,
                                      collect_meta_lengths ? &asm_n : NULL,
                                      collect_meta_lengths ? &asm_cap : NULL);
        if (opt->compute_meta)
            infile_meta_from_lengths(asm_lengths, asm_n, infile_fmt_from_path(path),
                                     create_type_from_opt(opt),
                                     infile_flags_from_path(path, opt),
                                     &sample_infile_meta[f]);
        free(asm_lengths);

        kmer_count_filter_result_t count_cutoff = {0};
        sketch_index[f+1] = finalize_filter_conflict_and_write(comb, comb_ab, opt, V_thr, nth, &count_cutoff);
        if (sample_qc_stats)
            sample_qc_stats[f] = count_cutoff.sample_qc;
        if (sketch_reports_count_filter(opt)) {
            print_kmer_count_filter_result(opt, path, &count_cutoff, true);
        }

        fprintf(stderr, "\r%d/%d genomes sketched", f+1, nfiles);
        fflush(stderr);
    }
    fprintf(stderr, "\n");

    for (int i = 0; i < nfiles; ++i) sketch_index[i+1] += sketch_index[i];
    write_to_file(format_string("%s/%s", opt->outdir, idx_sketch_suffix),
                  sketch_index, (size_t)(nfiles + 1) * sizeof(uint64_t));

    fclose(comb);
    if (opt->abundance) fclose(comb_ab);

    free(sketch_index);

    for (int t = 0; t < nth; ++t) v_free(&V_thr[t]);  // safe even if already freed if v_free handles NULL
    free(V_thr);
    V_thr = NULL;

    write_sketch_stat_ex(opt->outdir, tab, opt->anno && !sketch_has_stream_input(opt, tab), false);
    if (opt->compute_meta)
        write_sketch_infile_meta_stats(opt->outdir, sample_infile_meta, sample_qc_count);
    write_sketch_qc_stats(opt->outdir, sample_qc_stats, sample_qc_count);
    free(sample_infile_meta);
    free(sample_qc_stats);
}


/* old version not support asone mode
static void sketch_few_files_with_intrafile_parallel(sketch_opt_t *opt, infile_tab_t *tab, int BATCH_READS)
{
    const int nfiles = tab->infile_num;
    const int kmerocrs = opt->kmerocrs;
    FILE *comb = fopen(format_string("%s/%s", opt->outdir, combined_sketch_suffix), "wb");
    if (!comb) err(errno, "%s() open file error: %s/%s", __func__, opt->outdir, combined_sketch_suffix);
    setvbuf(comb, NULL, _IOFBF, 8u << 20);

    FILE *comb_ab = NULL;
    if (opt->abundance) {
        comb_ab = fopen(format_string("%s/%s", opt->outdir, combined_ab_suffix), "wb");
        if (!comb_ab) err(errno, "%s() open file error: %s/%s", __func__, opt->outdir, combined_ab_suffix);
        setvbuf(comb_ab, NULL, _IOFBF, 8u<<20);
    }

    uint64_t *sketch_index = (uint64_t*)calloc((size_t)nfiles + 1, sizeof(uint64_t));
    if (!sketch_index) err(errno, "%s(): OOM sketch_index", __func__);

    for (int f=0; f<nfiles; ++f){
        const char *path = tab->organized_infile_tab[f].fpath;
        const int nth = (opt->p > 0 ? opt->p : 1);

        // per-thread accumulators
        u64vec *V_thr = (u64vec*)malloc((size_t)nth*sizeof(u64vec));
        if (!V_thr) err(errno, "%s(): OOM V_thr", __func__);
        for (int t=0;t<nth;++t) v_init(&V_thr[t], 1u<<15);

        // 1) collect (mmap for plain; kseq for .gz)
        if (uses_stream_input(opt, path) || isCompressfile((char*)path) || isOK_fmt_infile((char*)path, fasta_fmt, FAS_FMT_SZ))
            collect_gz_into_vectors(path, opt->pipecmd, nth, opt->p, V_thr, ctxmask, tupmask, Bitslen.obj, klen, BATCH_READS,
                                    NULL, NULL, NULL);
        else // collect_plain_mmap_into_vectors() only for plain FASTQ not fast
            collect_plain_mmap_into_vectors(path, nth, opt->p, V_thr, ctxmask, tupmask, Bitslen.obj, klen,
                                            NULL, NULL, NULL);
        

        // 2) finalize once (shared shrink + bucketed merge)
        SortedKV_Arrays_t kv = finalize_vectors_bucketed(V_thr, nth, opt->p);
        if(kmerocrs > 1) // filter by min ocrs
                    filter_n_SortedKV_Arrays(&kv, kmerocrs);
        // 3) conflict filter + 4) write
        if (kv.len){
            if (!opt->conflict){
                    if(opt->abundance) remove_ctx_with_conflict_obj(&kv, Bitslen.obj);
                    else remove_ctx_with_conflict_obj_noabund (kv.keys, &(kv.len), Bitslen.obj);
            }
            write_index_payload_and_free(comb, (opt->abundance?comb_ab:NULL),
                                         opt->abundance, &sketch_index[f+1], &kv);
        } else {
            sketch_index[f+1] = 0;
        }

        fprintf(stderr, "\r%d/%d genomes sketched", f+1, nfiles);
        fflush(stderr);
    }

    for (int i=0;i<nfiles;++i) sketch_index[i+1] += sketch_index[i];
    write_to_file(format_string("%s/%s", opt->outdir, idx_sketch_suffix),
                  sketch_index, (size_t)(nfiles + 1) * sizeof(uint64_t));
    fclose(comb);
    if (opt->abundance) fclose(comb_ab);
    free(sketch_index);
    write_sketch_stat(opt->outdir, tab, opt->anno);
}
*/
// simple API for kssd3 ani use directly
// Faster, memory-only, keys-only sketcher (vector + radix + dedup + conflict filter)
// API kept identical to your original.

simple_sketch_t *simple_genomes2mem2sortedctxobj64_mem(infile_tab_t *infile_stat, int drfold)
{
    // ---- 1) configure sketch filter (same semantics as your original) ----
    FILTER = UINT32_MAX >> drfold;

    const int nfiles = infile_stat->infile_num;
    if (nfiles <= 0) {
        simple_sketch_t *ret = (simple_sketch_t*)calloc(1, sizeof(simple_sketch_t));
        if (!ret) err(EXIT_FAILURE, "%s(): OOM ret", __func__);
        ret->comb_sketch = NULL;
        ret->sketch_index = (uint64_t*)calloc(1, sizeof(uint64_t)); // [0] = 0
        ret->infile_num = 0;
        return ret;
    }

    // ---- 2) per-file outputs (keys only) ----
    uint64_t **per_keys = (uint64_t**)calloc((size_t)nfiles, sizeof(uint64_t*));
    uint64_t  *per_len  = (uint64_t*) calloc((size_t)nfiles, sizeof(uint64_t));
    if (!per_keys || !per_len) err(EXIT_FAILURE, "%s(): OOM per-file arrays", __func__);

    // ---- 3) parallel across files: load → vectorize (ctxobj) → sort → dedup → conflict-filter ----
    // Set number of threads via OMP env or your build flags (-fopenmp)
    #pragma omp parallel for schedule(dynamic,1) num_threads(nfiles)
    for (int f = 0; f < nfiles; ++f) {
        char *path = infile_stat->organized_infile_tab[f].fpath;

        // Collect kept ctxobj keys into a single vector for this file
        u64vec vec; v_init(&vec, 1u << 15); // heuristic starting capacity
        load_genome_into_single_vec(path, NULL, &vec, ctxmask, tupmask, Bitslen.obj, klen,
                                    NULL, NULL, NULL);
        if (vec.n == 0) {
            per_keys[f] = NULL;
            per_len[f]  = 0;
            v_free(&vec);
            continue;
        }

        // Sort + dedup (keys only)
        radix_sort_u64(vec.a, vec.n);
        vec.n = dedup_sorted_uint64(vec.a, vec.n);
        // Conflict filter (no abundance version) — keeps array sorted, shrinks length in-place
        remove_ctx_with_conflict_obj_noabund(vec.a, &vec.n, Bitslen.obj);

        // Transfer ownership: shrink to exact size and hand out
        uint64_t *outk = (uint64_t*)malloc(vec.n * sizeof(uint64_t));
        if (!outk) err(EXIT_FAILURE, "%s(): OOM outk", __func__);
        memcpy(outk, vec.a, vec.n * sizeof(uint64_t));
        per_keys[f] = outk;
        per_len[f]  = (uint64_t)vec.n;

        v_free(&vec);
    }

    // ---- 4) prefix-sum index & concatenate all into one combined sketch buffer ----
    uint64_t *sketch_index = (uint64_t*)calloc((size_t)nfiles + 1, sizeof(uint64_t));
    if (!sketch_index) err(EXIT_FAILURE, "%s(): OOM sketch_index", __func__);

    for (int i = 0; i < nfiles; ++i) sketch_index[i + 1] = sketch_index[i] + per_len[i];
    const uint64_t total_keys = sketch_index[nfiles];

    uint64_t *combined = NULL;
    if (total_keys) {
        combined = (uint64_t*)malloc((size_t)total_keys * sizeof(uint64_t));
        if (!combined) err(EXIT_FAILURE, "%s(): OOM combined", __func__);

        // copy each file's keys to its slot
        for (int i = 0; i < nfiles; ++i) {
            const uint64_t n = per_len[i];
            if (!n) continue;
            memcpy(combined + sketch_index[i], per_keys[i], (size_t)n * sizeof(uint64_t));
            free(per_keys[i]); // done with this piece
        }
    }
    free(per_keys);
    free(per_len);

    // ---- 5) build return object ----
    simple_sketch_t *ret = (simple_sketch_t*)malloc(sizeof(simple_sketch_t));
    if (!ret) err(EXIT_FAILURE, "%s(): OOM ret", __func__);
    ret->comb_sketch = combined;      // length = total_keys (may be 0)
    ret->sketch_index = sketch_index; // length = nfiles+1
    ret->infile_num = nfiles;
    return ret;
}



typedef struct {
    char *name;
    char *seq;
    int   len;
    char  header[PATHLEN];
} mfa_seq_t;

#define MFA_SEQ_BATCH 128

void mfa2sortedctxobj64_v2 (sketch_opt_t *sketch_opt_val, infile_tab_t *infile_stat)
{
    const bool resolve_conf = !sketch_opt_val->conflict;
    const int  nthreads     = sketch_opt_val->p;

    uint64_t totle_sketch_size = 0;   // cumulative across all sequences

    // ---- global index: one entry per "genome" (sequence), plus 0 at the front ----
    u64vec sketch_index;
    v_init(&sketch_index, 1024);      // initial capacity
    v_push(&sketch_index, 0ULL);      // index[0] = 0

    // tmpname: one name per genome (sequence)
    int tmpname_size_alloc = (int)sketch_index.cap;
    char (*tmpname)[PATHLEN] = malloc((size_t)tmpname_size_alloc * PATHLEN);
    if (!tmpname)
        err(errno, "%s(): Memory allocation failed for tmpname", __func__);
    char (*tmpanno)[PATHLEN] = sketch_opt_val->anno ? calloc((size_t)tmpname_size_alloc, PATHLEN) : NULL;
    if (sketch_opt_val->anno && !tmpanno)
        err(errno, "%s(): Memory allocation failed for tmpanno", __func__);

    dim_sketch_qc_stat_t *sample_qc_stats = NULL;
    size_t sample_qc_len = 0, sample_qc_cap = 0;
    infile_meta_t *sample_infile_meta = NULL;
    size_t sample_meta_len = 0, sample_meta_cap = 0;

    // combined sketch file (keys only)
    FILE *comb_sketch_fp =
        fopen(format_string("%s/%s", sketch_opt_val->outdir, combined_sketch_suffix), "wb");
    if (!comb_sketch_fp)
        err(errno, "%s() open file error: %s/%s",
            __func__, sketch_opt_val->outdir, combined_sketch_suffix);
    setvbuf(comb_sketch_fp, NULL, _IOFBF, 8u << 20);

    FILE *comb_pos_fp = NULL;
    if (sketch_opt_val->position) {
        comb_pos_fp =
            fopen(format_string("%s/%s", sketch_opt_val->outdir, sketch_position_suffix), "wb");
        if (!comb_pos_fp)
            err(errno, "%s() open file error: %s/%s",
                __func__, sketch_opt_val->outdir, sketch_position_suffix);
        setvbuf(comb_pos_fp, NULL, _IOFBF, 8u << 20);
    }

    // =========================
    //  loop over MFA *files*
    // =========================
    for (int fi = 0; fi < infile_stat->infile_num; ++fi) {
        char *seqfname = infile_stat->organized_infile_tab[fi].fpath;

        gzFile infile = gzopen(seqfname, "r");
        if (!infile)
            err(errno, "mfa2sortedctxobj64(): Cannot open file %s", seqfname);

        kseq_t *seq = kseq_init(infile);

        for (;;) {
            mfa_seq_t seqs[MFA_SEQ_BATCH];
            int       nseq = 0;
            int       ret;
            bool      eof = false;

            // -------- read up to MFA_SEQ_BATCH usable sequences into this batch --------
            while (nseq < MFA_SEQ_BATCH) {
                ret = kseq_read(seq);
                if (ret < 0) {  // EOF or error
                    eof = true;
                    break;
                }

                if (seq->seq.l <= klen)
                    continue;   // too short for any k-mer; don't treat as genome

                seqs[nseq].len = (int)seq->seq.l;

                // copy sequence
                seqs[nseq].seq = (char *)malloc((size_t)seq->seq.l);
                if (!seqs[nseq].seq)
                    err(errno, "%s(): OOM for seq string", __func__);
                memcpy(seqs[nseq].seq, seq->seq.s, (size_t)seq->seq.l);

                // copy name
                seqs[nseq].name = strdup(seq->name.s);
                if (!seqs[nseq].name)
                    err(errno, "%s(): OOM for seq name", __func__);
                if (sketch_opt_val->anno)
                    sketch_annotation_copy(seqs[nseq].header, seq->name.s,
                                           seq->comment.l ? seq->comment.s : NULL);

                ++nseq;
            }

            if (nseq == 0) {
                // no usable sequences in this batch
                if (eof) break;  // finished this file
                else continue;   // only short sequences encountered; keep reading
            }

            // -------- parallel sketch & KV build per sequence in this batch --------
            uint64_t **keys = (uint64_t **)calloc((size_t)nseq, sizeof(uint64_t *));
            uint64_t **positions = sketch_opt_val->position
                ? (uint64_t **)calloc((size_t)nseq, sizeof(uint64_t *))
                : NULL;
            uint64_t  *lens = (uint64_t  *)calloc((size_t)nseq, sizeof(uint64_t));
            if (!keys || !lens || (sketch_opt_val->position && !positions))
                err(errno, "%s(): OOM for per-sequence KV arrays", __func__);
            kmer_count_filter_result_t *count_cutoffs = NULL;
            if (sketch_reports_count_filter(sketch_opt_val) || sketch_records_sample_qc(sketch_opt_val)) {
                count_cutoffs = (kmer_count_filter_result_t *)calloc((size_t)nseq, sizeof(kmer_count_filter_result_t));
                if (!count_cutoffs)
                    err(errno, "%s(): OOM for per-sequence count cutoffs", __func__);
            }

            #pragma omp parallel for num_threads(nthreads) schedule(dynamic, 1)
            for (int si = 0; si < nseq; ++si) {
                SortedKV_Arrays_t kv = (SortedKV_Arrays_t){0};
                if (sketch_opt_val->position) {
                    u64posvec pvec;
                    pv_init(&pvec, 1u << 15);
                    sketch_read_into_posvec(seqs[si].seq, seqs[si].len, &pvec, 0,
                                            ctxmask, tupmask, Bitslen.obj, klen);
                    if (pvec.n)
                        kv = build_kv_from_posvec(&pvec, sketch_needs_kmer_counts(sketch_opt_val));
                    pv_free(&pvec);
                } else {
                    u64vec vec;
                    v_init(&vec, 1u << 15);
                    sketch_read_into_vec(seqs[si].seq,seqs[si].len,&vec,
                                        ctxmask,tupmask,Bitslen.obj,klen);
                    if (vec.n)
                        kv = build_kv_from_vec(&vec, sketch_needs_kmer_counts(sketch_opt_val));
                    v_free(&vec);
                }

                if (kv.len) {
                    const kmer_count_filter_result_t count_cutoff = apply_kmer_count_filters(&kv, sketch_opt_val);
                    if (count_cutoffs) count_cutoffs[si] = count_cutoff;

                    if (resolve_conf) {
                        if (kv.values)
                            remove_ctx_with_conflict_obj(&kv, Bitslen.obj);
                        else if (kv.positions)
                            remove_ctx_with_conflict_obj_noabund_pos(&kv, Bitslen.obj);
                        else
                            remove_ctx_with_conflict_obj_noabund(kv.keys, &kv.len, Bitslen.obj);
                    }

                    keys[si] = kv.keys;
                    if (positions)
                        positions[si] = kv.positions;
                    lens[si] = kv.len;

                    free(kv.values);  // counts not kept in MFA mode
                } else {
                    lens[si] = 0;
                    free(kv.keys);
                    free(kv.values);
                    free(kv.positions);
                }
            }

            // -------- serial write in original sequence order for this batch --------
            for (int si = 0; si < nseq; ++si) {
                // ensure tmpname capacity (sketch_index.n == #genomes_so_far + 1)
                if ((int)sketch_index.n >= tmpname_size_alloc) {
                    tmpname_size_alloc += 1000;
                    char (*newtmp)[PATHLEN] =
                        realloc(tmpname, (size_t)tmpname_size_alloc * PATHLEN);
                    if (!newtmp)
                        err(errno, "%s(): Realloc failed for tmpname", __func__);
                    tmpname = newtmp;
                    if (tmpanno) {
                        char (*newanno)[PATHLEN] =
                            realloc(tmpanno, (size_t)tmpname_size_alloc * PATHLEN);
                        if (!newanno)
                            err(errno, "%s(): Realloc failed for tmpanno", __func__);
                        tmpanno = newanno;
                    }
                }

                // record this sequence name (one per genome)
               // replace_special_chars_with_underscore(seqs[si].name);
                snprintf(tmpname[sketch_index.n - 1], PATHLEN, "%s", seqs[si].name);
                if (tmpanno)
                    memcpy(tmpanno[sketch_index.n - 1], seqs[si].header, PATHLEN);

                // write keys for this sequence
                if (count_cutoffs && sketch_reports_count_filter(sketch_opt_val)) {
                    print_kmer_count_filter_result(sketch_opt_val, seqs[si].name,
                                                   &count_cutoffs[si], true);
                }
                if (count_cutoffs && sketch_records_sample_qc(sketch_opt_val)) {
                    append_sketch_qc_stat(&sample_qc_stats, &sample_qc_len, &sample_qc_cap,
                                          count_cutoffs[si].sample_qc);
                }
                if (sketch_opt_val->compute_meta) {
                    infile_meta_t infile_meta = {
                        .total_length_bp = (uint64_t)seqs[si].len,
                        .record_count = 1,
                        .median_length_bp = (uint32_t)seqs[si].len,
                        .asm_level = 1.0f,
                        .length_cv = 0.0f,
                        .meta_fmt_version = KSSD_INFILE_META_VERSION,
                        .infile_fmt = KSSD_INFILE_FMT_FASTA,
                        .create_type = KSSD_CREATE_SPLITMFA,
                        .infile_flags = infile_flags_from_path(seqfname, sketch_opt_val),
                    };
                    append_sketch_infile_meta_stat(&sample_infile_meta, &sample_meta_len, &sample_meta_cap,
                                                   infile_meta);
                }
                if (lens[si]) {
                    fwrite(keys[si], sizeof(uint64_t), lens[si], comb_sketch_fp);
                    if (comb_pos_fp)
                        fwrite(positions[si], sizeof(uint64_t), lens[si], comb_pos_fp);
                }

                totle_sketch_size += lens[si];
                v_push(&sketch_index, totle_sketch_size);

                free(keys[si]);
                if (positions)
                    free(positions[si]);
            }

            free(keys);
            free(positions);
            free(lens);
            free(count_cutoffs);

            // free per-sequence buffers for this batch
            for (int si = 0; si < nseq; ++si) {
                free(seqs[si].seq);
                free(seqs[si].name);
            }

            if (eof) break;  // we've reached end of file
        } // end per-file batches

        kseq_destroy(seq);
        gzclose(infile);

        fprintf(stderr, "\r%dth/%d multifasta file %s completed!\t #genomes=%lu",
                fi + 1,
                infile_stat->infile_num,
                infile_stat->organized_infile_tab[fi].fpath,
                (unsigned long)(sketch_index.n - 1));
        if (fi == infile_stat->infile_num - 1) fprintf(stderr, "\n");
    } // for each file

    fclose(comb_sketch_fp);
    if (comb_pos_fp)
        fclose(comb_pos_fp);

    // write global index & stat
    write_to_file(test_create_fullpath(sketch_opt_val->outdir, idx_sketch_suffix),
                  sketch_index.a,sketch_index.n * sizeof(uint64_t));

    comblco_stat_one.infile_num = (uint32_t)(sketch_index.n - 1);
    concat_and_write_to_file(test_create_fullpath(sketch_opt_val->outdir, sketch_stat),
                             &comblco_stat_one,sizeof(comblco_stat_one),tmpname,
                             (size_t)comblco_stat_one.infile_num * PATHLEN);
    if (tmpanno)
        write_sketch_annotations(sketch_opt_val->outdir, tmpanno,
                                 (size_t)comblco_stat_one.infile_num);
    else
        remove_sketch_annotations(sketch_opt_val->outdir);
    write_sketch_qc_stats(sketch_opt_val->outdir, sample_qc_stats, sample_qc_len);
    if (sketch_opt_val->compute_meta)
        write_sketch_infile_meta_stats(sketch_opt_val->outdir, sample_infile_meta, sample_meta_len);

    v_free(&sketch_index);
    free(tmpname);
    free(tmpanno);
    free(sample_qc_stats);
    free(sample_infile_meta);
}
