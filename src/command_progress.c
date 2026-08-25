#include "command_progress.h"

#include <inttypes.h>
#include <stdio.h>

static void kssd_progress_format_duration(uint64_t seconds, char *buf, size_t buf_size)
{
    const uint64_t hours = seconds / 3600;
    const uint64_t minutes = (seconds % 3600) / 60;
    const uint64_t secs = seconds % 60;
    if (hours > 9999) {
        snprintf(buf, buf_size, ">9999h");
    } else if (hours > 0) {
        snprintf(buf, buf_size, "%02" PRIu64 ":%02" PRIu64 ":%02" PRIu64,
                 hours, minutes, secs);
    } else {
        snprintf(buf, buf_size, "%02" PRIu64 ":%02" PRIu64, minutes, secs);
    }
}

kssd_progress_t kssd_progress_start(bool enabled,
                                    const char *prefix,
                                    const char *label,
                                    const char *unit,
                                    uint64_t total)
{
    kssd_progress_t progress = {
        .enabled = enabled,
        .prefix = prefix,
        .label = label,
        .unit = unit,
        .total = total,
        .started_at = time(NULL),
        .last_at = 0,
    };
    if (progress.enabled)
        fprintf(stderr, "%s: %s started; total=%" PRIu64 " %s\n",
                progress.prefix, progress.label, progress.total, progress.unit);
    return progress;
}

void kssd_progress_update(kssd_progress_t *progress, uint64_t done, bool force)
{
    if (!progress || !progress->enabled)
        return;
    if (done > progress->total)
        done = progress->total;
    const time_t now = time(NULL);
    if (!force && now <= progress->last_at)
        return;
    progress->last_at = now;

    const double pct = progress->total > 0
                           ? (100.0 * (double)done / (double)progress->total)
                           : 100.0;
    const uint64_t elapsed = now >= progress->started_at
                                 ? (uint64_t)(now - progress->started_at)
                                 : 0;
    uint64_t eta = 0;
    if (done > 0 && progress->total > done) {
        const double rate = (double)done / (double)(elapsed > 0 ? elapsed : 1);
        if (rate > 0.0)
            eta = (uint64_t)(((double)(progress->total - done) / rate) + 0.5);
    }

    char elapsed_buf[32];
    char eta_buf[32];
    kssd_progress_format_duration(elapsed, elapsed_buf, sizeof(elapsed_buf));
    kssd_progress_format_duration(eta, eta_buf, sizeof(eta_buf));
    fprintf(stderr,
            "%s: %s %" PRIu64 "/%" PRIu64 " %s (%.2f%%) elapsed %s ETA %s\n",
            progress->prefix, progress->label, done, progress->total,
            progress->unit, pct, elapsed_buf,
            progress->total > done ? eta_buf : "00:00");
}

void kssd_progress_done(kssd_progress_t *progress)
{
    if (!progress || !progress->enabled)
        return;
    kssd_progress_update(progress, progress->total, true);
    fprintf(stderr, "%s: %s complete\n", progress->prefix, progress->label);
}

void kssd_progress_cb(void *ctx, uint64_t done, uint64_t total)
{
    (void)total;
    kssd_progress_update((kssd_progress_t *)ctx, done, false);
}
