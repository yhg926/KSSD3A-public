#ifndef COMMAND_PROGRESS_H
#define COMMAND_PROGRESS_H

#include <stdbool.h>
#include <stdint.h>
#include <time.h>

typedef struct kssd_progress
{
    bool enabled;
    const char *prefix;
    const char *label;
    const char *unit;
    uint64_t total;
    time_t started_at;
    time_t last_at;
} kssd_progress_t;

kssd_progress_t kssd_progress_start(bool enabled,
                                    const char *prefix,
                                    const char *label,
                                    const char *unit,
                                    uint64_t total);
void kssd_progress_update(kssd_progress_t *progress, uint64_t done, bool force);
void kssd_progress_done(kssd_progress_t *progress);
void kssd_progress_cb(void *ctx, uint64_t done, uint64_t total);

#endif
