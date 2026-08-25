# KSSD3A Assembly ANI Guard Rules

This note records the current practical assembly ANI guard and the input metadata sidecar used to make it cheap.

## Input Metadata Sidecar

Sketches may contain `lcofiles.infilemeta`, one `infile_meta_t` record per input in `lcofiles.stat` order.

Schema version 1:

```c
typedef struct infile_meta {
  uint64_t total_length_bp;
  uint32_t record_count;
  uint32_t median_length_bp;
  float asm_level;
  float length_cv;
  uint8_t meta_fmt_version;
  int8_t infile_fmt;
  int8_t create_type;
  uint8_t infile_flags;
} infile_meta_t;
```

Field meanings:

- `meta_fmt_version`: `1` for valid version-1 metadata; `0` means metadata unavailable or invalid for guard use.
- `infile_fmt`: `0` unknown/mixed, `1` FASTA, `2` FASTQ.
- `create_type`: `0` normal one sketch per file, `1` `--asone`, `2` `--splitmfa`, `3` `--pipecmd`.
- `infile_flags`: bit flags for pipe command, stdin, compression, approximate median, or mixed input format.
- `asm_level`: for FASTA, longest record length divided by total length; for FASTQ, `0`; unknown inputs default to `-1`.
- `length_cv`: standard deviation divided by mean record length.

`sketch --nocomputemeta` suppresses `lcofiles.infilemeta`.

To keep sketching cost unchanged for read inputs, version 1 only collects
length vectors for FASTA inputs. FASTQ and unknown-format inputs may still get
a sidecar record identifying format/create flags, but per-read total, median,
and CV are left zero. The assembly ANI guard ignores non-FASTA metadata.

When merging sketches:

- All inputs have metadata: merge metadata records normally.
- All inputs lack metadata: do not write `lcofiles.infilemeta`.
- Mixed metadata/no-metadata inputs: write `lcofiles.infilemeta`, but zero every record so `meta_fmt_version == 0` and downstream guards are disabled.

## Complete-Like Query Test

The current guard version treats the query as complete-like only when query metadata is version 1 FASTA and:

```text
qry_total_length_bp >= 3.5 Mb
qry_record_count < 5
qry_median_length_bp > 100 kb
qry_total_length_bp / qry_median_length_bp < 5
```

The rule is query-side only. Reference-side metadata is used when available, but older assembly reference sketches without `lcofiles.infilemeta` do not block the query-side guard.

## ANI Best Guard Version 1

The guard applies only to normal assembly-vs-assembly detail ANI output:

- disabled for `--qraw`/unassembled paths;
- disabled for manual `-v`;
- disabled when the query side lacks valid FASTA metadata.

Default selected values:

```text
best_ANI = recalibrated_ANI
best_confidence = recalibrated_confidence
```

Guard override:

```text
if qry_complete_like
and raw_ANI >= 0.958
and Calibrated_ANI >= 0.962
and 0.25 <= context_minAF < 0.50
and raw_ANI - context_exact_mean_ANI >= 0.015:
    best_ANI = context exact AAF ANI
    best_confidence = guarded_low_confidence
```

In the current selected-metric output schema, `-s 1` reports this guarded
best estimate as `ANI = best_ANI`, `Distance = 1 - best_ANI`,
`Confidence = best_confidence`, and `Selected_metric = BestDist`.

`context_minAF` is computed from KSSD3A context align fractions: `min(Qry_align_fraction, Ref_align_fraction)`.

`context_exact_mean_ANI` is the mean of context Mash ANI, context AAF ANI, exact-object Mash ANI, and exact-object AAF ANI derived from the row counts.

This version intentionally does not track matched positions. Query-side clustered match deserts can diagnose incomplete overlap, but position tracking is too expensive for the current practical guard.

Native KSSD3A preserves the user-visible `Qry` and `Ref` orientation in output even when it uses an internal small-query streaming path for speed or memory efficiency.
