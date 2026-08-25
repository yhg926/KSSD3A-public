# KSSD3A ANI Execution Strategy

This note records the intended ANI dispatch rules for the native command-line
tool, especially the small-query/ref-stream path used for large references.

## Design Scope

KSSD3A is intended as a unified ANI/distance estimation framework for both
assembled and unassembled nucleotide sequence comparisons. The ANI layer should
support the main comparison shapes:

- assembly-to-assembly;
- assembly-to-reads;
- reads-to-assembly;
- reads-to-reads.

Input mode selects biological/statistical semantics, while execution strategy
selects a physical algorithm. These should stay separate: `-q` versus `--qraw`
sets normal assembled-style versus raw-read/unassembled ANI semantics; the
stream-versus-index decision should follow query/reference size and explicit
force/cache policy.

## Input Boundary

The ANI kernels consume sketch directories only.

Direct sequence inputs such as:

```text
kssd3a ani -r ref_sketch -q query.fastq.gz
kssd3a ani -r ref_sketch --qraw reads.fastq.gz
```

are input adapters. They sketch the sequence into a temporary query sketch,
then call the same sketch-vs-sketch ANI dispatcher as normal sketch inputs.
The core ANI code should not need to know whether a query originally came from
FASTA, FASTQ, gzip FASTQ, stdin, or a pre-existing sketch.

## Sketch Loading Boundary

ANI, matrix, distance, reverse, and inverted-index construction should load
only the sketch payload needed for comparison:

- `comblco` and `comblco.index` for direct sketch scans;
- `sortedcomb_ctxgid64obj32` when the sorted index path is selected;
- small sample-level sidecars such as annotations and input metadata when they
  are needed for reporting or best-ANI guard logic.

They should not load large per-entry sidecars such as `comblco.position` or
`comblco.a`/abundance unless the selected algorithm explicitly uses positions
or abundance. Use `generic_sketch_parse()` and request only the
needed sidecars. Normal dense comparison paths should use `SKETCH_PARSE_NONE`;
matrix graph reports request `SKETCH_PARSE_INFILE_META` for representative
selection; composition and sketch QC request abundance only when they actually
consume abundance.

## Swap Rule

For one or a few query samples against a larger reference sketch,
prefer the small-query/ref-stream path:

```text
if output is detail mode
and reference has more samples than query
and KSSD3A_FORCE_REF_INDEX is not set:
    keep the small query lookup in memory
    stream the reference comblco
```

This is the ref/query swap rule. It preserves output orientation: reported
`Qry` and `Ref` stay the user-visible query and reference, even though the
implementation scan streams the reference side.

The sorted reference index should not bypass this rule just because the index
file exists. A large cold index can be slower than streaming and can consume
much more memory/page cache.

`KSSD3A_FORCE_REF_INDEX=1` is the explicit override for using the sorted
reference index when the index file exists.

This rule is a native `kssd3a ani` execution rule, not a promise that callers
should physically reverse their command-line reference and query arguments. The
native rule keeps output oriented as the user's `Qry` and `Ref`.

## Dispatch Order

After any direct sequence input has been converted to a temporary sketch,
execution mode is selected from query/reference size and index policy:

```text
if detail output
and query sample count is smaller than reference sample count
and KSSD3A_FORCE_REF_INDEX is not set:
    use the ref/query swap rule:
        query lookup in memory, stream reference comblco
else if a sorted reference index file exists:
    use sorted reference index
else if normal assembled detail ANI should auto-build an index:
    build sorted reference index and use it
else if output is ANI full matrix mode (-m1):
    if matrix is large:
        build/use sorted reference index
    else:
        use common-index sketch-vs-sketch matrix fallback
else if output is ANI triangle mode (-m2):
    if self matrix is large:
        build/use sorted index and write lower triangle
    else:
        use one-sketch self-matrix fallback
else:
    use the pairwise sketch scan
```

For small one-sketch self output (`-q sketch -m1`, `-q sketch -m2`, or
identical `-r/-q` in matrix mode), use the self-matrix path even when a sorted
reference index exists. This keeps small full output symmetric by mirroring the
canonical lower-triangle orientation and avoids treating the same sketch as a
rectangular reference-query job.

For large `ani -m1`, KSSD3A may build and use `sortedcomb_ctxgid64obj32` even
for one-sketch or identical `-r/-q` input. This is a physical execution choice:
the output remains query rows by reference columns, but the implementation avoids
the old per-cell common-index scan. For large `ani -m2`, KSSD3A uses the same
sorted-index counting engine but only scans lower-triangle pairs and writes a
one-sketch triangle report.

The important user-facing consequence is:

```text
-q query.fastq.gz      -> normal sketch ANI after temporary sketching
--qraw query.fastq.gz  -> raw/unassembled ANI after temporary sketching
both                   -> same stream/index dispatch rule after sketching
```

So the stream-versus-index choice does not directly depend on the file
extension or on `-q` versus `--qraw`. It depends on query/reference size and
whether the sorted reference index is explicitly forced. `-q` versus `--qraw`
sets ANI semantics: normal calibrated sketch ANI versus raw/unassembled ANI.

## Direct `-q` Versus `--qraw`

Direct FASTQ query inputs are only unassembled query inputs when `--qraw` is
used. Plain `-q` is a sketch-mode adapter.

These options define ANI semantics, not the physical stream/index execution
strategy.

## ANI Matrix Output

`ani -m1` and `ani -m2` are matrix-style output modes for the same ANI
reference/query comparison. They do not use the ref/query swap rule because
matrix output needs a stable query-row by reference-column table.

Matrix-mode dispatch is:

```text
if output is -m1 and KSSD3A_FORCE_REF_INDEX is set:
    require and use sorted reference index
else if output is -m1 and max(ref_samples, query_samples) > KSSD3A_ANI_MATRIX_DIRECT_THRESHOLD:
    build sorted reference index if missing, then use it
else if output is -m1:
    use common-index sketch-vs-sketch matrix fallback
else if output is -m2 and KSSD3A_FORCE_REF_INDEX is set:
    require and use sorted reference index for lower-triangle counting
else if output is -m2 and sample_count > KSSD3A_ANI_MATRIX_DIRECT_THRESHOLD:
    build sorted index in memory if missing, then use lower-triangle counting
else if output is -m2:
    use one-sketch self-matrix fallback
else:
    use detail-mode dispatch
```

`KSSD3A_ANI_MATRIX_DIRECT_THRESHOLD` defaults to `1000`. Setting it to `0`
forces the indexed route for all ANI full and triangle matrices. The old
`KSSD3A_ANI_M1_DIRECT_THRESHOLD` name is still accepted as a compatibility
alias when the new name is not set. The threshold is independent of
`KSSD3A_AUTO_REF_INDEX_THRESHOLD`, which controls normal detail output.

## Direct Input Semantics

`-q query.fastq.gz`:

- sketches the query sequence before ANI;
- does not keep conflicting query context-objects unless `--conflict` is also
  set;
- runs normal sketch ANI, the same as running `kssd3a sketch` first and then
  `kssd3a ani -q query_sketch`;
- is intended for assembled genome sketches/sequences by default;
- does not use naive distance unless `--naive` is explicitly set.

`--qraw query.fastq.gz`:

- sketches the query sequence before ANI;
- keeps conflicting query context-objects by default;
- uses naive distance by default, equivalent to `--naive`;
- uses unassembled/raw ANI defaults;
- is intended for raw reads, metagenomic reads, and other unassembled query
  inputs.

Naive distance is zero-anchored in the low-distance region. KSSD3A first
computes the raw context-object naive distance and then applies the empirical
broad-range scale plus an intercept ramp: the intercept weight is
`min(1, raw_naive / 0.01)`. This avoids an artificial positive distance floor
for identical or near-identical raw-read comparisons while preserving the
original calibration once raw naive distance is at least `0.01`.

## Reference Index Cache Discovery

Benchmark on 2026-06-13 with public_ref and one Salmonella FASTQ-derived query
sketch, identical TSV output in all rows below:

| Mode | Cache state | Time | Max RSS |
| --- | --- | ---: | ---: |
| Streaming, no sorted ref index | cold-ish | `2:52.39` | `175 MB` |
| Streaming, no sorted ref index | warm | `15.07s` | `175 MB` |
| Sorted ref index | cold-ish | `4:19.74` | `33.3 GB` |
| Sorted ref index after stream warmup | still cold for index | `4:18.90` | `33.3 GB` |
| Sorted ref index immediately after indexed run | warm | `7.75s` | `33.3 GB` |

The stream path warms `comblco`. It does not warm
`sortedcomb_ctxgid64obj32`. Therefore a previous streaming run does not make
the sorted-index path warm.

Practical interpretation:

- Warm sorted index is fastest.
- Warm streaming is slower but still acceptable and low memory.
- Cold streaming can beat cold sorted index.
- Cold sorted index is the bad case for large-reference responsiveness.

Future adaptive optimization should check page-cache residency of
`sortedcomb_ctxgid64obj32` with Linux `mincore()` before choosing the sorted
reference index automatically. Until then, the small unassembled query stream
rule is the safer default.

Additional large-reference benchmark on 2026-06-14 with `public_ref` and a
three-sample FASTA query sketch (`165 KiB`, 3 query samples, 200,709 reference
samples):

| Mode | Cache state | Time | Max RSS | Notes |
| --- | --- | ---: | ---: | --- |
| Native detail path | cold-ish | `2:12` | `176 MB` | Small-query/ref-stream scan. |
| Native detail path | warm | `28s` | `176 MB` | Same code path after OS page cache warmed. |
| Forced reference index | cold-ish | `1:33` | `33.3 GB` | Faster here but high memory/page-cache cost. |
| Physical reverse raw path | cold-ish | `1:11` | `3.0 GB` | Drops calibrated/best output semantics. |

Interpretation:

- A first detail run against a large public reference may be slow if the
  reference sketch is cold in page cache.
- Physically reversing query/reference arguments is not a clean speed fix
  because it changes query-oriented calibrated/best and guard semantics.
- For repeated large-reference analyses, cache or reuse result files at the
  workflow level instead of reversing the biological query/reference roles.
