# KSSD3A User Manual

This manual is an initial standalone user guide for `kssd3a`. It is written as
an independent Markdown file so it can later be copied to a GitHub Wiki, a
documentation site, or packaged as online help.

For exact option spelling in a specific build, the built-in help remains the
most direct reference:

```bash
kssd3a --help
kssd3a sketch --help
kssd3a ani --help
kssd3a set --help
```

## Contents

- [1. Overview](#1-overview)
- [2. Installation And Build](#2-installation-and-build)
- [3. Main Commands](#3-main-commands)
- [4. Core Concepts](#4-core-concepts)
- [5. Sketching](#5-sketching)
- [6. ANI Workflows](#6-ani-workflows)
- [7. ANI Detail Output Columns](#7-ani-detail-output-columns)
- [8. Matrix Command](#8-matrix-command)
- [9. Set Command](#9-set-command)
- [10. Environment Variables](#10-environment-variables)
- [11. Practical Recipes](#11-practical-recipes)
- [12. Troubleshooting](#12-troubleshooting)

## 1. Overview

KSSD3A is an alignment-free framework for nucleotide sequence sketching,
distance estimation, and ANI estimation.

This public repository is a minimal native CLI distribution. It intentionally
omits server code, browser/wasm builds, manuscript files, research notebooks,
and private maintenance scripts.

It is intended to support both assembled and unassembled sequence comparisons:

- assembly-to-assembly;
- assembly-to-reads;
- reads-to-assembly;
- reads-to-reads.

The core workflow is:

```text
FASTA/FASTQ input -> sketch -> ANI/distance analysis
```

Some commands can sketch sequence files temporarily, but the stable internal
boundary is still the sketch directory. Pre-sketching is recommended for large
or repeated analyses.

## 2. Installation And Build

Build from source:

```bash
make
```

The binary is written to:

```text
bin/kssd3a
```

For a CPU-specific optimized build:

```bash
make native
```

`make native` uses `-march=native` and is meant for the local CPU. Use plain
`make` when building a binary for unknown machines. On x86_64 systems where
AVX2 and BMI2 are known to be available, `make avx2` enables those flags
explicitly.

To check the build/runtime environment:

```bash
bin/kssd3a doctor
```

To run tests:

```bash
make test
```

The public test target runs the standalone CLI smoke test suite.

To inspect smoke-test outputs:

```bash
KSSD3_KEEP_TEST_WORK=1 make test-smoke
```

Optional install:

```bash
sudo make PREFIX=/usr/local install
```

Bash completion:

```bash
source etc/kssd3a.bash
```

or system-wide:

```bash
sudo make PREFIX=/usr/local install_completion
```

## 3. Main Commands

`kssd3a` uses subcommands:

| Command | Purpose |
| --- | --- |
| `sketch` | Create sketches from FASTA/FASTQ, modify sketches, or build indexes. |
| `ani` | Estimate average nucleotide identity. |
| `matrix` | Report pairwise sketch distance matrices, sparse edges, clusters, and dedup plans. |
| `set` | Run set operations, grouping, and marker database creation. |
| `examples` | Print common command workflows. |
| `doctor` | Check build/runtime environment basics. |

Compatibility/advanced subcommands such as `dist`, `shuffle`, and `reverse`
are retained for older workflows. For new sketches produced by this release,
use `ani` for ANI estimates and `matrix` for distance matrices, sparse edges,
clusters, and deduplication plans. The legacy `dist` command expects older
cofile-style inputs and should not be used as the normal command for current
KSSD3A sketch directories.

## 4. Core Concepts

### 4.1 Sketch Directories

A sketch directory stores sampled context-object data and metadata. Common
files include:

| File | Meaning |
| --- | --- |
| `lcofiles.stat` | Sketch parameters and sample names. |
| `comblco` | Combined sketch content. |
| `comblco.index` | Per-sample offsets into `comblco`. |
| `lcofiles.infilemeta` | Optional per-input metadata sidecar. |
| `lcofiles.anno` | Optional per-sample annotations from FASTA/FASTQ headers. |
| `comblco.position` | Optional per-entry sequence positions from `--position`. |
| `sortedcomb_ctxgid64obj32` | Optional sorted reference inverted index. |

Sketches compared against each other should use compatible sketch parameters,
especially the same `--DimRdcFold` and context/object layout.

### 4.2 Normal Query Mode And Raw Query Mode

KSSD3A separates biological/statistical mode from physical execution strategy.

`-q` means normal query mode:

- intended for assembled genome sketches/sequences by default;
- does not keep conflicting query context-objects unless `--conflict` is set;
- does not use naive distance unless `--naive` is set;
- direct sequence input is temporarily sketched and then treated like `-q
  query_sketch`.

`--qraw` means raw-read/unassembled query mode:

- intended for raw reads, metagenomic reads, and unassembled query input;
- keeps conflicting query context-objects by default for direct sequence input;
- enables naive/raw-query ANI behavior by default;
- direct sequence input is temporarily sketched and then treated like `--qraw
  query_sketch`.

The stream-versus-index execution rule is independent of `-q` versus `--qraw`.
After inputs are normalized to sketches, KSSD3A chooses physical execution
from query/reference size and explicit force policy.

### 4.3 Query/Reference Orientation

Output columns are oriented as `Qry` and `Ref`. Internal ref/query swap or
streaming optimization must not change the reported orientation.

For raw-read ANI, prefer putting reads on the query side with `--qraw`:

```bash
kssd3a ani -r assembly_ref_sketch --qraw reads.fastq.gz -o reads_vs_ref.tsv
```

If reads are on the reference side and assemblies are on the query side,
KSSD3A can still compare the sketches, but raw-query semantics such as
`--qraw` are query-side semantics.

### 4.4 Stream Versus Inverted Index

For detail ANI, the intended execution strategy is:

```text
if query sample count < reference sample count
and KSSD3A_FORCE_REF_INDEX is not set:
    keep query lookup in memory and stream reference comblco
else if a complete sorted reference index exists:
    use sorted reference index
else if normal assembled detail ANI should auto-build an index:
    build sorted reference index and use it
else:
    use pairwise sketch scan
```

Use the explicit override when you know the reference index is warm/resident
and want the indexed path:

```bash
KSSD3A_FORCE_REF_INDEX=1 kssd3a ani -r ref_sketches -q qry_sketches
```

The small-query stream path can still be slow on a cold large reference because
it streams the reference sketch from disk. A repeated run may be much faster
after the OS page cache is warm. Forcing a sorted reference index can be faster
in some cases, but it can use tens of GiB of RSS/page cache and is not the
default low-memory strategy.

The focused execution note is in [`docs/kssd3a_ani_execution_strategy.md`](docs/kssd3a_ani_execution_strategy.md).

## 5. Sketching

Basic sketching:

```bash
kssd3a sketch -f8 -p8 -o ref_sketches refs/*.fasta
kssd3a sketch -f8 -p8 -o qry_sketches queries/*.fasta
```

Use an input list for large datasets:

```bash
find /data/refs -name '*.fna.gz' | sort > refs.list
kssd3a sketch -f8 -p16 -l refs.list -o ref_sketches
```

Sketch FASTQ reads:

```bash
kssd3a sketch -f8 -p8 -o read_sketches reads/*.fastq.gz
```

Sketch raw reads while keeping conflicting context-objects:

```bash
kssd3a sketch --conflict -f8 -p8 -o raw_read_sketches reads/*.fastq.gz
```

Use stdin:

```bash
zcat sample.fastq.gz | kssd3a sketch --conflict -f8 -o sample_reads -
```

Use a converter command:

```bash
kssd3a sketch --pipecmd 'samtools fastq {}' --conflict \
  -f8 -o bam_read_sketches reads.bam
```

### 5.1 Sketch Parameters

Common options:

| Option | Meaning |
| --- | --- |
| `-f, --DimRdcFold INT` | K-mer space downsampling rate `1/2^f`; default `8`. |
| `-p, --threads INT` | Number of threads. |
| `-T, --use_coden_ctxobj` | Use default coden context-object pattern. |
| `-C, --ctxlen INT` | Manual half context length. |
| `-O, --outerobjlen INT` | Manual half outer object length. |
| `-I, --innerobjlen INT` | Manual inner object length. |

For most users, the default coden context-object mode plus `-f8` is the normal
choice. Use manual `-C/-O/-I` only when you intentionally need a custom
context/object layout.

### 5.2 Read Count Filtering And QC

FASTQ sketching can filter k-mers by count:

```bash
kssd3a sketch -f8 -n2 -o reads_n2 reads.fastq.gz
```

Useful options:

| Option | Meaning |
| --- | --- |
| `-n, --LstKmerOcrs INT` | Minimum k-mer count in FASTQ input. |
| `--npercentile P` | Weighted lower-tail count percentile among k-mers passing `-n`. |
| `--readsQC` | Infer and apply a count range for read QC. |
| `--sketchQC` | Apply stored QC ranges to an existing abundance sketch. |

Example:

```bash
kssd3a sketch --readsQC -f8 -p8 \
  -o qc_reads reads/*.fastq.gz
```

### 5.3 Abundance Sketching

Build an abundance sketch:

```bash
kssd3a sketch -A -f8 -p8 -o sample_abundance samples/*.fastq.gz
```

### 5.4 Sample Layout

Treat multiple inputs as one final sample:

```bash
kssd3a sketch --asone -f8 -o combined_sample lane1.fq.gz lane2.fq.gz
```

Treat each FASTA record in a multi-FASTA file as a sample:

```bash
kssd3a sketch --splitmfa -f8 -o contig_sketch assembly.fna
```

### 5.5 Annotations, Metadata, And Positions

Write FASTA/FASTQ header annotations:

```bash
kssd3a sketch --anno -f8 -o annotated_refs refs/*.fna
```

Suppress per-input metadata sidecar:

```bash
kssd3a sketch --nocomputemeta -f8 -o no_meta refs/*.fna
```

Track zero-based sequence-stream positions:

```bash
kssd3a sketch --position -f8 -o positioned_sketch refs/*.fna
kssd3a sketch --ppos positioned_sketch > positions.tsv
```

Positions are aligned one-for-one with `comblco` entries and are zero-based
starts in the biological sequence stream for each sample, ignoring headers,
quality lines, and line breaks.

### 5.6 Building A Reference Inverted Index

Build the sorted reference index:

```bash
kssd3a sketch -i ref_sketches
```

The index can accelerate repeated many-query comparisons when it is warm in
page cache, but it can use substantial memory. The default detail ANI strategy
may stream the reference instead when the query side is much smaller.

### 5.7 Appending And Removing Samples

Copy mode: append compatible sketch directories into a new sketch:

```bash
kssd3a sketch --append -o merged_sketches base_sketch add_sketch
```

In-place mode: append into the first sketch:

```bash
kssd3a sketch --append base_sketch add_sketch
```

Remove listed samples into a new sketch:

```bash
kssd3a sketch --remove remove_names.txt -o filtered_sketch ref_sketches
```

Remove listed samples in place:

```bash
kssd3a sketch --remove remove_names.txt ref_sketches
```

Keep only listed samples into a new sketch:

```bash
kssd3a sketch --keep keep_names.txt -o kept_sketch ref_sketches
```

Keep only listed samples in place:

```bash
kssd3a sketch --keep keep_names.txt ref_sketches
```

The remove/keep list should contain sample names as printed by:

```bash
kssd3a sketch --psmp ref_sketches
```

Deduplicate a sketch into a new copy:

```bash
kssd3a sketch --dedup 0.001 --metric ctx-moe -o dedup_ref ref_sketches
```

Compare with stricter full-linkage/clique grouping:

```bash
kssd3a sketch --dedup 0.001 --metric ctx-naive \
  --dedup-strategy full-linkage -o dedup_full ref_sketches
```

Conservative near-identical cleanup:

```bash
kssd3a sketch --dedup 0.0001 --metric ctx-naive --dedup-max-afcut 0.8 ref_sketches
```

Deduplicate in place:

```bash
kssd3a sketch --dedup 0.001 ref_sketches
```

For large combined sketches, first build the sorted index and pass
`--dedup-index`. Without `--dedup-index`, `sketch --dedup` uses exhaustive
all-vs-all pair scoring; this is fine for small sketches but can be impractical
for tens of thousands of genomes.

```bash
kssd3a sketch -i ref_sketches
kssd3a sketch --dedup 0.001 --metric ctx-naive --dedup-index \
  --drop-position -p16 -o dedup_ref ref_sketches
```

The indexed path requires a non-conflict LCO sketch and
`sortedcomb_ctxgid64obj32`, which is created by `kssd3a sketch -i`. The
`--dedup-index-max-ctx-freq`, `--dedup-index-min-votes`, and
`--dedup-index-sample-step` options tune candidate nomination; candidate pairs
are still exact-scored before removal. Large indexed dedup runs report row
progress to stderr. For auditability on large references, prefer
`matrix --format dedup-plan --keep-out keep.txt --remove-out remove.txt` and
then apply the reviewed list with `sketch --keep keep.txt`.

Build a deduplicated sketch directly from FASTA/FASTQ inputs:

```bash
kssd3a sketch --dedup 0.001 --metric ctx-moe -o dedup_sketch refs/*.fna
```

If `dedup_sketch` does not exist, or exists as an empty directory, KSSD3A
sketches the inputs, deduplicates them, and promotes the result to
`dedup_sketch`. If `dedup_sketch` is already a valid sketch, the command stops:
exact deduplication requires all candidate samples in one sketch. To add new
raw inputs to an existing deduplicated reference, sketch the new inputs,
append the old and new sketches into a temporary combined sketch, and run
`kssd3a sketch --dedup DIST combined_sketch` on that full combined sketch.

With `-o`, a single argument that is already a sketch directory keeps the older
copy-mode meaning:

```bash
kssd3a sketch --dedup 0.001 -o dedup_copy existing_sketch
```

For sketches built with `--position`, the position sidecar is often as large
as `comblco`. Use `--drop-position` with `--keep`, `--remove`, or `--dedup`
when the filtered output does not need per-entry positions:

```bash
kssd3a sketch --dedup 0.001 --drop-position -o dedup_copy existing_sketch
kssd3a sketch --keep keep.txt --drop-position -o kept_copy existing_sketch
```

`--dedup DIST` compares samples in one sketch and marks a pair as duplicate
when all dedup guards pass:

```text
selected distance < DIST
max(Qry_align_fraction, Ref_align_fraction) >= --dedup-max-afcut
XnY_ctx >= --dedup-ctxcut
```

The default guard is `--dedup-max-afcut 0.8 --dedup-ctxcut 0`. This rejects
low-overlap distance edges while still allowing contained or fragmented
duplicates to be removed. Duplicate links are resolved by
`--dedup-strategy greedy` by default: KSSD3A ranks candidate representatives by
metadata quality, keeps the best remaining representative, and removes only
its direct duplicate neighbors. Every removed sample is therefore directly
within the cutoff of the representative reported for it; transitive chains are
not collapsed unless each removed sample is directly linked to its
representative. For conservative comparisons, `--dedup-strategy full-linkage`
requires each removed sample to be within cutoff of all samples already
assigned to that representative group, producing clique-like groups that may
keep more representatives.

Supported `--metric` values are:

- `ctx-moe`: context-object MoE/linear-model distance. This is the default.
- `ctx-naive`: context-object naive distance.
- `p_dist`: uncalibrated low-divergence point-mutation proxy, computed as
  `N_diff_obj_section / (XnY_ctx * O)`, where `O = Bitslen.obj / 2`.
- `mash`: Mash-style context-overlap distance for LCO sketches.
- `aaf`: AAF-style context-overlap distance for LCO sketches.

Quoted combined metrics are also accepted by `--dedup`. `A&B` requires all
listed distances to pass the cutoff; `A|B` accepts the pair if any listed
distance passes. Because `&` and `|` are shell operators, quote the expression:

```bash
kssd3a sketch --dedup 0.001 --metric 'ctx-naive&aaf' ref_sketches
kssd3a sketch --dedup 0.001 --metric 'ctx-moe|mash' ref_sketches
```

In practical tests, `ctx-naive 0.0001` behaved like near-identical cleanup,
while `ctx-naive 0.001` removed broader same-species redundancy. `ctx-moe` is
more aggressive and should be used deliberately.

When a duplicate group has more than one sample, KSSD3A keeps the most
complete representative rather than the first sample encountered. The
representative ranking is deterministic:

1. sample with valid `lcofiles.infilemeta` and positive total length
2. FASTA before FASTQ/unknown input
3. higher assembly level (`asm_level`, currently largest FASTA record divided
   by total length)
4. higher total length
5. fewer records
6. higher median record length
7. lower record-length coefficient of variation
8. more sketch entries
9. earlier sample order

After append, remove, keep, or dedup, any stale sorted reference index is
invalidated and should be rebuilt with `kssd3a sketch -i` if needed.

## 6. ANI Workflows

### 6.1 Assembly-To-Assembly

```bash
kssd3a sketch -f8 -p8 -o ref_sketches refs/*.fasta
kssd3a sketch -f8 -p8 -o qry_sketches queries/*.fasta
kssd3a sketch -i ref_sketches

kssd3a ani -r ref_sketches -q qry_sketches \
  -m 0 -p8 -o assembly_ani.tsv
```

Direct sequence input is also supported:

```bash
kssd3a ani -r ref_sketches -q query.fasta.gz \
  -m 0 -p8 -o query_vs_ref.tsv
```

Direct `-q query.fasta.gz` is equivalent to making a temporary normal query
sketch and then running `ani -q temp_sketch`.

### 6.2 Reads-To-Assembly

Reads are usually best placed on the query side:

```bash
kssd3a sketch -f8 -p8 -o ref_sketches refs/*.fasta
kssd3a sketch --conflict -f8 -p8 -o read_qry_sketches reads/*.fastq.gz

kssd3a ani -r ref_sketches --qraw read_qry_sketches \
  -f0.2 -n0.9 -p8 -o reads_to_assembly.tsv
```

One-step direct raw query:

```bash
kssd3a ani -r ref_sketches --qraw reads.fastq.gz \
  -f0.2 -n0.9 -p8 -o reads_to_assembly.tsv
```

`--qraw` keeps query conflicts by default for direct sequence input and uses
raw/unassembled ANI semantics.

If the query sketch was built with abundance counts, detail output can append
native coverage/depth estimates:

```bash
kssd3a sketch --conflict -A -f8 -p8 -o read_qry_sketches reads/*.fastq.gz
kssd3a ani -r ref_sketches --qraw read_qry_sketches \
  --estimate-coverage -m0 -f0.2 -n0.9 -p8 -o reads_to_assembly.coverage.tsv
```

`--estimate-coverage` requires `--qraw` and detail output (`-m0`). It uses the
query `comblco.a` counts and exact context-object markers that are unique among
the reported reference rows to append coverage/depth and relative abundance
columns. If the query sketch lacks `comblco.a`, the added columns report
`NA:query_missing_comblco.a`.

### 6.3 Assembly-To-Reads

If reads are stored on the reference side and assemblies on the query side,
normal sketch comparison is possible:

```bash
kssd3a sketch -f8 -p8 -o read_ref_sketches read_samples/*.fastq.gz
kssd3a sketch -f8 -p8 -o assembly_qry_sketches assemblies/*.fna

kssd3a ani -r read_ref_sketches -q assembly_qry_sketches \
  -f0 -n0 -p8 -o assembly_to_reads.tsv
```

For raw-read ANI semantics, orient reads as query and assemblies as reference
when practical. That keeps `--qraw` semantics aligned with the unassembled
side.

### 6.4 Reads-To-Reads

```bash
kssd3a sketch -f8 -p8 -o ref_read_sketches ref_reads/*.fastq.gz
kssd3a sketch --conflict -f8 -p8 -o qry_read_sketches query_reads/*.fastq.gz

kssd3a ani -r ref_read_sketches --qraw qry_read_sketches \
  -f0.2 -n0.9 -p8 -o reads_to_reads.tsv
```

You can also use direct sequence input for a single query:

```bash
kssd3a ani -r ref_read_sketches --qraw query_reads.fastq.gz \
  -f0.2 -n0.9 -p8 -o reads_to_reads.tsv
```

### 6.5 Positional Auto ANI

When `-r`, `-q`, and `--qraw` are omitted, the first positional input is the
reference and the later inputs are queries:

```bash
kssd3a ani -f0 -n0 -o pair.tsv ref.fasta query1.fasta query2.fasta
```

`--pair` is accepted for compatibility but is optional:

```bash
kssd3a ani --pair -f0 -n0 -o pair.tsv ref.fasta query.fasta
```

### 6.6 List-Based Auto ANI

Use lists when reference or query paths are too numerous for the shell:

```bash
find refs -name '*.fna.gz' | sort > refs.list
find queries -name '*.fastq.gz' | sort > qrys.list

kssd3a ani --DimRdcFold 8 \
  --reflist refs.list --qrylist qrys.list \
  -f0 -n0 -p16 -o list_ani.tsv
```

If any input in auto ANI is already a sketch, direct sequence inputs are
sketched with that sketch's parameters. If no input is a sketch, KSSD3A uses
default sketch parameters.

### 6.7 Streamed Or Converted ANI Inputs

Use stdin for one raw sequence input:

```bash
samtools fastq reads.bam | kssd3a ani \
  --conflict -f0.2 -n0.9 -o bam_vs_ref.tsv \
  ref_sketches -
```

Use `--pipecmd` when each input path needs conversion:

```bash
kssd3a ani --conflict --pipecmd 'samtools fastq {}' \
  -f0.2 -n0.9 -o bam_vs_ref.tsv \
  ref_sketches reads.bam
```

### 6.8 Filtering And Reporting

Common ANI filters:

| Option | Meaning |
| --- | --- |
| `-f, --afcut FLOAT` | Minimum `max(Qry_align_fraction, Ref_align_fraction)`. |
| `-n, --anicut FLOAT` | Minimum ANI. |
| `-t, --ctxcut INT` | Minimum overlapped context count. |
| `-N, --top INT` | Report at most top N references per query. |

Default `--afcut` is:

```text
0.5 for normal assembled-style ANI
0.2 for unassembled/qraw ANI
```

Use permissive filters for exploratory output:

```bash
kssd3a ani -r ref_sketches -q qry_sketches \
  -f0 -n0 -t0 -o all_hits.tsv
```

Select metrics:

| `-s` value | Metric |
| --- | --- |
| `1` | Best distance: best practical ANI for assembled inputs; raw-read/unassembled mode falls back to naive distance. |
| `2` | Recalibrated distance: model-calibrated ANI for assembled inputs; raw-read/unassembled mode falls back to naive distance. |
| `3` | Context-object MoE distance for assembled inputs; raw-read/unassembled mode falls back to naive distance. |
| `4` | Naive context-object distance. |
| `5` | Mash distance from shared contexts. |
| `6` | AAF distance from shared contexts. |
| `7` | MashD_if_far: MoE distance when AF passes, otherwise MashD. |
| `8` | AafD_if_far: MoE distance when AF passes, otherwise AafD. |
| `9` | `p_dist`: uncalibrated low-divergence point-mutation proxy, `N_diff_obj_section / (XnY_ctx * O)`, where `O = Bitslen.obj / 2`. |

The naive context-object distance is zero-anchored for near-identical pairs.
Below raw naive distance `0.01`, KSSD3A ramps in the empirical intercept
linearly; at raw distance `0`, the reported naive distance is also `0`, and at
raw distance `>=0.01`, the original broad-range calibration is unchanged.

For raw-read/unassembled mode, `-s 1`, `-s 2`, `-s 3`, and `-s 4` all
report the naive context-object distance by default. Metrics `-s 5` through
`-s 9` stay available as unified raw distance choices across assembled and
unassembled inputs.

Use `--unified-metric` when you intentionally want unassembled/qraw mode to
honor the requested `-s` metric instead of forcing `-s 1..4` to naive. In this
mode `-s 3` reports `CtxMoE` and `-s 4` reports `Naive`; `Best` and
`Recalibrated` still fall back when calibrated/best estimates are unavailable,
so `-s 1` and `-s 2` usually become `CtxMoE` for qraw input.

For detail output, KSSD3A always prints both `ANI` and `Distance` for the selected metric. For matrix or triangle output, positive `-s` values print distance and negative `-s` values print ANI/similarity. For example, `-s 1 -m1` prints a best-distance matrix, while `-s -1 -m1` prints the corresponding ANI matrix.

### 6.9 Output Formats

Detail output is default:

```bash
kssd3a ani -r ref_sketches -q qry_sketches -m0 -o ani_detail.tsv
```

Distance matrix output:

```bash
kssd3a ani -r ref_sketches -q qry_sketches -m1 -o distance_matrix.tsv
```

ANI matrix output:

```bash
kssd3a ani -r ref_sketches -q qry_sketches -m1 -s -1 -o ani_matrix.tsv
```

One-sketch self ANI matrix:

```bash
kssd3a ani -q qry_sketches -m1 -s -1 -o self_ani_matrix.tsv
```

For `ani -m1` and `ani -m2`, small sketch matrices use the common-index
sketch-vs-sketch comparison. Large matrices use the sorted-index counting engine
because the old direct per-cell matrix path is too slow at large sample counts.
`KSSD3A_ANI_MATRIX_DIRECT_THRESHOLD` controls the switch point and defaults to
`1000` samples. Set it to `0` to force the indexed matrix route.

When reference and query are the same small sketch, KSSD3A uses the self-matrix
path so the full matrix mirrors one canonical lower-triangle orientation. For
large same-sketch full ANI matrices, KSSD3A may use the indexed query-by-ref
path for speed. Large same-sketch triangle ANI matrices use indexed lower-pair
counting and write only the lower triangle. Use `matrix --format full` when you
need a strict symmetric sketch-distance matrix or PHYLIP output.

Triangle output:

```bash
kssd3a ani -q qry_sketches -m2 -s -1 -d -o ani_triangle.tsv
```

Skip calibrated/best ANI computation:

```bash
kssd3a ani -r ref_sketches -q qry_sketches \
  --raw-output -o raw_counts.tsv
```

## 7. ANI Detail Output Columns

Default detail output is tab-separated. Current columns are:

| Column | Meaning |
| --- | --- |
| `Qry` | Query sample name. |
| `Ref` | Reference sample name. |
| `ANI` | `1 - Distance` for the selected metric. |
| `Distance` | Distance for the selected metric. |
| `Confidence` | Confidence/status for the selected metric; `unassembled` for raw-read/unassembled mode and `raw` for raw distance metrics. |
| `Selected_metric` | Effective metric selected by `-s`, after fallback when a calibrated/best estimate is unavailable. |
| `XnY_ctx` | Shared/overlapped context count. |
| `Qry_align_fraction` | Query-side context alignment fraction. |
| `blastn_Qry_align_fraction` | BLAST-like query AF estimate, clamped to `<= 1`. |
| `Ref_align_fraction` | Reference-side context alignment fraction. |
| `blastn_Ref_align_fraction` | BLAST-like reference AF estimate, clamped to `<= 1`. |
| `N_diff_obj` | Contexts with object differences. |
| `N_diff_obj_section` | Section-level object difference count. |
| `N_mut2_ctx` | Contexts with more than one mutation section. |
| `Ref_annotation` | Reference annotation, or `NA`. |

The AF filter uses:

```text
max(Qry_align_fraction, Ref_align_fraction)
```

The focused guard-rule note is in:

```text
docs/kssd3a_ani_guard_rules.md
```

## 8. Matrix Command

`matrix` is the report-only command for pairwise sketch relation work. It can
write dense matrices, lower triangles, sparse edge lists, connected-component
cluster reports, and deduplication plans. It does not modify sketches.

One-sketch lower triangle:

```bash
kssd3a matrix --format triangle \
  -o distance_triangle.tsv ref_sketches
```

Rectangular matrix with query rows and reference columns:

```bash
kssd3a matrix -r ref_sketches -q qry_sketches \
  --format full --metric ctx-moe -o query_by_ref.tsv
```

One-sketch full PHYLIP distance matrix for tree tools:

```bash
kssd3a matrix --format full --metric ctx-naive \
  --matrix-format phylip --matrix-idmap matrix.idmap.tsv \
  -o matrix.phy ref_sketches
```

Sparse edge list:

```bash
kssd3a matrix --sparse --metric 'ctx-naive&aaf' \
  --cut 0.05 --max-afcut 0.8 --ctxcut 3 \
  -o edges.tsv ref_sketches
```

`--sparse` is an alias for `--format edges`. It prints one passing pair per
line with query, reference, distance, similarity, metric name, shared context
count, AF values, and difference counts.

Cluster report from threshold-connected components:

```bash
kssd3a matrix --format clusters --metric ctx-moe \
  --cut 0.05 -o clusters.tsv ref_sketches
```

Deduplication plan, using the same representative ranking as
`sketch --dedup`:

```bash
kssd3a matrix --format dedup-plan --metric ctx-naive \
  --cut 0.001 -o dedup_plan.tsv \
  --keep-out keep.txt --remove-out remove.txt \
  --keep-matrix-out keep_matrix.tsv ref_sketches
```

Stricter full-linkage/clique deduplication plan:

```bash
kssd3a matrix --format dedup-plan --metric ctx-naive \
  --dedup-strategy full-linkage --cut 0.001 -o dedup_plan_full.tsv \
  --keep-out keep_full.txt ref_sketches
```

PHYLIP kept-representative matrix for tree tools:

```bash
kssd3a matrix --format dedup-plan --metric ctx-naive \
  --cut 0.001 -o dedup_plan.tsv --keep-out keep.txt \
  --keep-matrix-out keep_matrix.phy --keep-matrix-format phylip \
  --keep-matrix-idmap keep_matrix.idmap.tsv ref_sketches
```

Build a distance tree from the kept-representative matrix:

```bash
rapidnj keep_matrix.phy -i pd -x keep_tree.nwk
fastme -i keep_matrix.phy -o keep_tree.nwk -T 1
```

Large indexed sketch:

```bash
kssd3a sketch -i ref_sketches
kssd3a matrix --format dedup-plan --metric ctx-naive --cut 0.001 \
  --index-max-ctx-freq 256 --index-min-votes 14 --index-sample-step 64 \
  -p8 -o dedup_plan.tsv --keep-out keep.txt \
  --keep-matrix-out keep_matrix.tsv ref_sketches
```

For cluster and dedup-plan reports, also write the reviewed edge list:

```bash
kssd3a matrix --format dedup-plan --metric ctx-naive \
  --cut 0.001 -o dedup_plan.tsv \
  --edge-out dedup_edges.tsv ref_sketches
```

Supported `--metric` values are `ctx-moe`, `ctx-naive`, `p_dist`, `mash`, and `aaf`.
The default is `ctx-naive`. For compatibility, `-m 0` means `mash` and
`-m 1` means `aaf`. For LCO sketches, `mash` and `aaf` use shared context
counts, not exact ctx-object overlap counts.

Sparse graph formats also accept quoted combined metric expressions:

```bash
kssd3a matrix --sparse --metric 'ctx-naive&aaf' --cut 0.001 ref_sketches
kssd3a matrix --format clusters --metric 'ctx-moe|mash' --cut 0.05 ref_sketches
```

`A&B` requires every listed distance to pass `--cut`; the reported edge
distance is the worst passing distance, `max(A, B, ...)`. `A|B` accepts a pair
when any listed distance passes `--cut`; the reported edge distance is the best
passing distance, `min(A, B, ...)`. Dense `full` and `triangle` matrix reports
require one metric because each cell is one numeric distance.

Sparse graph formats use this edge rule:

```text
distance < --cut
XnY_ctx >= --ctxcut
max(Qry_align_fraction, Ref_align_fraction) >= --max-afcut
```

`--format dedup-plan` defaults to `--max-afcut 0.8`, matching
`sketch --dedup`. Other sparse formats default to `--max-afcut 0`.
`--format clusters` reports connected components/single-linkage groups from
the accepted edge graph. `--format dedup-plan` uses the same default greedy
representative-centered rule as `sketch --dedup`, so component-like IDs in
that report identify representative groups and every removed sample has a
direct accepted edge to its reported representative. Pass
`--dedup-strategy full-linkage` to require clique-like duplicate groups where
each removed sample is linked to all samples already assigned to its
representative group.

When `sortedcomb_ctxgid64obj32` is present, one-sketch sparse formats
(`edges`, `clusters`, and `dedup-plan`) use the inverted index for candidate
nomination. Candidate nomination is controlled by:

```text
--index-max-ctx-freq INT   skip index context groups larger than INT; 0 disables this cap
--index-min-votes INT      exact-score only pairs with at least INT candidate votes
--index-sample-step INT    use every INT-th query context for nomination
```

These settings only affect candidate nomination. Reported edges are still
exact-scored with the selected metric and then filtered by `--cut`, `--ctxcut`,
and `--max-afcut`. Use `--index-max-ctx-freq 0 --index-min-votes 1
--index-sample-step 1` for exhaustive indexed nomination; this can still be
slow on dense same-species sketches.

`--format dedup-plan` writes the stable leading columns:

```text
sample  component_id  action  representative  reason
```

It then appends representative-quality details:

```text
rep_rule  component_size  sample_sketch_entries  rep_sketch_entries
sample_fmt  sample_asm_level  sample_total_length_bp  sample_record_count
sample_median_length_bp  sample_length_cv
rep_fmt  rep_asm_level  rep_total_length_bp  rep_record_count
rep_median_length_bp  rep_length_cv
```

`reason` describes the row role: `singleton`, `best_representative`, or
`duplicate`. `rep_rule` describes the ranking criterion that explains why the
representative wins over that row, such as `higher_asm_level`,
`longer_total_length`, `fewer_records`, `larger_sketch_entries`, or
`earlier_sample_order`. Metadata fields are `NA` when the sketch has no valid
`lcofiles.infilemeta` for that sample.

`--keep-out FILE`, `--remove-out FILE`, `--keep-matrix-out FILE`,
`--keep-matrix-format tsv|phylip`, and `--keep-matrix-idmap FILE` are
available only with `--format dedup-plan`. The first two write plain
one-sample-per-line lists without headers. The keep list is directly
consumable by `sketch --keep`:

```bash
kssd3a matrix --format dedup-plan --cut 0.001 \
  -o dedup_plan.tsv --keep-out keep.txt \
  --keep-matrix-out keep_matrix.tsv ref_sketches
kssd3a sketch --keep keep.txt -o dedup_ref ref_sketches
```

`--keep-matrix-out` writes the dense self-distance matrix for only the
representatives whose plan action is `keep`, using the same selected metric as
the dedup plan. In indexed candidate mode, the duplicate graph still uses the
index for representative selection; the final kept-representative matrix is
then scored directly so the user does not need a second sketch/dedup/matrix
workflow.

By default, `--keep-matrix-out` writes the KSSD3A TSV matrix. With
`--keep-matrix-format phylip`, it writes a PHYLIP distance matrix suitable for
distance-tree tools such as FastME and RapidNJ. PHYLIP output uses stable
10-character short IDs (`S000000001`, `S000000002`, ...), because original
sample names are often long paths and are unsafe as tree labels. Use
`--keep-matrix-idmap FILE` to write the two-column `id`/`sample` mapping.
Translate Newick labels back to original sample names with that id map after
tree construction.

The PHYLIP matrix contains distances for the selected matrix metric. It does
not contain ANI values. Standard FastME/RapidNJ workflows use dense distance
matrices, not sparse edge lists, so the PHYLIP file is `O(k^2)` in output size
for `k` kept representatives. KSSD3A writes it with one row buffer instead of
keeping the full matrix in memory, but downstream tree tools still need to load
or process the dense matrix.

For a full one-sketch matrix, `--matrix-format phylip` writes PHYLIP directly
to `-o`, and `--matrix-idmap FILE` writes the same `id`/`sample` mapping for
the full sample set. For LCO context metrics (`ctx-naive`, `ctx-moe`, `mash`,
and `aaf`), KSSD3A computes one-sketch full matrices through an indexed
lower-triangle path and then mirrors the stored distances while writing the
full matrix. This avoids the old per-cell pairwise scan for large square
matrices.

For long matrix jobs, use `--progress auto|on|off` to control progress
reporting. Progress is written only to `stderr`, and only when the main matrix
report is written with `-o`. If the matrix report is streamed to `stdout`,
KSSD3A suppresses progress even with `--progress on`, so pipes and redirected
matrix data remain clean.

`--diagonal` prints diagonal values in triangle mode. `--diagonal-value FLOAT`
sets that value. `--exception FLOAT` is the value used when a pair has no
overlap or an undefined distance.

The intended workflow is:

```bash
kssd3a matrix --format dedup-plan ... -o plan.tsv ref_sketches
# review plan.tsv
# later: kssd3a sketch --apply-dedup plan.tsv -o dedup_ref ref_sketches
```

`sketch --apply-dedup` and `sketch --split` are planned apply modes; current
`sketch --dedup` remains the immediate in-place/copy convenience command.

## 9. Set Command

`set` is for operations that derive or change sketch content across samples.
Read-only inspection lives under `sketch --p...`:

```bash
kssd3a sketch --psmp ref_sketches
kssd3a sketch --psketch ref_sketches
kssd3a sketch --pindex ref_sketches
kssd3a sketch --ppos positioned_sketches > positions.tsv
```

Union:

```bash
kssd3a set --union -o union_sketch input_sketches
```

Unique union:

```bash
kssd3a set --uniq_union -o uniq_union input_sketches
```

Marker database:

```bash
kssd3a set --uniq_union --markerdb -o markerdb input_sketches
```

Intersect or subtract with a pan sketch:

```bash
kssd3a set --intersect pan_sketch -o intersected input_sketches
kssd3a set --subtract pan_sketch -o subtracted input_sketches
```

The default set key is `--key full`, which compares the full encoded
context-object record. For `-T` long sketches, `--intersect` and `--subtract`
also support context-only matching:

```bash
kssd3a set --intersect pan_sketch --key ctx -o intersected_by_context input_sketches
kssd3a set --subtract pan_sketch --key ctx -o subtracted_by_context input_sketches
```

With `--key ctx`, only the context part is used for the membership test. The
output still stores the original full context-object records from
`input_sketches`, so the object part is preserved for downstream ANI or object
difference analysis. The older spelling `--intsect` remains available as an
alias.

For abundance sketches produced with `sketch -A`, `--intersect` and
`--subtract` also write a filtered `comblco.a` file. Each retained count remains
aligned with the retained context-object record from the original input sketch.

Group by a metadata table:

```bash
kssd3a set --grouping groups.tsv -o grouped_sketches input_sketches
```

## 10. Environment Variables

Advanced environment controls:

| Variable | Meaning |
| --- | --- |
| `KSSD3A_FORCE_REF_INDEX=1` | Force sorted reference index when the index file exists. |
| `KSSD3A_AUTO_REF_INDEX_THRESHOLD=N` | Auto-build ref index for normal detail ANI when ref sample count exceeds `N`; `0` disables. |
| `KSSD3A_ANI_MATRIX_DIRECT_THRESHOLD=N` | Keep `ani -m1`/`ani -m2` on the direct sketch-vs-sketch path up to `N` samples; above that, auto-build/use the sorted-index matrix path. Default `1000`; `0` always uses the indexed path. |
| `KSSD3A_ANI_M1_DIRECT_THRESHOLD=N` | Compatibility alias for `KSSD3A_ANI_MATRIX_DIRECT_THRESHOLD` when the new variable is unset. |
| `KSSD3A_QRAW_MULTI_THRESHOLD=N` | Maximum query sample count for multi-query small-query streaming. |
| `KSSD3A_QRAW_MULTI=sortedindex|refindex|legacy` | Override multi-query small-query strategy for unassembled/qraw mode. |
| `KSSD3A_QRAW_LOOKUP=hash|sorted` | Select lookup implementation for single-query streaming. |

Use these only when you understand the memory/speed tradeoff. In particular,
forcing a large sorted reference index can be very fast when the index is warm
but can use tens of GB of RSS/page cache.


## 11. Practical Recipes

### 11.1 Fast Search Of One Genome Against A Large Reference

```bash
kssd3a sketch -f8 -p8 -o one_genome query.fna
kssd3a ani -r public_ref -q one_genome -N10 -o top_hits.tsv
```

Force hot sorted-index path:

```bash
KSSD3A_FORCE_REF_INDEX=1 \
  kssd3a ani -r public_ref -q one_genome -N10 -o top_hits.tsv
```

### 11.2 One FASTQ Against A Large Reference

Normal sketch semantics:

```bash
kssd3a ani -r public_ref -q reads.fastq.gz -N10 -o reads_normal.tsv
```

Raw-read semantics:

```bash
kssd3a ani -r public_ref --qraw reads.fastq.gz \
  -f0.2 -n0.9 -N10 -o reads_raw.tsv
```

### 11.3 Add New Samples To A Reference Sketch

```bash
kssd3a sketch -f8 -p8 -o new_samples new/*.fna
kssd3a sketch --append ref_sketches new_samples
kssd3a sketch -i ref_sketches
```

### 11.4 Remove Samples From A Reference Sketch

```bash
kssd3a sketch --psmp ref_sketches > sample_names.tsv
cut -f2 sample_names.tsv | grep unwanted_pattern > remove.txt
kssd3a sketch --remove remove.txt ref_sketches
kssd3a sketch -i ref_sketches
```

Keep a whitelist instead:

```bash
kssd3a sketch --psmp ref_sketches > sample_names.tsv
cut -f2 sample_names.tsv | grep wanted_pattern > keep.txt
kssd3a sketch --keep keep.txt ref_sketches
kssd3a sketch -i ref_sketches
```

Deduplicate near-identical samples while keeping the most complete
representative per duplicate cluster:

```bash
kssd3a sketch --dedup 0.001 --metric ctx-moe ref_sketches
kssd3a sketch -i ref_sketches
```

Use the stricter full-linkage/clique strategy for comparison runs:

```bash
kssd3a sketch --dedup 0.001 --metric ctx-naive \
  --dedup-strategy full-linkage ref_sketches
```

For stricter near-identical cleanup:

```bash
kssd3a sketch --dedup 0.0001 --metric ctx-naive --dedup-max-afcut 0.8 ref_sketches
kssd3a sketch -i ref_sketches
```

Build a deduplicated sketch from raw references:

```bash
kssd3a sketch --dedup 0.001 --metric ctx-moe -o dedup_sketch refs/*.fna
kssd3a sketch -i dedup_sketch
```

### 11.5 Build A Human-Readable Sample Table

```bash
kssd3a sketch --psmp ref_sketches > samples.tsv
```

If the sketch was built with annotations:

```bash
kssd3a sketch --anno -f8 -o annotated refs/*.fna
kssd3a ani -r annotated -q qry_sketches -o annotated_ani.tsv
```

## 12. Troubleshooting

### 12.1 `is not a director` Or Missing `lcofiles.stat`

Older command paths treated `-q` as a sketch directory only. Current `ani`
supports direct sequence input for `-q` and `--qraw`; if you see this error,
check that you are running the intended binary:

```bash
which kssd3a
kssd3a doctor
```

### 12.2 Direct FASTQ Is Slow Compared With A Sketch

Make sure the commands have the same semantics:

```bash
kssd3a ani -r ref -q reads.fastq.gz
kssd3a sketch -o tmp_reads reads.fastq.gz
kssd3a ani -r ref -q tmp_reads
```

These should use normal sketch ANI semantics.

For raw reads:

```bash
kssd3a ani -r ref --qraw reads.fastq.gz
kssd3a sketch --conflict -o tmp_reads reads.fastq.gz
kssd3a ani -r ref --qraw tmp_reads
```

These should use raw/unassembled ANI semantics.

### 12.3 Memory Use Is High

The sorted reference index can use substantial memory. For low-memory one-query
analysis, let the default small-query stream rule run. For speed when a warm resident reference index is already available, use:

```bash
KSSD3A_FORCE_REF_INDEX=1 kssd3a ani ...
```

### 12.4 No Hits Reported

Relax filters:

```bash
kssd3a ani -r ref_sketches -q qry_sketches -f0 -n0 -t0 -o debug.tsv
```

Then inspect `XnY_ctx`, ANI, and alignment fractions.

### 12.5 Incompatible Sketches

If query and reference sketches use different sketch parameters, ANI may fail
or produce invalid comparisons. Rebuild query and reference sketches with the
same `-f` and context/object settings.
