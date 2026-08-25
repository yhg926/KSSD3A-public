# KSSD3: K-mer Space Sampling and Decomposition

KSSD3 is an alignment-free toolkit for fast comparison of genomic sequences. It uses probabilistic sampling of the k-mer space to build compact sketches and compute distances, containment metrics, and ANI estimation.

## Building from Source

The repository ships with a Makefile. To build the `kssd3a` binary run:

```bash
make
```

The build requires GCC or Clang with OpenMP support and `zlib`. When the build completes, the executable is placed in `bin/kssd3a`.

For a CPU-specific optimized build on the local machine, run:

```bash
make native
```

To print the `PATH` line for the local build output, run:

```bash
make install_env
```

To enable bash completion in the current shell:

```bash
source etc/kssd3a.bash
```

To print that completion command from the Makefile:

```bash
make completion
```

To install the completion file system-wide:

```bash
sudo make PREFIX=/usr/local install_completion
```

Optionally install the binary into `/usr/local/bin` with:

```bash
sudo make PREFIX=/usr/local install
```

## Testing

Run the Python unit tests and the main CLI smoke tests with:

```bash
make test
```

The smoke tests build small temporary sketches under `/tmp` and remove them on
success. To inspect the generated fixtures and outputs after a run:

```bash
KSSD3_KEEP_TEST_WORK=1 make test-smoke
```

## Basic Usage

Common workflows use these subcommands:

- `sketch` – create sketches from FASTA/FASTQ sequences.
- `ani` – estimate average nucleotide identity (ANI).
- `matrix` – report pairwise sketch distance matrices and graph reports.
- `set` – build union, unique-union, grouping, and marker sketches.
- `examples` – print common command workflows.
- `doctor` – check basic build/runtime environment details.

Display help for any subcommand with `--help` or `--usage`.

```bash
./bin/kssd3a sketch --help
./bin/kssd3a ani --help
./bin/kssd3a examples
./bin/kssd3a doctor
```

## Examples

The examples below show the intended command shapes. `sketch` uses the coden
context-object pattern by default, so `-T` is usually optional. Use `-p` to set
threads and keep sketch settings consistent between reference and query
sketches. Use `-C/-O/-I` only when you want to manually define the
context/object pattern instead of the coden default.

Use `-n` to require a fixed minimum k-mer count. Use `--npercentile P` to derive
an additional minimum from each input's count-weighted k-mer distribution. Use
`--readsQC` for raw reads when you want KSSD3 to infer a genomic read-count
range from the count histogram.

Use `--position` when you need `sketch` to write `comblco.position`, a
`uint64_t` sidecar aligned one-for-one with `comblco`. Coordinates are
zero-based starts in the concatenated biological sequence stream for each
sample, ignoring FASTA/FASTQ headers, quality lines, and line breaks.

### Assembled Genomes: Reference vs Query ANI

```bash
# Sketch reference genomes.
./bin/kssd3a sketch -f8 -p8 \
  -o ref_sketches refs/*.fasta

# Sketch query genomes with the same sketch parameters.
./bin/kssd3a sketch -f8 -p8 \
  -o qry_sketches queries/*.fasta

# Build the reference index for faster many-query ANI.
./bin/kssd3a sketch -i ref_sketches

# Detail output.
./bin/kssd3a ani -r ref_sketches -q qry_sketches \
  -m 0 -p8 -o ani.tsv
```

### ANI Matrix

```bash
# Matrix format from ANI.
./bin/kssd3a ani -r ref_sketches -q qry_sketches \
  -m 1 -p8 -o ani_matrix.tsv
```

### Sketch Distance Matrix And Graph Reports

```bash
# Lower-triangle sketch distance matrix among samples in one sketch.
./bin/kssd3a matrix --format triangle \
  -o distance_triangle.tsv ref_sketches

# Rectangular sketch distance matrix: query rows by reference columns.
./bin/kssd3a matrix -r ref_sketches -q qry_sketches \
  --format full --metric ctx-moe -o query_by_ref.tsv

# Full one-sketch PHYLIP distance matrix for tree tools.
./bin/kssd3a matrix --format full --metric ctx-naive \
  --matrix-format phylip --matrix-idmap matrix.idmap.tsv \
  -o matrix.phy ref_sketches

# Sparse edge list for pairs below a distance cutoff.
./bin/kssd3a matrix --format edges --metric 'ctx-naive&aaf' \
  --cut 0.05 --max-afcut 0.8 -o edges.tsv ref_sketches

# Connected-component cluster report, no sketch modification.
./bin/kssd3a matrix --format clusters --metric ctx-moe \
  --cut 0.05 -o clusters.tsv ref_sketches

# Reviewable deduplication plan using the same representative rule as sketch --dedup.
./bin/kssd3a matrix --format dedup-plan --metric ctx-naive \
  --cut 0.001 -o dedup_plan.tsv \
  --keep-out keep.txt --keep-matrix-out keep_matrix.tsv ref_sketches

# Stricter full-linkage/clique deduplication plan for conservative comparisons.
./bin/kssd3a matrix --format dedup-plan --metric ctx-naive \
  --dedup-strategy full-linkage --cut 0.001 -o dedup_plan_full.tsv \
  --keep-out keep_full.txt ref_sketches

# PHYLIP distance matrix for tree tools, with a short-ID map.
./bin/kssd3a matrix --format dedup-plan --metric ctx-naive \
  --cut 0.001 -o dedup_plan.tsv --keep-out keep.txt \
  --keep-matrix-out keep_matrix.phy --keep-matrix-format phylip \
  --keep-matrix-idmap keep_matrix.idmap.tsv ref_sketches

# The PHYLIP file can be passed to distance-tree tools.
rapidnj keep_matrix.phy -i pd -x keep_tree.nwk
fastme -i keep_matrix.phy -o keep_tree.nwk -T 1

# Large indexed sketch: use sampled candidate nomination, then exact-score candidates.
./bin/kssd3a sketch -i ref_sketches
./bin/kssd3a matrix --format dedup-plan --metric ctx-naive --cut 0.001 \
  --index-max-ctx-freq 256 --index-min-votes 14 --index-sample-step 64 \
  -p8 -o dedup_plan.tsv --keep-out keep.txt \
  --keep-matrix-out keep_matrix.tsv ref_sketches
```

### Direct ANI From FASTA/FASTQ Or Sketch Inputs

When `-r`, `-q`, and `--qraw` are omitted, the first positional input is the
reference and every later input is compared as a query. `--pair` is still
accepted for compatibility but is no longer required.
With `-r`, `-q` can also point directly to a FASTA/FASTQ file, including
`.gz` inputs. Direct `-q` sequence inputs are first sketched into a temporary
sketch and then run through the same normal sketch ANI path as `-q
query_sketch`. Conflicting context-objects are removed unless `--conflict` is
set. Use `--qraw` when the query should use raw-read/unassembled ANI.

```bash
./bin/kssd3a ani -f0 -n0 -o pair_ani.tsv \
  ref.fasta query1.fasta query2.fasta

./bin/kssd3a ani -r ref_sketches --qraw reads.fastq.gz \
  -f0.2 -n0.9 -o raw_read_ani.tsv
```

### Streamed or Tool-Converted Inputs

Use `-` as one FASTA/FASTQ input from stdin. Use `--pipecmd CMD` when each input
path must be converted by another tool; `{}` is replaced by the quoted input
path, otherwise the path is appended to the command.

```bash
samtools fastq reads.bam | ./bin/kssd3a sketch \
  --conflict -o reads_sketch -

samtools fastq reads.bam | ./bin/kssd3a ani \
  --conflict ref_read_sketch - -o raw_read_ani.tsv

./bin/kssd3a ani --conflict --pipecmd 'samtools fastq {}' \
  ref_read_sketch reads.bam -o raw_read_ani.tsv
```

### Raw-Read Query ANI

For unassembled raw reads, use `ani --qraw` with either a raw-read query sketch
or a FASTA/FASTQ sequence file, including `.gz` compressed inputs. Sequence
files passed to `--qraw` are sketched temporarily with conflicting
context-objects kept. Reference sketches are
usually built without `--conflict`.
Conflict-containing reference sketches are also supported when you intentionally need to keep reference-side conflicting context-objects.
Use `--ignoreconflict` with `ani` to skip reference-side conflict contexts during comparison.

```bash
# Reference genomes: conflict-free sketches.
./bin/kssd3a sketch -f8 -p8 \
  -o ref_sketches refs/*.fasta

# Raw-read queries: keep conflicting context-objects.
./bin/kssd3a sketch --conflict -f8 -p8 \
  -o raw_qry_sketches reads/*.fastq.gz

./bin/kssd3a sketch -i ref_sketches

./bin/kssd3a ani -r ref_sketches --qraw raw_qry_sketches \
  -m 0 -p8 -o raw_read_ani.tsv

# One-step raw query: sketch reads.fastq.gz temporarily, then run qraw ANI.
./bin/kssd3a ani -r ref_sketches --qraw reads.fastq.gz \
  -m 0 -p8 -o raw_read_ani.tsv
```

### Reads-to-Reads ANI

For reads-to-reads ANI, sketch reference reads without `--conflict`. Query reads can be sketched either without `--conflict` or with `--conflict`; use `ani --qraw` in both cases because the query is unassembled reads.

```bash
# Reference reads: conflict-free sketches.
./bin/kssd3a sketch -f8 -p8 \
  -o ref_read_sketches ref_reads/*.fastq.gz

# Query reads, conflict-free.
./bin/kssd3a sketch -f8 -p8 \
  -o qry_read_sketches query_reads/*.fastq.gz

# Alternative query reads: keep conflicting context-objects.
./bin/kssd3a sketch --conflict -f8 -p8 \
  -o qry_read_sketches_conflict query_reads/*.fastq.gz

./bin/kssd3a sketch -i ref_read_sketches

./bin/kssd3a ani -r ref_read_sketches --qraw qry_read_sketches \
  -m 0 -p8 -o reads_to_reads_ani.tsv

./bin/kssd3a ani -r ref_read_sketches --qraw qry_read_sketches_conflict \
  -m 0 -p8 -o reads_to_reads_conflict_ani.tsv
```

### Using Input Lists

For large datasets, pass a file list instead of expanding many paths in the shell.

```bash
find /data/refs -name '*.fna' | sort > refs.txt
find /data/queries -name '*.fna' | sort > queries.txt

./bin/kssd3a sketch -f8 -p16 \
  -l refs.txt -o ref_sketches

./bin/kssd3a sketch -f8 -p16 \
  -l queries.txt -o qry_sketches
```

### Group or Modify Sketch Sets

`set` can also group sketches by a metadata table or compute set operations over sketch directories.

```bash
# Print sample names stored in a sketch directory.
./bin/kssd3a sketch --psmp ref_sketches

# Print sketch positions from a sketch built with --position.
./bin/kssd3a sketch --ppos positioned_sketches

# Group genomes by a two-column metadata file.
./bin/kssd3a set -g groups.tsv -o grouped_sketches ref_sketches

# Copy mode: append compatible sketch directories into a new output sketch.
./bin/kssd3a sketch --append -o merged_sketches base_sketch new_sketch

# In-place mode: append into the first sketch argument.
./bin/kssd3a sketch --append existing_ref_sketch new_sketch

# Copy mode: remove listed samples into a new output sketch.
./bin/kssd3a sketch --remove remove_names.txt -o filtered_sketch ref_sketch

# In-place mode: remove listed samples from the first sketch argument.
./bin/kssd3a sketch --remove remove_names.txt existing_ref_sketch

# Copy mode: keep only listed samples in a new output sketch.
./bin/kssd3a sketch --keep keep_names.txt -o kept_sketch ref_sketch

# In-place mode: keep only listed samples in the first sketch argument.
./bin/kssd3a sketch --keep keep_names.txt existing_ref_sketch

# Copy mode: deduplicate and keep the most complete representative per cluster.
./bin/kssd3a sketch --dedup 0.001 --metric ctx-moe \
  -o dedup_ref_sketch ref_sketch

# Stricter full-linkage/clique grouping; useful for comparing dedup rules.
./bin/kssd3a sketch --dedup 0.001 --metric ctx-naive \
  --dedup-strategy full-linkage -o dedup_full_sketch ref_sketch

# Conservative near-identical cleanup.
./bin/kssd3a sketch --dedup 0.0001 --metric ctx-naive \
  --dedup-max-afcut 0.8 existing_ref_sketch

# In-place mode: deduplicate the first sketch argument.
./bin/kssd3a sketch --dedup 0.001 existing_ref_sketch

# Large sketch: build the sorted index and use indexed candidate nomination.
# Without --dedup-index, sketch --dedup uses exhaustive all-vs-all scoring.
./bin/kssd3a sketch -i existing_ref_sketch
./bin/kssd3a sketch --dedup 0.001 --metric ctx-naive --dedup-index \
  --drop-position -p16 -o dedup_ref_sketch existing_ref_sketch

# Build a deduplicated sketch from raw FASTA/FASTQ inputs.
./bin/kssd3a sketch --dedup 0.001 --metric ctx-moe \
  -o dedup_sketch refs/*.fna

# For sketches built with --position, omit the large comblco.position sidecar
# from filtered output.
./bin/kssd3a sketch --dedup 0.001 --drop-position \
  -o dedup_ref_sketch ref_sketch
```

## Output Notes

Default `ani -m 0` output is a tab-separated detail table. The main columns are query name, reference name, shared context count, query/reference alignment fractions, BLAST-like query/reference alignment fractions, mutation counters, and ANI. `ani -f/--afcut` filters on `max(Qry_align_fraction, Ref_align_fraction)`; the default is `0.5` for assembled queries and `0.2` for unassembled/qraw queries. Indexed ANI output can append one selected metric column after ANI.

`sketch --dedup DIST` marks duplicate pairs when selected distance `< DIST`,
`max(Qry_align_fraction, Ref_align_fraction) >= --dedup-max-afcut`, and
`XnY_ctx >= --dedup-ctxcut`. The default guard is
`--dedup-max-afcut 0.8 --dedup-ctxcut 0`. The default `--metric ctx-moe` uses
the context-object MoE/linear-model distance; `ctx-naive`, `mash`, and `aaf`
are also supported. Quoted combined metrics are supported for graph decisions:
`--metric 'ctx-naive&aaf'` requires both distances to pass, while
`--metric 'ctx-moe|mash'` accepts either distance. If metadata is available,
duplicate groups prefer valid
FASTA records with higher assembly level and total length before falling back
to sketch-entry count and original order. By default, deduplication uses
`--dedup-strategy greedy`: KSSD3A processes candidate representatives in
quality order and removes only their direct duplicate neighbors, so every
removed sample is directly within the cutoff of its reported representative.
`--dedup-strategy full-linkage` is stricter: each removed sample must be
within cutoff of all samples already assigned to that representative group,
producing clique-like duplicate groups that may keep more representatives.

For large combined sketches, build the sorted index with `sketch -i` and pass
`--dedup-index`. Otherwise `sketch --dedup` intentionally uses exhaustive
all-vs-all scoring, which can be too slow for tens of thousands of genomes.
The indexed path requires a non-conflict LCO sketch and exact-scores nominated
candidate pairs before removing samples. Large indexed dedup runs report row
progress to stderr. For reviewable large-reference work, use
`matrix --format dedup-plan --keep-out keep.txt --remove-out remove.txt` and
then apply the reviewed keep list with `sketch --keep keep.txt`.

With `-o`, `sketch --dedup` can also build a deduplicated sketch directly from
FASTA/FASTQ inputs when the output path does not exist or is an empty
directory. It refuses to update an existing sketch from raw inputs because
exact deduplication requires all candidate samples in one sketch, not only the
representatives retained by a previous dedup run.

`--drop-position` can be combined with `--keep`, `--remove`, or `--dedup` to
omit `comblco.position` from the filtered output. This is useful when a
position-tracking sketch is used for selection but the deduplicated reference
does not need per-entry positions.

`matrix` is the generic pairwise sketch relation reporter. `--format full`
writes a dense rectangular table with query rows and reference columns,
`--format triangle` writes one-sketch
lower-triangle output, `--format edges` writes a sparse TSV of pairs passing
`--cut`, `--ctxcut`, and `--max-afcut`, and `--format clusters` /
`--format dedup-plan` write reviewable graph reports without modifying the
input sketch. `dedup-plan` defaults to `--max-afcut 0.8`, matching
`sketch --dedup`. The default matrix metric is `ctx-naive`; for LCO sketches,
`mash` and `aaf` are context-overlap distances using shared context counts, not
ctx-object exact-overlap counts. Sparse
matrix graph formats also accept quoted combined metrics: `A&B` requires all
listed distances to pass `--cut`, and `A|B` accepts the pair if any listed
distance passes. Dense `full` and `triangle` reports require one metric.
`clusters` reports connected components/single-linkage groups. `dedup-plan`
uses the same default greedy representative-centered rule as `sketch --dedup`;
pass `--dedup-strategy full-linkage` to require clique-like duplicate groups.
Component-like IDs in that report identify representative groups, not
transitive connected components. `dedup-plan` keeps the stable leading columns
`sample`, `component_id`,
`action`, `representative`, and `reason`, then appends `rep_rule` plus
sample/representative metadata so the chosen representative can be audited.
Use `--keep-out keep.txt` and `--remove-out remove.txt` to write one-column
sample lists; `keep.txt` can be passed directly to `sketch --keep`.
For a one-sketch full dense matrix, `--matrix-format phylip` writes PHYLIP
directly and `--matrix-idmap matrix.idmap.tsv` records the short-ID to sample
mapping. For LCO context metrics, one-sketch full matrices use a lower-triangle
indexed path and then mirror the stored distances when writing the full table.
Use `--keep-matrix-out keep_matrix.tsv` with `dedup-plan` to also write the
dense self-distance matrix for the kept representatives, avoiding a separate
deduplication and matrix workflow.
Use `--keep-matrix-format phylip` to write the same kept-representative matrix
in PHYLIP distance-matrix format for tools such as FastME and RapidNJ.
PHYLIP output uses short IDs (`S000000001`, `S000000002`, ...); use
`--keep-matrix-idmap keep_matrix.idmap.tsv` to preserve the mapping to original
sample names. The matrix values are distances for the selected metric, not ANI.
FastME and RapidNJ expect dense PHYLIP matrices, so large retained sets still
produce `O(k^2)` output even though KSSD3A writes the file one row at a time.
When `sortedcomb_ctxgid64obj32` is present, sparse one-sketch matrix reports
use the inverted index for candidate nomination. `--index-max-ctx-freq`,
`--index-min-votes`, and `--index-sample-step` tune that nomination stage;
nominated pairs are still exact-scored before an edge is reported.
For long matrix jobs, `--progress auto|on|off` reports row/pair progress to
`stderr` only when the main matrix report is written with `-o`; if matrix output
is streamed to `stdout`, progress is suppressed so the data stream stays clean.

`ani -m 1` and `ani -m 2` write ANI-oriented matrix-style tables. The first row
contains column sample names, and each later row starts with the row sample name
followed by numeric values. Small ANI matrix jobs use the direct
sketch-vs-sketch path; large full and triangle ANI matrices use the sorted-index
counting path. Tune that switch with `KSSD3A_ANI_MATRIX_DIRECT_THRESHOLD`
(`1000` by default, `0` to always use the indexed matrix path). The old
`KSSD3A_ANI_M1_DIRECT_THRESHOLD` name remains accepted as an alias.

## Further Documentation

The standalone user manual is in [`docs/kssd3a_user_manual.md`](docs/kssd3a_user_manual.md).

Command line options for each subcommand are provided by the built-in help system:

```bash
./bin/kssd3a <subcommand> --help
```

For license details see [`LICENSE.txt`](LICENSE.txt).
