# Choose A KSSD3A Workflow

Commands below assume `kssd3a` is on PATH and input paths are your own.
For a runnable example with included inputs, use the [quickstart](quickstart.md).

| Task | Start here | Important choice |
| --- | --- | --- |
| Compare assemblies | `sketch`, then `ani -q` | Compatible fold and sketch pattern. |
| Search with raw reads | `sketch --asone`, then `ani --qraw` | One biological sample per sketch sample. |
| Compare very close genomes | `ani --metric p_dist` | An uncalibrated proxy, not an exact SNP count. |
| Make an ANI table | `ani --format matrix --values ani` | Explicitly choose similarity rather than distance. |
| Make a tree input matrix | `matrix --matrix-format phylip` | Keep the ID map for sample labels. |
| Remove near duplicates | `matrix --format dedup-plan`, then `sketch --keep` | Review first; keep original inputs. |
| Estimate raw-read coverage | `sketch -A`, then `ani --estimate-coverage` | Estimates depend on marker specificity and overlap. |
| Tree placement | `place` | Experimental; requires a compatible backbone. |

## Assemblies And ANI Matrices

```bash
kssd3a sketch -f8 -p8 -o refs refs/*.fna
kssd3a sketch -f8 -p8 -o queries queries/*.fna
kssd3a ani -r refs -q queries --metric best --format detail -o matches.tsv
kssd3a ani -r refs -q queries --metric best \
  --format matrix --values ani -o ani_matrix.tsv
```

Detail output includes both ANI and distance plus overlap information. A high
ANI with little shared sequence is not equivalent to whole-genome identity.
`--anicut 0.95` means 95% ANI; this is a filter, not a universal species
definition. In particular, do not assume one threshold fits all fungi.

## Paired Reads, ONT, And Quality Control

Two files from one biological sample:

```bash
kssd3a sketch --asone --conflict -A -f8 -p8 -o reads R1.fq.gz R2.fq.gz
kssd3a ani -r refs --qraw reads --format detail -o reads_matches.tsv
```

For multiple biological samples, sketch each pair into its own directory,
then use the manual's append/list workflows. Do not merge unrelated samples
with `--asone`. Separate files otherwise become separate samples by default;
`--separate` instead writes independent sketch directories.

An explicit abundance-filtered ONT sketch example:

```bash
kssd3a sketch --readsQC -n4 --conflict -A -f4 -p8 -o ont ont.fastq.gz
kssd3a sketch -f4 -p8 -o refs_f4 refs/*.fna
kssd3a ani -r refs_f4 --qraw ont --format detail -o ont_matches.tsv
```

`--readsQC -n4` applies the sketcher's abundance threshold. It does not
guarantee whole-genome coverage, remove every sequencing error, or replace
base-quality assessment. `-A` stores counts without enabling this filter.
Sketching without `--readsQC` preserves more low-count evidence and noise.
Choose filtering based on data quality; do not silently compare incompatible
folds or patterns. See the manual for later `--sketchQC` processing.

## Low-Divergence Distance

```bash
kssd3a ani -r refs -q queries --metric p_dist --format detail -o close.tsv
kssd3a matrix --metric p_dist --format full -o close_matrix.tsv refs
```

`p_dist = N_diff_obj_section / (XnY_ctx * O)`, where `O` is the number of
object positions actually encoded by the sketch (`Bitslen.obj / 2`). It is
uncalibrated and intended for low-divergence comparisons with adequate shared
context. Sequencing errors, missing sequence, repeats, and sparse sampling
can all affect it. It is not a validated exact SNP counter. For raw queries,
use `--qraw` and consult the raw-query guards before forcing other metrics.

`ctx-naive` uses the calibrated context-object distance, including its
zero-anchored low-distance correction. `best`, `recalibrated`, and `ctx-moe`
retain their existing model/guard behavior. No one metric is asserted best
for all datasets; compare against suitable truth for your intended task.

## Review Duplicates Before Applying

```bash
kssd3a matrix --format dedup-plan --metric p_dist --cut 0.0001 \
  --keep-out keep.txt --remove-out remove.txt -o plan.tsv refs
kssd3a sketch --keep keep.txt -o representative_refs refs
kssd3a matrix --format full --metric p_dist --matrix-format phylip \
  --matrix-idmap ids.tsv -o representatives.phy representative_refs
```

Inspect `plan.tsv` before the second command. The cutoff is a distance, not
ANI. Deduplication is metric/overlap-based, not taxonomy-normalized species
selection. A cutoff alone does not guarantee one representative per species.
For very large references, see indexed candidate nomination and its
speed/recall tradeoffs in the manual. Dense matrices grow quadratically.

## Coverage From Abundance

```bash
kssd3a sketch --asone --conflict -A -f8 -p8 -o counts R1.fq.gz R2.fq.gz
kssd3a sketch -i refs
kssd3a ani -r refs --qraw counts --estimate-coverage -o coverage.tsv
```

Keep the `.a` abundance sidecars. Existing sketches lacking counts cannot
recover read abundance. Coverage is an inference, not known truth: divergent
queries, overlapping species markers, and mixtures can bias it. Consult the
manual's coverage columns and uniqueness-index behavior before interpreting
relative abundances. This workflow does not claim parity with any benchmark
tool without a corresponding evaluation.

## Status And Reproducibility

The routinely tested path is the native `sketch`/`ani`/`matrix`/`set` CLI.
`place` is experimental. `dist`, `shuffle`, and `reverse` are legacy or
advanced commands, not the starting point for current sketches. The minimal
public distribution excludes the web server, WASM, and research records.

Keep `kssd3a --version`, `kssd3a doctor`, exact commands, input checksums,
sketch parameters, and output paths with each run. Public source exports
also include `SOURCE_COMMIT` and a SHA-256 file manifest. Development builds
must not be represented as frozen releases.
