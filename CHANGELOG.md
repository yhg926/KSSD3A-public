# Changelog

## 3.1.1-dev (Unreleased)

- Reject missing mandatory CLI arguments with exit status 64 and diagnostics
  on stderr; propagate `place` failure status through command dispatch.
  Empty sequence lists and invalid `set` inputs no longer silently succeed.
- Preserve ANI 1 / distance 0 at the full-overlap, zero-observed-mismatch
  boundary for Best/Recalibrated output. Previously the learned recalibration
  could add a small residual even to identical sketches, including at fold 0.
  Nonzero-mismatch and partial-overlap calibration is unchanged.
- Reject modern sketch directories passed to legacy `dist` with an explicit
  format explanation and directions to `ani`/`matrix`.
- Add `sketch --asone --sample-name NAME`. The default first-input label is
  retained for compatibility; naming does not alter either mate's sketch data.
- Document missing matrix values and the existing `--exception 2` option:
  distance 2 / ANI -1 are unambiguous sentinels. Numeric defaults are unchanged.

The published `v3.1.0` tag and its benchmarking records remain unchanged.

## 3.1.0 (2026-09-23)

First tagged release of the minimal public native CLI. The version follows
the existing KSSD3A version series; it is not a claim of paper validation.

### Included Capabilities

- Native sketching, ANI search, distance matrices, representative selection,
  set operations, and abundance-aware composite analysis.
- The zero-anchored low-distance `ctx-naive` calibration and the uncalibrated
  `p_dist` micro-object rate. `p_dist` is not an exact SNP counter.
- Raw-read abundance storage with `sketch -A` and CLI coverage estimation
  with `ani --estimate-coverage`.

### CLI And Packaging

- Add readable ANI selectors: `--metric`, `--format`, and `--values`.
  Existing `-s`/`-m` options and distance formulas remain unchanged.
- Reject invalid ANI/alignment-fraction thresholds, including percentage
  inputs such as `--anicut 95`.
- Report release label, source commit/state, and build details through
  `--version` and `doctor`; retain exact public export checksums.
- Consolidate the CLI manual and add a bundled, checked artificial-data
  tutorial plus task-oriented workflow guide.
- Check docs/completion against the executable, preserve required includes
  with overridden `CFLAGS`, and add automated public-export regression checks.
- Harden public export validation and refuse to overwrite a dirty Git checkout.
- Publish source archives, checksums, and build/source provenance from tested
  annotated release tags. Check default and externally supplied compiler flags.

### Scope And Limitations

- Tested release platform: Linux with GCC, OpenMP, GNU Make, and zlib.
  Other platforms are not qualified by this release.
- `place` remains experimental. Legacy `dist` does not accept current-format
  sketch directories; use `ani` or `matrix` for those inputs.
- This distribution excludes the web server, WASM, private deployment
  configuration, research data, and manuscript/experiment records.
- The separate SketchHub saved-hit coverage executable (`--coverage-from`)
  is not part of this release.
- Release checks verify software behavior and packaging, not biological
  accuracy. See Git history and `SOURCE_COMMIT` for source provenance.
