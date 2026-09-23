# Changelog

## 3.1.0-dev (Unreleased)

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

This is a development change set, not a published/tagged release. See Git
history and each export's provenance for earlier changes.
