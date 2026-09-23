# KSSD3A

KSSD3A is a native command-line tool for sequence sketching, ANI estimation,
distance matrices, and representative selection from assemblies or reads.

## Start Here

- [Quickstart with bundled example data](docs/quickstart.md)
- [Choose a workflow](docs/workflows.md)
- [Full CLI manual](docs/kssd3a_user_manual.md)
- [Changes in this development version](CHANGELOG.md)

## Build And Try

The tested build is Linux with GCC/OpenMP, GNU Make, and zlib development
headers. Python 3 is needed for the tests and tutorial checks.

```bash
make -j2
bin/kssd3a --version
bin/kssd3a doctor
bash examples/tutorial.sh
```

The tutorial uses tiny artificial sequences included in this repository.
It checks assembly comparisons, paired reads, matrices, and non-destructive
representative selection, then prints its output directory.

Plain `make` is the portable choice. `make native` targets your current CPU;
do not distribute that binary to unknown CPUs. `make avx2` is only for x86
with AVX2/BMI2. Run `make clean` when changing compilers or build flags.

## Common Tasks

Add the binary to your current PATH:

```bash
export PATH="$PWD/bin:$PATH"
```

Compare assemblies and write an **ANI**, not distance, matrix:

```bash
kssd3a sketch -f8 -p8 -o refs refs/*.fna
kssd3a sketch -f8 -p8 -o queries queries/*.fna
kssd3a ani -r refs -q queries --format detail -o matches.tsv
kssd3a ani -r refs -q queries --format matrix --values ani -o ani_matrix.tsv
```

Search using paired reads from **one biological sample**:

```bash
kssd3a sketch --asone --conflict -A -f8 -p8 -o reads R1.fq.gz R2.fq.gz
kssd3a ani -r refs --qraw reads --anicut 0.95 -o reads_matches.tsv
```

Both sides must use compatible fold/pattern settings. `-A` preserves counts,
not read quality. ANI and overlap filters are fractions: **0.95, not 95**.

Named metrics are available in `ani`, for example `--metric p_dist` or
`--metric ctx-naive`. Existing numeric selectors still work.
`p_dist` is a low-divergence proxy, **not an exact SNP counter**.

For distance matrices, reviewed dedup plans, coverage, and experimental
placement, see the [workflow guide](docs/workflows.md). Use `ani` or `matrix`
for current sketches; `dist` expects legacy-format files.

## Help, Tests, And Provenance

```bash
kssd3a ani --help
kssd3a sketch --help
source etc/kssd3a.bash
```

Completion requires Bash 4+. The manual and completion are checked against
the built executable. `--version` reports release/source identity; `doctor`
adds compiler/build information. Keep both with exact commands and input
checksums for reproducible analyses.

Run `make test` for native smoke tests, CLI regression tests,
documentation/completion checks, and the tutorial.

## Distribution

This is the minimal public CLI distribution. It excludes the web server,
WASM, research/manuscript records, and internal maintenance scripts.
`place` is experimental; see its manual section before use.
`SOURCE_COMMIT` records the mother source commit and clean/dirty state.
`PUBLIC_EXPORT_MANIFEST.tsv` contains SHA-256 hashes of exported files.
A dirty development snapshot is not a frozen release.

## License

See [LICENSE.txt](LICENSE.txt), [NOTICE](NOTICE), and
[commercial use guidance](COMMERCIAL.md) before distribution or commercial use.
