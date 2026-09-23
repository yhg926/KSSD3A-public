# KSSD3A Quickstart

## Build

The tested platform is Linux with GCC, OpenMP, GNU Make, and zlib headers.
Python 3 is needed for the tutorial checks and tests, not ordinary CLI use.
On Debian/Ubuntu, install the build/test prerequisites:

```bash
sudo apt-get install build-essential zlib1g-dev python3
make -j2
bin/kssd3a --version
bin/kssd3a doctor
```

Plain `make` is the portable build choice. `make native` targets the current
CPU; do not distribute that binary to other CPU models. `make avx2` is only
for x86 CPUs with AVX2 and BMI2. After switching compilers or compiler flags,
run `make clean` before rebuilding. OpenMP and required include paths remain
enabled when Conda or a caller supplies `CFLAGS`.

## Run A Complete Example

From the repository root:

```bash
bash examples/tutorial.sh
```

This uses bundled, deterministic artificial DNA and paired reads. No downloads
or public genome data are needed. It writes to a new temporary directory,
prints that directory, and prints a `PASS` summary of the checked results.
To choose an output directory (it must be empty or absent):

```bash
bash examples/tutorial.sh bin/kssd3a /tmp/my-kssd3a-tutorial
```

The tutorial makes assembly sketches, compares two sequences, writes distance
and ANI matrices, sketches paired reads with abundance, and applies a dedup
plan to a new output sketch. It does not modify its inputs.

| Output | What to check |
| --- | --- |
| `pairs.tsv` | Four comparisons; self-distance zero, nonself distance positive. |
| `distance.tsv` | A 2-by-2 distance matrix, diagonal zero. |
| `similarity.tsv` | A 2-by-2 ANI matrix, diagonal one; complements the distances. |
| `reads.tsv` | One biological query sample, compared to two references. |
| `keep.txt`, `plan.tsv` | One retained representative at the example cutoff. |
| `representatives/` | A new sketch containing the retained representative. |

The 8,192-base query differs from the reference at four artificial positions.
This checks software behavior, not biological accuracy. It deliberately uses
`-f0` to avoid losing most of a tiny sequence to sketch subsampling. Real
genomes normally use a larger fold, such as `-f8`, and the same fold/pattern
must be used on both sides. Lower folds retain more information at greater
memory and storage cost. Do not multiply the example's sketch distance by
genome length and call it a validated SNP count.

## Use Your Own Data

Add `bin` to your current shell's PATH:

```bash
export PATH="$PWD/bin:$PATH"
```

Each assembly file becomes one reference sample:

```bash
kssd3a sketch -f8 -p4 -o references refs/*.fna
kssd3a sketch -f8 -p4 -o query assembly.fna
kssd3a ani -r references -q query --format detail -o matches.tsv
```

Paired FASTQ files from one sample must be merged into one sketch sample:

```bash
kssd3a sketch --asone --conflict -A -f8 -p4 -o reads R1.fq.gz R2.fq.gz
kssd3a ani -r references --qraw reads --anicut 0.95 -o reads_matches.tsv
```

`-A` records abundance; it is not read-quality or depth filtering.
`--qraw` chooses the raw-query behavior; it is not a FASTQ quality-control
pipeline. See the [workflow guide](workflows.md) for read filtering and
low-divergence metric choices.

## Common Surprises

- ANI and alignment fractions use 0..1, not 0..100. Use `--anicut 0.95`.
- `ani --format matrix` alone writes distances by default. Add
  `--values ani` for similarity; `matrix` always reports distances.
- `dist` is a legacy-format command. Use `ani` or `matrix` for new sketches.
- An empty result may mean low overlap or a filter removed all hits. Check
  input layout, compatible sketch parameters, and detail-output overlap.
- `sketch --dedup` without `-o` can modify a sketch in place. Prefer a reviewed
  plan and a new output directory.

Next: [choose a workflow](workflows.md), [full manual](kssd3a_user_manual.md),
or `kssd3a ani --help`.
