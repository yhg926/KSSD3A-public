#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BIN=${1:-"$ROOT/bin/kssd3a"}
OUT=${2:-$(mktemp -d "${TMPDIR:-/tmp}/kssd3a-tutorial.XXXXXX")}
if [[ -e "$OUT" && ! -d "$OUT" ]]; then
    printf 'Output must be a directory: %s\n' "$OUT" >&2
    exit 1
fi
mkdir -p "$OUT"
if [[ -n "$(find "$OUT" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
    printf 'Choose an empty output directory: %s\n' "$OUT" >&2
    exit 1
fi
DATA=$ROOT/examples/data
"$BIN" --version
# f0 retains enough contexts in the deliberately tiny artificial genomes.
"$BIN" sketch -f0 -p2 -o "$OUT/assemblies" "$DATA/reference.fa" "$DATA/query.fa"
"$BIN" ani -r "$OUT/assemblies" -q "$OUT/assemblies" --metric p_dist \
    --anicut 0 --afcut 0 -o "$OUT/pairs.tsv"
"$BIN" ani -q "$OUT/assemblies" --metric p_dist --format matrix --values ani \
    --anicut 0 --afcut 0 -o "$OUT/similarity.tsv"
"$BIN" matrix --format full --metric p_dist -o "$OUT/distance.tsv" "$OUT/assemblies"
# Both mates belong to ONE biological sample, so combine them with --asone.
"$BIN" sketch --asone --conflict -A -f0 -p2 -o "$OUT/reads" \
    "$DATA/reads_R1.fq" "$DATA/reads_R2.fq"
"$BIN" ani -r "$OUT/assemblies" --qraw "$OUT/reads" --anicut 0 --afcut 0 \
    -o "$OUT/reads.tsv"
"$BIN" matrix --format dedup-plan --metric p_dist --cut 0.001 \
    --keep-out "$OUT/keep.txt" -o "$OUT/plan.tsv" "$OUT/assemblies"
"$BIN" sketch --keep "$OUT/keep.txt" -o "$OUT/representatives" "$OUT/assemblies"
python3 "$ROOT/examples/check_tutorial.py" "$OUT"
printf 'Tutorial results: %s\n' "$OUT"
