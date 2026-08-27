#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BIN=${KSSD3_BIN:-"$ROOT/bin/kssd3a"}
TMPBASE=${TMPDIR:-/tmp}
WORKDIR=$(mktemp -d "$TMPBASE/kssd3_smoke.XXXXXX")

cleanup() {
  if [ "${KSSD3_KEEP_TEST_WORK:-0}" = "1" ]; then
    printf 'smoke: kept workdir %s\n' "$WORKDIR" >&2
  else
    rm -rf "$WORKDIR"
  fi
}
trap cleanup EXIT

run() {
  printf 'smoke: %s\n' "$*" >&2
  "$@"
}

require_nonempty() {
  if [ ! -s "$1" ]; then
    printf 'smoke: expected non-empty file: %s\n' "$1" >&2
    exit 1
  fi
}

require_exists() {
  if [ ! -e "$1" ]; then
    printf 'smoke: expected file or directory to exist: %s\n' "$1" >&2
    exit 1
  fi
}

require_same_size() {
  local a=$1
  local b=$2
  local asize
  local bsize
  asize=$(stat -c %s "$a")
  bsize=$(stat -c %s "$b")
  if [ "$asize" != "$bsize" ]; then
    printf 'smoke: expected %s and %s to have same size, got %s and %s\n' "$a" "$b" "$asize" "$bsize" >&2
    exit 1
  fi
}

require_abundance_matches_comblco() {
  local sketch=$1
  require_exists "$sketch/comblco"
  require_exists "$sketch/comblco.a"
  local csize
  local asize
  csize=$(stat -c %s "$sketch/comblco")
  asize=$(stat -c %s "$sketch/comblco.a")
  if [ $((csize % 8)) -ne 0 ] || [ $((asize % 4)) -ne 0 ] || [ $((csize / 8)) -ne $((asize / 4)) ]; then
    printf 'smoke: expected %s/comblco and %s/comblco.a to have matching record counts, got %s and %s bytes\n' "$sketch" "$sketch" "$csize" "$asize" >&2
    exit 1
  fi
}

require_multiset_subset() {
  local subset=$1
  local superset=$2
  awk -v subset="$subset" -v superset="$superset" '
    NR == FNR { seen[$0]++; next }
    {
      if (seen[$0] > 0) {
        seen[$0]--
      } else {
        printf "smoke: expected every line in %s to be retained from %s; missing: %s\n", subset, superset, $0 > "/dev/stderr"
        exit 1
      }
    }
  ' "$superset" "$subset"
}

require_same_file() {
  local a=$1
  local b=$2
  if ! cmp -s "$a" "$b"; then
    printf 'smoke: expected %s and %s to be identical\n' "$a" "$b" >&2
    exit 1
  fi
}

require_same_tsv_except_first_col() {
  local a=$1
  local b=$2
  local aa="$WORKDIR/cmp_a.tsv"
  local bb="$WORKDIR/cmp_b.tsv"
  cut -f2- "$a" > "$aa"
  cut -f2- "$b" > "$bb"
  if ! cmp -s "$aa" "$bb"; then
    printf 'smoke: expected %s and %s to match after removing the query-name column\n' "$a" "$b" >&2
    exit 1
  fi
}

require_blastn_af_le_one() {
  local file=$1
  awk -F '\t' -v file="$file" '
    NR == 1 { next }
    (($9 != "NULL" && $9 + 0 > 1.000001) ||
     ($11 != "NULL" && $11 + 0 > 1.000001)) {
      printf "smoke: %s has blastn AF > 1 at line %d: col9=%s col11=%s\n", file, NR, $9, $11 > "/dev/stderr"
      exit 1
    }
  ' "$file"
}

require_header_has() {
  local file=$1
  local pattern=$2
  if ! head -n 1 "$file" | grep -q "$pattern"; then
    printf 'smoke: expected header in %s to contain %s\n' "$file" "$pattern" >&2
    exit 1
  fi
}

require_header_lacks() {
  local file=$1
  local pattern=$2
  if head -n 1 "$file" | grep -q "$pattern"; then
    printf 'smoke: expected header in %s to omit %s\n' "$file" "$pattern" >&2
    exit 1
  fi
}

require_file_has() {
  local file=$1
  local pattern=$2
  if ! grep -Fqx "$pattern" "$file"; then
    printf 'smoke: expected %s to contain exact line %s\n' "$file" "$pattern" >&2
    exit 1
  fi
}

require_file_contains() {
  local file=$1
  local pattern=$2
  if ! grep -Fq "$pattern" "$file"; then
    printf 'smoke: expected %s to contain %s\n' "$file" "$pattern" >&2
    exit 1
  fi
}

require_file_not_contains() {
  local file=$1
  local pattern=$2
  if grep -Fq "$pattern" "$file"; then
    printf 'smoke: expected %s not to contain %s\n' "$file" "$pattern" >&2
    exit 1
  fi
}

require_file_lacks() {
  local file=$1
  local pattern=$2
  if grep -Fqx "$pattern" "$file"; then
    printf 'smoke: expected %s to omit exact line %s\n' "$file" "$pattern" >&2
    exit 1
  fi
}

require_line_count() {
  local file=$1
  local expected=$2
  local observed
  observed=$(wc -l < "$file")
  if [ "$observed" != "$expected" ]; then
    printf 'smoke: expected %s lines in %s, got %s\n' "$expected" "$file" "$observed" >&2
    exit 1
  fi
}

write_repeated_fasta() {
  local path=$1
  local name=$2
  local motif=$3
  local repeats=$4
  {
    printf '>%s\n' "$name"
    local i
    for ((i = 0; i < repeats; i++)); do
      printf '%s' "$motif"
    done
    printf '\n'
  } > "$path"
}

write_concat_fasta() {
  local path=$1
  local name=$2
  shift 2
  {
    printf '>%s\n' "$name"
    while [ "$#" -gt 0 ]; do
      local motif=$1
      local repeats=$2
      local i
      for ((i = 0; i < repeats; i++)); do
        printf '%s' "$motif"
      done
      shift 2
    done
    printf '\n'
  } > "$path"
}

make_lcg_seq() {
  local length=$1
  local seed=$2
  local x=$seed
  local i
  for ((i = 0; i < length; i++)); do
    x=$(((x * 1103515245 + 12345) & 2147483647))
    case $(((x >> 16) & 3)) in
      0) printf 'A' ;;
      1) printf 'C' ;;
      2) printf 'G' ;;
      *) printf 'T' ;;
    esac
  done
}

write_seq_fasta() {
  local path=$1
  local name=$2
  local seq=$3
  {
    printf '>%s\n' "$name"
    printf '%s\n' "$seq"
  } > "$path"
}

write_fastq_record() {
  local name=$1
  local seq=$2
  printf '@%s\n%s\n+\n' "$name" "$seq"
  local i
  for ((i = 0; i < ${#seq}; i++)); do
    printf 'I'
  done
  printf '\n'
}

REFS="$WORKDIR/refs"
QRYS="$WORKDIR/qrys"
READS="$WORKDIR/reads"
AFREFS="$WORKDIR/afcut_refs"
AFQRYS="$WORKDIR/afcut_qrys"
OUT="$WORKDIR/out"
LISTS="$WORKDIR/lists"
mkdir -p "$REFS" "$QRYS" "$READS" "$AFREFS" "$AFQRYS" "$OUT" "$LISTS"

write_repeated_fasta "$REFS/refA.fna" refA ACGTTGCA 30
{
  printf '>refA_frag_1\n'
  for ((i = 0; i < 15; i++)); do printf 'ACGTTGCA'; done
  printf '\n>refA_frag_2\n'
  for ((i = 0; i < 15; i++)); do printf 'ACGTTGCA'; done
  printf '\n'
} > "$REFS/refA_frag.fna"
write_repeated_fasta "$REFS/refB.fna" refB GATTACAG 30
write_repeated_fasta "$QRYS/qryA.fna" qryA ACGTTGCA 30
write_repeated_fasta "$QRYS/qryA2.fna" qryA2 ACGTTGCA 30
write_concat_fasta "$QRYS/qryMix.fna" qryMix ACGTTGCA 20 GATTACAG 5
{
  write_fastq_record read1 ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA
  write_fastq_record read2 GATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAG
} > "$READS/reads.fastq"
gzip -c "$READS/reads.fastq" > "$READS/reads.fastq.gz"
AF_REF_SEQ=$(make_lcg_seq 800 17)
AF_EXTRA_SEQ=$(make_lcg_seq 2400 923)
write_seq_fasta "$AFREFS/refRefPlusExtra.fna" refRefPlusExtra "${AF_REF_SEQ}${AF_EXTRA_SEQ}"
write_seq_fasta "$AFQRYS/qrySubset.fna" qrySubset "$AF_REF_SEQ"

run "$BIN" doctor

run "$BIN" sketch -f0 -p2 -o "$OUT/ref_sketch" "$REFS/refA.fna" "$REFS/refB.fna"
run "$BIN" sketch -f0 -p2 -o "$OUT/qry_sketch" "$QRYS/qryA.fna" "$QRYS/qryA2.fna"
run "$BIN" sketch --conflict -f0 -p2 -o "$OUT/read_sketch" "$READS/reads.fastq"
run "$BIN" sketch -f0 -p2 -o "$OUT/read_sketch_noconflict" "$READS/reads.fastq"
run "$BIN" sketch --position -f0 -p2 -o "$OUT/pos_sketch_many" "$REFS/refA.fna" "$REFS/refB.fna"
run "$BIN" sketch --position --conflict -f0 -p4 -o "$OUT/pos_sketch_one" "$READS/reads.fastq"
run "$BIN" sketch -i "$OUT/ref_sketch"
require_nonempty "$OUT/ref_sketch/lcofiles.stat"
require_nonempty "$OUT/ref_sketch/sortedcomb_ctxgid64obj32"
run "$BIN" set -P "$OUT/ref_sketch" > "$OUT/ref_sketch.set_psmp"
run "$BIN" sketch --psmp "$OUT/ref_sketch" > "$OUT/ref_sketch.sketch_psmp"
require_same_file "$OUT/ref_sketch.set_psmp" "$OUT/ref_sketch.sketch_psmp"
run "$BIN" sketch --psketch "$OUT/ref_sketch" > "$OUT/ref_sketch.psketch"
require_nonempty "$OUT/ref_sketch.psketch"
run "$BIN" sketch --pindex "$OUT/ref_sketch" > "$OUT/ref_sketch.pindex"
require_nonempty "$OUT/ref_sketch.pindex"
cp -a "$OUT/ref_sketch" "$OUT/remove_sketch"
printf '%s\n' "$REFS/refB.fna" > "$OUT/remove_names.txt"
run "$BIN" sketch --remove "$OUT/remove_names.txt" -o "$OUT/remove_copy_sketch" "$OUT/ref_sketch"
run "$BIN" set -P "$OUT/remove_copy_sketch" > "$OUT/remove_copy_sketch.names"
run "$BIN" sketch --psmp "$OUT/remove_copy_sketch" > "$OUT/remove_copy_sketch.sketch_names"
require_same_file "$OUT/remove_copy_sketch.names" "$OUT/remove_copy_sketch.sketch_names"
cut -f2- "$OUT/remove_copy_sketch.names" > "$OUT/remove_copy_sketch.sample_names"
require_line_count "$OUT/remove_copy_sketch.sample_names" 1
require_file_has "$OUT/remove_copy_sketch.sample_names" "$REFS/refA.fna"
require_file_lacks "$OUT/remove_copy_sketch.sample_names" "$REFS/refB.fna"
run "$BIN" sketch --remove "$OUT/remove_names.txt" "$OUT/remove_sketch"
run "$BIN" set -P "$OUT/remove_sketch" > "$OUT/remove_sketch.names"
cut -f2- "$OUT/remove_sketch.names" > "$OUT/remove_sketch.sample_names"
require_line_count "$OUT/remove_sketch.sample_names" 1
require_file_has "$OUT/remove_sketch.sample_names" "$REFS/refA.fna"
require_file_lacks "$OUT/remove_sketch.sample_names" "$REFS/refB.fna"
if [ -e "$OUT/remove_sketch/sortedcomb_ctxgid64obj32" ]; then
  printf 'smoke: --remove should delete stale sortedcomb_ctxgid64obj32\n' >&2
  exit 1
fi
cp -a "$OUT/ref_sketch" "$OUT/keep_sketch"
printf '%s\n' "$REFS/refA.fna" > "$OUT/keep_names.txt"
run "$BIN" sketch --keep "$OUT/keep_names.txt" -o "$OUT/keep_copy_sketch" "$OUT/ref_sketch"
run "$BIN" set -P "$OUT/keep_copy_sketch" > "$OUT/keep_copy_sketch.names"
cut -f2- "$OUT/keep_copy_sketch.names" > "$OUT/keep_copy_sketch.sample_names"
require_line_count "$OUT/keep_copy_sketch.sample_names" 1
require_file_has "$OUT/keep_copy_sketch.sample_names" "$REFS/refA.fna"
require_file_lacks "$OUT/keep_copy_sketch.sample_names" "$REFS/refB.fna"
printf '%s\n' "$REFS/refB.fna" > "$OUT/keep_names_refB.txt"
run "$BIN" sketch --keep "$OUT/keep_names_refB.txt" "$OUT/keep_sketch"
run "$BIN" set -P "$OUT/keep_sketch" > "$OUT/keep_sketch.names"
cut -f2- "$OUT/keep_sketch.names" > "$OUT/keep_sketch.sample_names"
require_line_count "$OUT/keep_sketch.sample_names" 1
require_file_has "$OUT/keep_sketch.sample_names" "$REFS/refB.fna"
require_file_lacks "$OUT/keep_sketch.sample_names" "$REFS/refA.fna"
if [ -e "$OUT/keep_sketch/sortedcomb_ctxgid64obj32" ]; then
  printf 'smoke: --keep should delete stale sortedcomb_ctxgid64obj32\n' >&2
  exit 1
fi
run "$BIN" sketch -f0 -p2 -o "$OUT/dedup_src" "$REFS/refA_frag.fna" "$REFS/refA.fna" "$REFS/refB.fna"
run "$BIN" sketch -i "$OUT/dedup_src"
run "$BIN" matrix --format edges --metric ctx-naive --cut 0.001 --max-afcut 0.8 -o "$OUT/matrix_edges.tsv" "$OUT/dedup_src"
require_line_count "$OUT/matrix_edges.tsv" 2
require_file_contains "$OUT/matrix_edges.tsv" "Distance"
require_file_contains "$OUT/matrix_edges.tsv" "$REFS/refA.fna"
require_file_contains "$OUT/matrix_edges.tsv" "$REFS/refA_frag.fna"
run "$BIN" matrix --format edges --metric 'ctx-naive&aaf' --cut 0.001 --max-afcut 0.8 -o "$OUT/matrix_edges_and.tsv" "$OUT/dedup_src"
require_line_count "$OUT/matrix_edges_and.tsv" 2
require_file_contains "$OUT/matrix_edges_and.tsv" "ctx-naive&aaf"
run "$BIN" matrix --format edges --metric 'ctx-moe|mash' --cut 0.001 --max-afcut 0.8 -o "$OUT/matrix_edges_or.tsv" "$OUT/dedup_src"
require_line_count "$OUT/matrix_edges_or.tsv" 2
require_file_contains "$OUT/matrix_edges_or.tsv" "ctx-moe|mash"
run "$BIN" matrix --format edges --metric p_dist --cut 1 --max-afcut 0 -o "$OUT/matrix_edges_pdist.tsv" "$OUT/dedup_src"
require_line_count "$OUT/matrix_edges_pdist.tsv" 2
require_file_contains "$OUT/matrix_edges_pdist.tsv" "p_dist"
if run "$BIN" matrix --format triangle --metric 'ctx-naive&aaf' -o "$OUT/matrix_bad_combined.tsv" "$OUT/dedup_src" 2> "$OUT/matrix_bad_combined.err"; then
  printf 'smoke: combined matrix metric should require graph format\n' >&2
  exit 1
fi
require_file_contains "$OUT/matrix_bad_combined.err" "combined --metric expressions"
run "$BIN" matrix --format clusters --metric ctx-naive --cut 0.001 --max-afcut 0.8 -o "$OUT/matrix_clusters.tsv" "$OUT/dedup_src"
require_line_count "$OUT/matrix_clusters.tsv" 4
awk -F '\t' -v a="$REFS/refA.fna" -v frag="$REFS/refA_frag.fna" '
  NR > 1 { cid[$1] = $2; rep[$1] = $4 }
  END {
    if (cid[a] == "" || cid[frag] == "" || cid[a] != cid[frag] ||
        rep[a] != a || rep[frag] != a) exit 1
  }
' "$OUT/matrix_clusters.tsv"
run "$BIN" matrix --format dedup-plan --metric ctx-naive --cut 0.001 -o "$OUT/matrix_dedup_plan.tsv" --edge-out "$OUT/matrix_dedup_edges.tsv" --keep-out "$OUT/matrix_keep.txt" --remove-out "$OUT/matrix_remove.txt" --keep-matrix-out "$OUT/matrix_keep_dist.tsv" "$OUT/dedup_src"
run "$BIN" matrix --format dedup-plan --metric ctx-naive --dedup-strategy full-linkage --cut 0.001 -o "$OUT/matrix_dedup_plan_full.tsv" --keep-out "$OUT/matrix_keep_full.txt" --remove-out "$OUT/matrix_remove_full.txt" "$OUT/dedup_src"
run "$BIN" matrix --format dedup-plan --metric ctx-naive --cut 0.001 -o "$OUT/matrix_dedup_plan_phylip.tsv" --keep-matrix-out "$OUT/matrix_keep_dist.phy" --keep-matrix-format phylip --keep-matrix-idmap "$OUT/matrix_keep_dist.idmap.tsv" "$OUT/dedup_src"
if run "$BIN" matrix --format edges --metric ctx-naive --dedup-strategy full-linkage --cut 0.001 -o "$OUT/matrix_bad_dedup_strategy.tsv" "$OUT/dedup_src" 2> "$OUT/matrix_bad_dedup_strategy.err"; then
  printf 'smoke: matrix --dedup-strategy should require dedup-plan format\n' >&2
  exit 1
fi
require_file_contains "$OUT/matrix_bad_dedup_strategy.err" "dedup-strategy requires --format dedup-plan"
if run "$BIN" matrix --format dedup-plan --metric ctx-naive --cut 0.001 -o "$OUT/matrix_same.tsv" --edge-out "$OUT/matrix_same.tsv" "$OUT/dedup_src" 2> "$OUT/matrix_same.err"; then
  printf 'smoke: matrix --edge-out should reject --outfile path reuse\n' >&2
  exit 1
fi
require_file_contains "$OUT/matrix_same.err" "edge-out must differ from --outfile"
if run "$BIN" matrix --format dedup-plan --metric ctx-naive --cut 0.001 -o "$OUT/matrix_same_keep.tsv" --keep-out "$OUT/matrix_same_keep_list.tsv" --keep-matrix-out "$OUT/matrix_same_keep_list.tsv" "$OUT/dedup_src" 2> "$OUT/matrix_same_keep.err"; then
  printf 'smoke: matrix --keep-matrix-out should reject --keep-out path reuse\n' >&2
  exit 1
fi
require_file_contains "$OUT/matrix_same_keep.err" "keep-matrix-out must differ from --keep-out"
if run "$BIN" matrix --format dedup-plan --metric ctx-naive --cut 0.001 -o "$OUT/matrix_bad_idmap.tsv" --keep-matrix-out "$OUT/matrix_bad_idmap.keep.tsv" --keep-matrix-idmap "$OUT/matrix_bad_idmap.idmap.tsv" "$OUT/dedup_src" 2> "$OUT/matrix_bad_idmap.err"; then
  printf 'smoke: matrix --keep-matrix-idmap should require phylip format\n' >&2
  exit 1
fi
require_file_contains "$OUT/matrix_bad_idmap.err" "keep-matrix-idmap requires --keep-matrix-format phylip"
if run "$BIN" matrix --format full --metric ctx-naive --matrix-idmap "$OUT/matrix_bad_full_idmap.tsv" -o "$OUT/matrix_bad_full_idmap.out" "$OUT/dedup_src" 2> "$OUT/matrix_bad_full_idmap.err"; then
  printf 'smoke: matrix --matrix-idmap should require phylip format\n' >&2
  exit 1
fi
require_file_contains "$OUT/matrix_bad_full_idmap.err" "matrix-idmap requires --matrix-format phylip"
require_line_count "$OUT/matrix_dedup_plan.tsv" 4
require_line_count "$OUT/matrix_dedup_plan_full.tsv" 4
require_line_count "$OUT/matrix_dedup_edges.tsv" 2
require_line_count "$OUT/matrix_keep.txt" 2
require_line_count "$OUT/matrix_keep_full.txt" 2
require_line_count "$OUT/matrix_remove.txt" 1
require_line_count "$OUT/matrix_remove_full.txt" 1
require_line_count "$OUT/matrix_keep_dist.tsv" 3
require_line_count "$OUT/matrix_keep_dist.phy" 3
require_line_count "$OUT/matrix_keep_dist.idmap.tsv" 3
require_header_has "$OUT/matrix_dedup_plan.tsv" rep_rule
require_header_has "$OUT/matrix_dedup_plan.tsv" sample_asm_level
require_header_has "$OUT/matrix_dedup_plan.tsv" rep_total_length_bp
require_file_has "$OUT/matrix_keep.txt" "$REFS/refA.fna"
require_file_has "$OUT/matrix_keep.txt" "$REFS/refB.fna"
require_file_has "$OUT/matrix_remove.txt" "$REFS/refA_frag.fna"
awk 'NR == 1 { if ($1 != 2) exit 1; next } { if (NF != 3 || $1 !~ /^S[0-9]{9}$/) exit 1 }' "$OUT/matrix_keep_dist.phy"
require_file_contains "$OUT/matrix_keep_dist.idmap.tsv" "id	sample"
require_file_contains "$OUT/matrix_keep_dist.idmap.tsv" "$REFS/refA.fna"
require_file_contains "$OUT/matrix_keep_dist.idmap.tsv" "$REFS/refB.fna"
awk -F '\t' -v a="$REFS/refA.fna" -v frag="$REFS/refA_frag.fna" '
  NR > 1 { action[$1] = $3; rep[$1] = $4; rule[$1] = $6 }
  END {
    if (action[a] != "keep" || action[frag] != "remove" || rep[frag] != a) exit 1
    if (rule[frag] == "") exit 1
  }
' "$OUT/matrix_dedup_plan.tsv"
run "$BIN" sketch --keep "$OUT/matrix_keep.txt" -o "$OUT/matrix_keep_sketch" "$OUT/dedup_src"
run "$BIN" matrix --format full --metric ctx-naive -o "$OUT/matrix_keep_dist_from_sketch.tsv" "$OUT/matrix_keep_sketch"
require_same_file "$OUT/matrix_keep_dist.tsv" "$OUT/matrix_keep_dist_from_sketch.tsv"
run "$BIN" matrix --format full --metric ctx-naive --matrix-format phylip --matrix-idmap "$OUT/matrix_full_dist.idmap.tsv" -o "$OUT/matrix_full_dist.phy" "$OUT/matrix_keep_sketch"
require_line_count "$OUT/matrix_full_dist.phy" 3
require_line_count "$OUT/matrix_full_dist.idmap.tsv" 3
awk 'NR == 1 { if ($1 != 2) exit 1; next } { if (NF != 3 || $1 !~ /^S[0-9]{9}$/) exit 1 }' "$OUT/matrix_full_dist.phy"
require_file_contains "$OUT/matrix_full_dist.idmap.tsv" "id	sample"
require_file_contains "$OUT/matrix_full_dist.idmap.tsv" "$REFS/refA.fna"
require_file_contains "$OUT/matrix_full_dist.idmap.tsv" "$REFS/refB.fna"
run "$BIN" matrix --format dedup-plan --metric aaf --cut 0.001 -o "$OUT/matrix_dedup_plan_aaf.tsv" --keep-out "$OUT/matrix_keep_aaf.txt" --keep-matrix-out "$OUT/matrix_keep_aaf.tsv" "$OUT/dedup_src"
run "$BIN" sketch --keep "$OUT/matrix_keep_aaf.txt" -o "$OUT/matrix_keep_aaf_sketch" "$OUT/dedup_src"
run "$BIN" matrix --format full --metric aaf -o "$OUT/matrix_keep_aaf_from_sketch.tsv" "$OUT/matrix_keep_aaf_sketch"
require_same_file "$OUT/matrix_keep_aaf.tsv" "$OUT/matrix_keep_aaf_from_sketch.tsv"
run "$BIN" matrix --format full --metric ctx-naive --progress on -o "$OUT/matrix_progress.tsv" "$OUT/matrix_keep_sketch" 2> "$OUT/matrix_progress.err"
require_file_contains "$OUT/matrix_progress.err" "matrix: indexed full"
require_file_not_contains "$OUT/matrix_progress.tsv" "matrix:"
run "$BIN" matrix --format full --metric ctx-naive --progress on "$OUT/matrix_keep_sketch" > "$OUT/matrix_progress_stdout.tsv" 2> "$OUT/matrix_progress_stdout.err"
require_nonempty "$OUT/matrix_progress_stdout.tsv"
require_file_not_contains "$OUT/matrix_progress_stdout.tsv" "matrix:"
require_file_not_contains "$OUT/matrix_progress_stdout.err" "matrix:"
run "$BIN" ani -r "$OUT/ref_sketch" -q "$OUT/qry_sketch" -m1 -s -4 -o "$OUT/ani_m1_rect_direct.tsv"
run env KSSD3A_ANI_MATRIX_DIRECT_THRESHOLD=0 "$BIN" ani -r "$OUT/ref_sketch" -q "$OUT/qry_sketch" -m1 -s -4 -o "$OUT/ani_m1_rect_indexed.tsv"
require_same_file "$OUT/ani_m1_rect_direct.tsv" "$OUT/ani_m1_rect_indexed.tsv"
run "$BIN" ani -q "$OUT/matrix_keep_sketch" -m1 -s -4 -o "$OUT/ani_m1_self_direct.tsv"
run env KSSD3A_ANI_MATRIX_DIRECT_THRESHOLD=0 "$BIN" ani -q "$OUT/matrix_keep_sketch" -m1 -s -4 -o "$OUT/ani_m1_self_indexed.tsv"
require_same_file "$OUT/ani_m1_self_direct.tsv" "$OUT/ani_m1_self_indexed.tsv"
run "$BIN" ani -q "$OUT/matrix_keep_sketch" -m2 -s -4 -o "$OUT/ani_m2_self_direct.tsv"
run env KSSD3A_ANI_MATRIX_DIRECT_THRESHOLD=0 "$BIN" ani -q "$OUT/matrix_keep_sketch" -m2 -s -4 -o "$OUT/ani_m2_self_indexed.tsv"
require_same_file "$OUT/ani_m2_self_direct.tsv" "$OUT/ani_m2_self_indexed.tsv"
run "$BIN" sketch --conflict -f0 -p2 -o "$OUT/conflict_pair_sketch" "$REFS/refA.fna" "$REFS/refA_frag.fna"
run "$BIN" ani -q "$OUT/conflict_pair_sketch" -m2 -s -4 -o "$OUT/ani_m2_conflict_direct.tsv"
run env KSSD3A_ANI_MATRIX_DIRECT_THRESHOLD=0 "$BIN" ani -q "$OUT/conflict_pair_sketch" -m2 -s -4 -o "$OUT/ani_m2_conflict_indexed.tsv"
require_same_file "$OUT/ani_m2_conflict_direct.tsv" "$OUT/ani_m2_conflict_indexed.tsv"
run "$BIN" matrix --format triangle --metric ctx-naive --progress on -o "$OUT/matrix_triangle_progress.tsv" "$OUT/matrix_keep_sketch" 2> "$OUT/matrix_triangle_progress.err"
require_file_contains "$OUT/matrix_triangle_progress.err" "matrix: indexed triangle lower pairs"
awk -F '\t' '
  NR == FNR {
    if (FNR == 1) next
    row = FNR - 1
    for (i = 2; i <= NF; ++i) full[row, i - 1] = $i
    next
  }
  {
    row = FNR
    for (i = 2; i <= NF; ++i)
      if ($i != full[row, i - 1]) exit 1
  }
' "$OUT/matrix_keep_dist_from_sketch.tsv" "$OUT/matrix_triangle_progress.tsv"
run "$BIN" sketch --psmp "$OUT/matrix_keep_sketch" > "$OUT/matrix_keep_sketch.names"
cut -f2- "$OUT/matrix_keep_sketch.names" > "$OUT/matrix_keep_sketch.sample_names"
require_line_count "$OUT/matrix_keep_sketch.sample_names" 2
require_file_has "$OUT/matrix_keep_sketch.sample_names" "$REFS/refA.fna"
require_file_has "$OUT/matrix_keep_sketch.sample_names" "$REFS/refB.fna"
require_file_lacks "$OUT/matrix_keep_sketch.sample_names" "$REFS/refA_frag.fna"
run "$BIN" sketch --dedup 0.001 --metric ctx-naive -o "$OUT/dedup_copy_sketch" "$OUT/dedup_src"
run "$BIN" sketch --psmp "$OUT/dedup_copy_sketch" > "$OUT/dedup_copy_sketch.names"
cut -f2- "$OUT/dedup_copy_sketch.names" > "$OUT/dedup_copy_sketch.sample_names"
require_line_count "$OUT/dedup_copy_sketch.sample_names" 2
require_file_has "$OUT/dedup_copy_sketch.sample_names" "$REFS/refA.fna"
require_file_has "$OUT/dedup_copy_sketch.sample_names" "$REFS/refB.fna"
require_file_lacks "$OUT/dedup_copy_sketch.sample_names" "$REFS/refA_frag.fna"
run "$BIN" sketch --dedup 0.001 --metric ctx-naive --dedup-strategy full-linkage -o "$OUT/dedup_full_sketch" "$OUT/dedup_src"
run "$BIN" sketch --psmp "$OUT/dedup_full_sketch" > "$OUT/dedup_full_sketch.names"
cut -f2- "$OUT/dedup_full_sketch.names" > "$OUT/dedup_full_sketch.sample_names"
require_line_count "$OUT/dedup_full_sketch.sample_names" 2
require_file_has "$OUT/dedup_full_sketch.sample_names" "$REFS/refA.fna"
require_file_has "$OUT/dedup_full_sketch.sample_names" "$REFS/refB.fna"
require_file_lacks "$OUT/dedup_full_sketch.sample_names" "$REFS/refA_frag.fna"
run "$BIN" sketch --dedup 0.001 --metric 'ctx-naive&aaf' -o "$OUT/dedup_and_sketch" "$OUT/dedup_src"
run "$BIN" sketch --psmp "$OUT/dedup_and_sketch" > "$OUT/dedup_and_sketch.names"
cut -f2- "$OUT/dedup_and_sketch.names" > "$OUT/dedup_and_sketch.sample_names"
require_line_count "$OUT/dedup_and_sketch.sample_names" 2
require_file_has "$OUT/dedup_and_sketch.sample_names" "$REFS/refA.fna"
require_file_has "$OUT/dedup_and_sketch.sample_names" "$REFS/refB.fna"
require_file_lacks "$OUT/dedup_and_sketch.sample_names" "$REFS/refA_frag.fna"
run "$BIN" sketch --dedup 1 --metric p_dist -o "$OUT/dedup_pdist_sketch" "$OUT/dedup_src"
run "$BIN" sketch --psmp "$OUT/dedup_pdist_sketch" > "$OUT/dedup_pdist_sketch.names"
cut -f2- "$OUT/dedup_pdist_sketch.names" > "$OUT/dedup_pdist_sketch.sample_names"
require_file_has "$OUT/dedup_pdist_sketch.sample_names" "$REFS/refA.fna"
run "$BIN" sketch --dedup 0.001 --metric ctx-naive --dedup-max-afcut 0.8 --dedup-ctxcut 1 -o "$OUT/dedup_guarded_sketch" "$OUT/dedup_src"
run "$BIN" sketch --psmp "$OUT/dedup_guarded_sketch" > "$OUT/dedup_guarded_sketch.names"
cut -f2- "$OUT/dedup_guarded_sketch.names" > "$OUT/dedup_guarded_sketch.sample_names"
require_line_count "$OUT/dedup_guarded_sketch.sample_names" 2
require_file_has "$OUT/dedup_guarded_sketch.sample_names" "$REFS/refA.fna"
require_file_has "$OUT/dedup_guarded_sketch.sample_names" "$REFS/refB.fna"
if run "$BIN" sketch --dedup-max-afcut 0.5 "$OUT/dedup_src" 2> "$OUT/dedup_guard_without_dedup.err"; then
  printf 'smoke: --dedup-max-afcut should require --dedup\n' >&2
  exit 1
fi
require_file_contains "$OUT/dedup_guard_without_dedup.err" "only used with --dedup"
if run "$BIN" sketch --dedup-strategy full-linkage "$OUT/dedup_src" 2> "$OUT/dedup_strategy_without_dedup.err"; then
  printf 'smoke: --dedup-strategy should require --dedup\n' >&2
  exit 1
fi
require_file_contains "$OUT/dedup_strategy_without_dedup.err" "only used with --dedup"
cp -a "$OUT/dedup_src" "$OUT/dedup_inplace"
run "$BIN" sketch -i "$OUT/dedup_inplace"
run "$BIN" sketch --dedup 0.001 --metric ctx-moe "$OUT/dedup_inplace"
run "$BIN" sketch --psmp "$OUT/dedup_inplace" > "$OUT/dedup_inplace.names"
cut -f2- "$OUT/dedup_inplace.names" > "$OUT/dedup_inplace.sample_names"
require_line_count "$OUT/dedup_inplace.sample_names" 2
require_file_has "$OUT/dedup_inplace.sample_names" "$REFS/refA.fna"
require_file_has "$OUT/dedup_inplace.sample_names" "$REFS/refB.fna"
require_file_lacks "$OUT/dedup_inplace.sample_names" "$REFS/refA_frag.fna"
if [ -e "$OUT/dedup_inplace/sortedcomb_ctxgid64obj32" ]; then
  printf 'smoke: --dedup should delete stale sortedcomb_ctxgid64obj32\n' >&2
  exit 1
fi
run "$BIN" sketch --dedup 0.001 --metric ctx-naive -f0 -p2 -o "$OUT/dedup_raw_build" "$REFS/refA_frag.fna" "$REFS/refA.fna" "$REFS/refB.fna"
run "$BIN" sketch --psmp "$OUT/dedup_raw_build" > "$OUT/dedup_raw_build.names"
cut -f2- "$OUT/dedup_raw_build.names" > "$OUT/dedup_raw_build.sample_names"
require_line_count "$OUT/dedup_raw_build.sample_names" 2
require_file_has "$OUT/dedup_raw_build.sample_names" "$REFS/refA.fna"
require_file_has "$OUT/dedup_raw_build.sample_names" "$REFS/refB.fna"
require_file_lacks "$OUT/dedup_raw_build.sample_names" "$REFS/refA_frag.fna"
mkdir -p "$OUT/dedup_raw_empty"
run "$BIN" sketch --dedup 0.001 --metric ctx-naive -f0 -p2 -o "$OUT/dedup_raw_empty" "$REFS/refA_frag.fna" "$REFS/refA.fna"
run "$BIN" sketch --psmp "$OUT/dedup_raw_empty" > "$OUT/dedup_raw_empty.names"
cut -f2- "$OUT/dedup_raw_empty.names" > "$OUT/dedup_raw_empty.sample_names"
require_line_count "$OUT/dedup_raw_empty.sample_names" 1
require_file_has "$OUT/dedup_raw_empty.sample_names" "$REFS/refA.fna"
require_file_lacks "$OUT/dedup_raw_empty.sample_names" "$REFS/refA_frag.fna"
run "$BIN" sketch --dedup 0.001 --metric ctx-naive -f0 -p2 -o "$OUT/dedup_raw_update" "$REFS/refA_frag.fna" "$REFS/refB.fna"
if run "$BIN" sketch --dedup 0.001 --metric ctx-naive -f0 -p2 -o "$OUT/dedup_raw_update" "$REFS/refA.fna" 2> "$OUT/dedup_raw_update.err"; then
  printf 'smoke: raw --dedup -o should refuse an existing sketch output\n' >&2
  exit 1
fi
require_file_contains "$OUT/dedup_raw_update.err" "cannot exactly update an existing sketch"
run "$BIN" sketch --psmp "$OUT/dedup_raw_update" > "$OUT/dedup_raw_update.names"
cut -f2- "$OUT/dedup_raw_update.names" > "$OUT/dedup_raw_update.sample_names"
require_line_count "$OUT/dedup_raw_update.sample_names" 2
require_file_has "$OUT/dedup_raw_update.sample_names" "$REFS/refA_frag.fna"
require_file_has "$OUT/dedup_raw_update.sample_names" "$REFS/refB.fna"
require_file_lacks "$OUT/dedup_raw_update.sample_names" "$REFS/refA.fna"
cp -a "$OUT/ref_sketch" "$OUT/append_sketch"
run "$BIN" sketch --append -o "$OUT/append_copy_sketch" "$OUT/ref_sketch" "$OUT/qry_sketch"
run "$BIN" set -P "$OUT/append_copy_sketch" > "$OUT/append_copy_sketch.names"
cut -f2- "$OUT/append_copy_sketch.names" > "$OUT/append_copy_sketch.sample_names"
require_line_count "$OUT/append_copy_sketch.sample_names" 4
require_file_has "$OUT/append_copy_sketch.sample_names" "$REFS/refA.fna"
require_file_has "$OUT/append_copy_sketch.sample_names" "$QRYS/qryA.fna"
run "$BIN" sketch --append "$OUT/append_sketch" "$OUT/qry_sketch"
run "$BIN" set -P "$OUT/append_sketch" > "$OUT/append_sketch.names"
cut -f2- "$OUT/append_sketch.names" > "$OUT/append_sketch.sample_names"
require_line_count "$OUT/append_sketch.sample_names" 4
require_file_has "$OUT/append_sketch.sample_names" "$REFS/refA.fna"
require_file_has "$OUT/append_sketch.sample_names" "$QRYS/qryA.fna"
if [ -e "$OUT/append_sketch/sortedcomb_ctxgid64obj32" ]; then
  printf 'smoke: --append should delete stale sortedcomb_ctxgid64obj32\n' >&2
  exit 1
fi
require_nonempty "$OUT/pos_sketch_many/comblco.position"
require_nonempty "$OUT/pos_sketch_one/comblco.position"
require_same_size "$OUT/pos_sketch_many/comblco" "$OUT/pos_sketch_many/comblco.position"
require_same_size "$OUT/pos_sketch_one/comblco" "$OUT/pos_sketch_one/comblco.position"
run "$BIN" set --ppos "$OUT/pos_sketch_many" > "$OUT/pos_sketch_many.positions.tsv"
run "$BIN" sketch --ppos "$OUT/pos_sketch_many" > "$OUT/pos_sketch_many.sketch_positions.tsv"
require_same_file "$OUT/pos_sketch_many.positions.tsv" "$OUT/pos_sketch_many.sketch_positions.tsv"
require_nonempty "$OUT/pos_sketch_many.positions.tsv"
cp -a "$OUT/pos_sketch_many" "$OUT/pos_remove_sketch"
printf '%s\n' "$REFS/refA.fna" > "$OUT/pos_remove_names.txt"
run "$BIN" sketch --remove "$OUT/pos_remove_names.txt" "$OUT/pos_remove_sketch"
require_same_size "$OUT/pos_remove_sketch/comblco" "$OUT/pos_remove_sketch/comblco.position"
run "$BIN" set --ppos "$OUT/pos_remove_sketch" > "$OUT/pos_remove_sketch.positions.tsv"
run "$BIN" sketch --ppos "$OUT/pos_remove_sketch" > "$OUT/pos_remove_sketch.sketch_positions.tsv"
require_same_file "$OUT/pos_remove_sketch.positions.tsv" "$OUT/pos_remove_sketch.sketch_positions.tsv"
require_nonempty "$OUT/pos_remove_sketch.positions.tsv"
cp -a "$OUT/pos_sketch_many" "$OUT/pos_keep_sketch"
printf '%s\n' "$REFS/refB.fna" > "$OUT/pos_keep_names.txt"
run "$BIN" sketch --keep "$OUT/pos_keep_names.txt" "$OUT/pos_keep_sketch"
require_same_size "$OUT/pos_keep_sketch/comblco" "$OUT/pos_keep_sketch/comblco.position"
run "$BIN" set --ppos "$OUT/pos_keep_sketch" > "$OUT/pos_keep_sketch.positions.tsv"
run "$BIN" sketch --ppos "$OUT/pos_keep_sketch" > "$OUT/pos_keep_sketch.sketch_positions.tsv"
require_same_file "$OUT/pos_keep_sketch.positions.tsv" "$OUT/pos_keep_sketch.sketch_positions.tsv"
require_nonempty "$OUT/pos_keep_sketch.positions.tsv"
run "$BIN" sketch --keep "$OUT/pos_keep_names.txt" --drop-position -o "$OUT/pos_keep_drop_sketch" "$OUT/pos_sketch_many"
require_nonempty "$OUT/pos_keep_drop_sketch/comblco"
if [ -e "$OUT/pos_keep_drop_sketch/comblco.position" ]; then
  printf 'smoke: --drop-position should omit comblco.position\n' >&2
  exit 1
fi

printf 'smoke: sketch stdin\n' >&2
"$BIN" sketch --conflict -f0 -o "$OUT/stdin_sketch" - < "$READS/reads.fastq"
run "$BIN" sketch --pipecmd "cat {}" --conflict -f0 -o "$OUT/pipecmd_sketch" "$READS/reads.fastq"
require_nonempty "$OUT/stdin_sketch/lcofiles.stat"
require_nonempty "$OUT/pipecmd_sketch/lcofiles.stat"

run "$BIN" ani -r "$OUT/ref_sketch" -q "$OUT/qry_sketch" -f0 -n0 -m0 -o "$OUT/ani_detail.tsv"
run "$BIN" ani -r "$OUT/ref_sketch" -q "$OUT/read_sketch" -f0 -n0 -m0 -o "$OUT/ani_query_conflict_sketch.tsv"
run "$BIN" ani -r "$OUT/ref_sketch" -q "$OUT/read_sketch_noconflict" -f0 -n0 -m0 -o "$OUT/ani_query_noconflict_sketch.tsv"
run "$BIN" ani -r "$OUT/ref_sketch" --qraw "$OUT/read_sketch" -f0 -n0 -m0 -o "$OUT/ani_qraw.tsv"
run "$BIN" ani -r "$OUT/ref_sketch" --qraw "$READS/reads.fastq" -f0 -n0 -m0 -o "$OUT/ani_qraw_direct_fastq.tsv"
run "$BIN" ani -r "$OUT/ref_sketch" --qraw "$READS/reads.fastq.gz" -f0 -n0 -m0 -o "$OUT/ani_qraw_direct_fastq_gz.tsv"
run "$BIN" ani -r "$OUT/ref_sketch" --qraw "$READS/reads.fastq.gz" -n0 -m0 -o "$OUT/ani_qraw_default_afcut.tsv"
run "$BIN" ani -r "$OUT/ref_sketch" --qraw "$READS/reads.fastq" -s3 -f0 -n0 -m0 -o "$OUT/ani_qraw_s3_default.tsv"
awk -F '\t' 'NR > 1 && $6 == "Naive" { found = 1 } END { if (!found) exit 1 }' "$OUT/ani_qraw_s3_default.tsv"
run "$BIN" ani -r "$OUT/ref_sketch" --qraw "$READS/reads.fastq" -s9 -f0 -n0 -m0 -o "$OUT/ani_qraw_pdist.tsv"
awk -F '\t' 'NR > 1 && $6 == "p_dist" { found = 1 } END { if (!found) exit 1 }' "$OUT/ani_qraw_pdist.tsv"
run "$BIN" ani -r "$OUT/ref_sketch" --qraw "$REFS/refA.fna" -f0 -n0 -m0 -o "$OUT/ani_qraw_identical.tsv"
awk -F '\t' 'NR > 1 && $2 ~ /refA\.fna$/ && $4 <= 1e-12 && $6 == "Naive" { found = 1 } END { if (!found) exit 1 }' "$OUT/ani_qraw_identical.tsv"
run "$BIN" ani -r "$OUT/ref_sketch" --qraw "$READS/reads.fastq" -s3 --unified-metric -f0 -n0 -m0 -o "$OUT/ani_qraw_s3_unified.tsv"
awk -F '\t' 'NR > 1 && $6 == "CtxMoE" { found = 1 } END { if (!found) exit 1 }' "$OUT/ani_qraw_s3_unified.tsv"
run "$BIN" ani -r "$OUT/ref_sketch" --qraw "$OUT/read_sketch_noconflict" -f0 -n0 -m0 -o "$OUT/ani_qraw_noconflict_sketch.tsv"
run "$BIN" ani -r "$OUT/ref_sketch" -q "$READS/reads.fastq.gz" -f0 -n0 -m0 -o "$OUT/ani_query_direct_fastq_gz_noconflict.tsv"
run "$BIN" ani --conflict -r "$OUT/ref_sketch" -q "$READS/reads.fastq.gz" -f0 -n0 -m0 -o "$OUT/ani_query_direct_fastq_gz.tsv"
run "$BIN" sketch -f0 -p2 -o "$OUT/afcut_ref_sketch" "$AFREFS/refRefPlusExtra.fna"
run "$BIN" sketch -f0 -p2 -o "$OUT/afcut_qry_sketch" "$AFQRYS/qrySubset.fna"
run "$BIN" ani -r "$OUT/afcut_ref_sketch" -q "$OUT/afcut_qry_sketch" -f0.9 -n0 -m0 -o "$OUT/ani_max_afcut.tsv"
awk -F '\t' 'NR > 1 && $8 >= 0.9 && $10 < 0.9 { found = 1 } END { if (!found) exit 1 }' "$OUT/ani_max_afcut.tsv"
run "$BIN" ani --DimRdcFold 0 -f0 -n0 -o "$OUT/ani_auto_fasta.tsv" "$OUT/ref_sketch" "$QRYS/qryA.fna"
run "$BIN" ani --DimRdcFold 0 -f0 -n0 -o "$OUT/ani_auto_fastq.tsv" "$OUT/ref_sketch" "$READS/reads.fastq"
printf 'smoke: ani stdin raw query\n' >&2
"$BIN" ani --DimRdcFold 0 --conflict -f0 -n0 -o "$OUT/ani_auto_stdin.tsv" "$OUT/ref_sketch" - < "$READS/reads.fastq"
run "$BIN" ani --DimRdcFold 0 --conflict --pipecmd "cat {}" -f0 -n0 -o "$OUT/ani_auto_pipecmd.tsv" "$OUT/ref_sketch" "$READS/reads.fastq"
printf '%s\n' "$OUT/ref_sketch" > "$LISTS/refs.lst"
printf '%s\n' "$QRYS/qryA.fna" "$OUT/qry_sketch" > "$LISTS/qrys.lst"
run "$BIN" ani --DimRdcFold 0 --reflist "$LISTS/refs.lst" --qrylist "$LISTS/qrys.lst" -f0 -n0 -o "$OUT/ani_list.tsv"
require_nonempty "$OUT/ani_detail.tsv"
require_nonempty "$OUT/ani_query_conflict_sketch.tsv"
require_nonempty "$OUT/ani_query_noconflict_sketch.tsv"
require_nonempty "$OUT/ani_qraw.tsv"
require_nonempty "$OUT/ani_qraw_direct_fastq.tsv"
require_nonempty "$OUT/ani_qraw_direct_fastq_gz.tsv"
require_nonempty "$OUT/ani_qraw_default_afcut.tsv"
require_nonempty "$OUT/ani_qraw_s3_default.tsv"
require_nonempty "$OUT/ani_qraw_pdist.tsv"
require_nonempty "$OUT/ani_qraw_s3_unified.tsv"
require_nonempty "$OUT/ani_qraw_noconflict_sketch.tsv"
require_nonempty "$OUT/ani_query_direct_fastq_gz_noconflict.tsv"
require_nonempty "$OUT/ani_query_direct_fastq_gz.tsv"
require_same_tsv_except_first_col "$OUT/ani_qraw.tsv" "$OUT/ani_qraw_direct_fastq_gz.tsv"
require_same_tsv_except_first_col "$OUT/ani_query_conflict_sketch.tsv" "$OUT/ani_query_direct_fastq_gz.tsv"
require_same_tsv_except_first_col "$OUT/ani_query_noconflict_sketch.tsv" "$OUT/ani_query_direct_fastq_gz_noconflict.tsv"
require_nonempty "$OUT/ani_max_afcut.tsv"
require_nonempty "$OUT/ani_auto_fasta.tsv"
require_nonempty "$OUT/ani_auto_fastq.tsv"
require_nonempty "$OUT/ani_auto_stdin.tsv"
require_nonempty "$OUT/ani_auto_pipecmd.tsv"
require_nonempty "$OUT/ani_list.tsv"
require_header_has "$OUT/ani_auto_fasta.tsv" Selected_metric
require_header_has "$OUT/ani_auto_fastq.tsv" Selected_metric
require_header_has "$OUT/ani_auto_stdin.tsv" Selected_metric
require_header_has "$OUT/ani_auto_pipecmd.tsv" Selected_metric
for ani_tsv in \
  "$OUT/ani_detail.tsv" \
  "$OUT/ani_query_conflict_sketch.tsv" \
  "$OUT/ani_query_noconflict_sketch.tsv" \
  "$OUT/ani_qraw.tsv" \
  "$OUT/ani_qraw_direct_fastq.tsv" \
  "$OUT/ani_qraw_direct_fastq_gz.tsv" \
  "$OUT/ani_qraw_default_afcut.tsv" \
  "$OUT/ani_qraw_s3_default.tsv" \
  "$OUT/ani_qraw_s3_unified.tsv" \
  "$OUT/ani_qraw_noconflict_sketch.tsv" \
  "$OUT/ani_query_direct_fastq_gz_noconflict.tsv" \
  "$OUT/ani_query_direct_fastq_gz.tsv" \
  "$OUT/ani_max_afcut.tsv" \
  "$OUT/ani_auto_fasta.tsv" \
  "$OUT/ani_auto_fastq.tsv" \
  "$OUT/ani_auto_stdin.tsv" \
  "$OUT/ani_auto_pipecmd.tsv" \
  "$OUT/ani_list.tsv"; do
  require_blastn_af_le_one "$ani_tsv"
done

run "$BIN" matrix -r "$OUT/ref_sketch" -q "$OUT/qry_sketch" -o "$OUT/matrix_dist.tsv" -p2
run "$BIN" matrix -q "$OUT/qry_sketch" -o "$OUT/matrix_triangle.tsv" -p2
require_nonempty "$OUT/matrix_dist.tsv"
require_nonempty "$OUT/matrix_triangle.tsv"
awk -F '\t' -v refa="$REFS/refA.fna" -v refb="$REFS/refB.fna" \
  -v qrya="$QRYS/qryA.fna" -v qrya2="$QRYS/qryA2.fna" '
  NR == 1 && ($2 != refa || $3 != refb) {
    printf "smoke: matrix rectangular header should be reference columns, got %s and %s\n", $2, $3 > "/dev/stderr"
    exit 1
  }
  NR == 2 && $1 != qrya {
    printf "smoke: matrix rectangular row 1 should be query %s, got %s\n", qrya, $1 > "/dev/stderr"
    exit 1
  }
  NR == 3 && $1 != qrya2 {
    printf "smoke: matrix rectangular row 2 should be query %s, got %s\n", qrya2, $1 > "/dev/stderr"
    exit 1
  }
' "$OUT/matrix_dist.tsv"

run "$BIN" set --union -o "$OUT/union_sketch" "$OUT/ref_sketch"
run "$BIN" set --uniq_union -o "$OUT/uniq_union" "$OUT/ref_sketch"
run "$BIN" set --uniq_union --markerdb -o "$OUT/markerdb" "$OUT/ref_sketch"
require_nonempty "$OUT/union_sketch/lpan"
require_nonempty "$OUT/uniq_union/luniq_pan"
require_nonempty "$OUT/markerdb/comblco"
require_nonempty "$OUT/markerdb/comblco.index"
run "$BIN" set --intersect "$OUT/union_sketch" --key ctx -o "$OUT/intersect_ctx" "$OUT/qry_sketch"
run "$BIN" set --subtract "$OUT/union_sketch" --key ctx -o "$OUT/subtract_ctx" "$OUT/qry_sketch"
require_nonempty "$OUT/intersect_ctx/comblco"
require_nonempty "$OUT/intersect_ctx/comblco.index"
require_exists "$OUT/subtract_ctx/comblco"
require_exists "$OUT/subtract_ctx/comblco.index"
if "$BIN" set --union --key ctx -o "$OUT/rejected_ctx_union" "$OUT/ref_sketch" > "$OUT/rejected_ctx_union.log" 2>&1; then
  printf 'smoke: set --union --key ctx should be rejected for now\n' >&2
  exit 1
fi
require_file_contains "$OUT/rejected_ctx_union.log" "currently supports only --intersect/--intsect and --subtract"

run "$BIN" sketch -A -f0 -o "$OUT/qry_abundance" "$QRYS/qryA.fna" "$QRYS/qryMix.fna"
run "$BIN" set --intersect "$OUT/union_sketch" --key ctx -o "$OUT/qry_abundance_intersect_ctx" "$OUT/qry_abundance"
run "$BIN" set --subtract "$OUT/union_sketch" --key ctx -o "$OUT/qry_abundance_subtract_ctx" "$OUT/qry_abundance"
require_abundance_matches_comblco "$OUT/qry_abundance_intersect_ctx"
require_abundance_matches_comblco "$OUT/qry_abundance_subtract_ctx"
run "$BIN" set --psketch "$OUT/qry_abundance" > "$OUT/qry_abundance.psketch.tsv"
run "$BIN" set --psketch "$OUT/qry_abundance_intersect_ctx" > "$OUT/qry_abundance_intersect_ctx.psketch.tsv"
run "$BIN" set --psketch "$OUT/qry_abundance_subtract_ctx" > "$OUT/qry_abundance_subtract_ctx.psketch.tsv"
require_multiset_subset "$OUT/qry_abundance_intersect_ctx.psketch.tsv" "$OUT/qry_abundance.psketch.tsv"
require_multiset_subset "$OUT/qry_abundance_subtract_ctx.psketch.tsv" "$OUT/qry_abundance.psketch.tsv"
mkdir -p "$OUT/qry_abundance_noqc"
cp "$OUT/qry_abundance/lcofiles.stat" "$OUT/qry_abundance_noqc/lcofiles.stat"
cp "$OUT/qry_abundance/comblco" "$OUT/qry_abundance_noqc/comblco"
cp "$OUT/qry_abundance/comblco.index" "$OUT/qry_abundance_noqc/comblco.index"
cp "$OUT/qry_abundance/comblco.a" "$OUT/qry_abundance_noqc/comblco.a"
cp "$OUT/qry_abundance/lcofiles.infilemeta" "$OUT/qry_abundance_noqc/lcofiles.infilemeta"
run "$BIN" sketch --sketchQC -o "$OUT/qry_abundance_inferred_qc" "$OUT/qry_abundance_noqc"
run "$BIN" composite -r "$OUT/markerdb" -q "$OUT/qry_abundance" -o "$OUT/composite_profile.tsv"
require_nonempty "$OUT/qry_abundance/comblco.a"
require_nonempty "$OUT/qry_abundance_inferred_qc/lcofiles.qc"
require_nonempty "$OUT/composite_profile.tsv"

printf 'smoke: all main CLI checks passed\n' >&2
