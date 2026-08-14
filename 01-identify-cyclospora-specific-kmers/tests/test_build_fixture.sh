#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

fixture_root="$(pwd)/tests/fixtures"
output_root="$(mktemp -d)"
trap 'rm -rf "$output_root"' EXIT

RRNA_BAIT_FIXTURE_DIR="$fixture_root" \
RRNA_BAIT_OUTPUT_ROOT="$output_root" \
RRNA_BAIT_THREADS=2 \
bash scripts/build_rrna_bait.sh

test -s "$output_root/kmers/cyclospora_cayetanensis_rrna_baits.fasta"
test -s "$output_root/cyclospora_cayetanensis_rrna_k31w1.idx"
test -s "$output_root/kmers/cyclospora_cayetanensis_rrna_specific.kmers.tsv"
test -s "$output_root/reports/cyclospora_cayetanensis_rrna_index_summary.tsv"
test -s "$output_root/reports/cyclospora_cayetanensis_rrna_index_report.md"
test -s "$output_root/work/final_index_info.txt"
test -s "$output_root/work/bait_roundtrip.fasta"

manifest="$output_root/kmers/cyclospora_cayetanensis_rrna_specific.kmers.tsv"
awk -F '\t' '
  NR == 1 {
    for (column = 1; column <= NF; column++) {
      if ($column == "rrna_classes") class_column = column
      if ($column == "status") status_column = column
    }
    next
  }
  $status_column == "PASS" {
    passing++
    if ($class_column != "5.8S") exit 1
  }
  END {
    if (passing != 34) exit 1
  }
' "$manifest"

bait_count="$(grep -c '^>' "$output_root/kmers/cyclospora_cayetanensis_rrna_baits.fasta")"
test "$bait_count" -eq 34

index_info="$(< "$output_root/work/final_index_info.txt")"
grep -q 'k.*31' <<< "$index_info"
grep -q 'w.*1' <<< "$index_info"
index_count="$(awk -F ': ' '/Distinct minimizer count/{print $2}' <<< "$index_info")"
test "$bait_count" -eq "$index_count"

cmp \
  <(seqkit seq -s "$output_root/kmers/cyclospora_cayetanensis_rrna_baits.fasta" | sort) \
  <(seqkit seq -s "$output_root/work/bait_roundtrip.fasta" | sort)

for background in silva rfam other_cyclospora; do
  roundtrip="$output_root/work/${background}_background_roundtrip.fasta"
  test -e "$roundtrip"
  test "$(awk '/^>/{count++} END{print count+0}' "$roundtrip")" -eq 0
done

summary="$output_root/reports/cyclospora_cayetanensis_rrna_index_summary.tsv"
report="$output_root/reports/cyclospora_cayetanensis_rrna_index_report.md"
grep -q $'^conclusion\tPASS: A species-specific mature-rRNA k=31 bait set was produced\\.$' "$summary"
grep -q '^PASS: A species-specific mature-rRNA k=31 bait set was produced\.$' "$report"
grep -q '^Wastewater performance has not been evaluated\.' "$report"
