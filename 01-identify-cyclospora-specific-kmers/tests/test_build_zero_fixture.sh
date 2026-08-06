#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

fixture_root="$(mktemp -d)"
output_root="$(mktemp -d)"
trap 'rm -rf "$fixture_root" "$output_root"' EXIT

cp tests/fixtures/build_curated_targets.fasta "$fixture_root/build_curated_targets.fasta"
cp tests/fixtures/background_silva.fasta "$fixture_root/background_silva.fasta"
cp tests/fixtures/background_rfam.fasta "$fixture_root/background_rfam.fasta"
cp tests/fixtures/build_curated_targets.fasta "$fixture_root/background_related.fasta"

printf 'stale index\n' > "$output_root/cyclospora_cayetanensis_rrna_k31w1.idx"

set +e
RRNA_BAIT_FIXTURE_DIR="$fixture_root" \
RRNA_BAIT_OUTPUT_ROOT="$output_root" \
RRNA_BAIT_THREADS=2 \
bash scripts/build_rrna_bait.sh
status=$?
set -e

test "$status" -eq 2
test ! -e "$output_root/cyclospora_cayetanensis_rrna_k31w1.idx"

report="$output_root/reports/build_summary.tsv"
test -s "$report"
grep -q $'^build_status\tINFEASIBLE$' "$report"
grep -q $'^conclusion\tavailable mature-rRNA references do not support a C. cayetanensis-specific k=31 assay$' "$report"
grep -q $'^genus_level_status\tGENUS_ONLY_CANDIDATES$' "$report"
awk -F '\t' '$1 == "genus_compatible_pre_entropy_count" && $2 > 0 {found=1} END {exit !found}' "$report"

verification_summary="$output_root/reports/cyclospora_cayetanensis_rrna_index_summary.tsv"
verification_report="$output_root/reports/cyclospora_cayetanensis_rrna_index_report.md"
test -s "$verification_summary"
test -s "$verification_report"
grep -q $'^conclusion\tINFEASIBLE: No species-specific mature-rRNA k=31 bait survived; no index was emitted\\.$' "$verification_summary"
grep -q '^INFEASIBLE: No species-specific mature-rRNA k=31 bait survived; no index was emitted\.$' "$verification_report"
grep -q '^Wastewater performance has not been evaluated\.' "$verification_report"
