#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

blast_tsv="${1:-results/core_nt_bait_exact_match_blast.tsv}"

if [[ ! -f "$blast_tsv" ]]; then
  printf 'Completed exact-match BLAST TSV is absent: %s\n' "$blast_tsv" >&2
  exit 1
fi

bash scripts/finalize_core_nt_validation.sh \
  --blast-tsv "$blast_tsv" \
  --input-baits kmers/cyclospora_cayetanensis_rrna_baits.fasta \
  --input-manifest kmers/cyclospora_cayetanensis_rrna_specific.kmers.tsv \
  --output-baits kmers/cyclospora_cayetanensis_rrna_core_nt_validated_baits.fasta \
  --output-manifest kmers/cyclospora_cayetanensis_rrna_core_nt_validated.kmers.tsv \
  --decision-tsv reports/cyclospora_cayetanensis_core_nt_validation.tsv \
  --near-hits-tsv reports/cyclospora_cayetanensis_core_nt_near_hits.tsv \
  --index cyclospora_cayetanensis_rrna_core_nt_validated_k31w1.idx
