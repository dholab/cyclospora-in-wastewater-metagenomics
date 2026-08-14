#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

if (( $# != 1 )); then
  printf 'Usage: %s <core-nt-database>\n' "$0" >&2
  exit 2
fi

database="$1"
baits=kmers/cyclospora_cayetanensis_rrna_baits.fasta
blast_tsv=work/core_nt_bait_exact_match_blast.tsv

if [[ ! -s "$baits" ]]; then
  printf 'Current bait FASTA is absent or empty: %s\n' "$baits" >&2
  exit 1
fi

mkdir -p work
blastn \
  -task blastn \
  -word_size 31 \
  -ungapped \
  -perc_identity 100 \
  -qcov_hsp_perc 100 \
  -dust no \
  -db "$database" \
  -query "$baits" \
  -outfmt '6 qseqid saccver staxids sscinames pident length mismatch gapopen qstart qend sstart send stitle' \
  -out "$blast_tsv"

bash scripts/build_core_nt_validated_index.sh "$blast_tsv"
