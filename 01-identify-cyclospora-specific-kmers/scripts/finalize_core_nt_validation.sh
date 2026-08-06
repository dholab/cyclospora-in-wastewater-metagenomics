#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

usage() {
  printf 'Usage: %s --blast-tsv TSV [--input-baits FASTA] [--input-manifest TSV] [--output-baits FASTA] [--output-manifest TSV] [--decision-tsv TSV] [--near-hits-tsv TSV] [--index IDX] [--background FASTA ...]\n' "$0" >&2
}

blast_tsv=
input_baits=kmers/cyclospora_cayetanensis_rrna_baits.fasta
input_manifest=kmers/cyclospora_cayetanensis_rrna_specific.kmers.tsv
output_baits=kmers/cyclospora_cayetanensis_rrna_baits.fasta
output_manifest=kmers/cyclospora_cayetanensis_rrna_specific.kmers.tsv
decision_tsv=reports/cyclospora_cayetanensis_core_nt_validation.tsv
near_hits_tsv=reports/cyclospora_cayetanensis_core_nt_near_hits.tsv
index=cyclospora_cayetanensis_rrna_k31w1.idx
threads="${RRNA_BAIT_THREADS:-8}"
backgrounds=(
  background/silva_non_target.fasta.gz
  background/rfam_non_target.fasta.gz
  background/other_cyclospora_rrna.fasta
)
explicit_backgrounds=0

while (( $# )); do
  case "$1" in
    --blast-tsv) blast_tsv="${2:-}"; shift 2 ;;
    --input-baits) input_baits="${2:-}"; shift 2 ;;
    --input-manifest) input_manifest="${2:-}"; shift 2 ;;
    --output-baits) output_baits="${2:-}"; shift 2 ;;
    --output-manifest) output_manifest="${2:-}"; shift 2 ;;
    --decision-tsv) decision_tsv="${2:-}"; shift 2 ;;
    --near-hits-tsv) near_hits_tsv="${2:-}"; shift 2 ;;
    --index) index="${2:-}"; shift 2 ;;
    --background)
      if (( explicit_backgrounds == 0 )); then
        backgrounds=()
        explicit_backgrounds=1
      fi
      backgrounds+=("${2:-}")
      shift 2
      ;;
    *) usage; exit 2 ;;
  esac
done

if [[ -z "$blast_tsv" || ! -f "$blast_tsv" ]]; then
  printf 'Completed all-bait BLAST TSV is absent: %s\n' "$blast_tsv" >&2
  exit 1
fi
for required in "$input_baits" "$input_manifest"; do
  if [[ ! -s "$required" ]]; then
    printf 'Required finalization input is absent or empty: %s\n' "$required" >&2
    exit 1
  fi
done
for background in "${backgrounds[@]}"; do
  if [[ ! -s "$background" ]]; then
    printf 'Static background is absent or empty: %s\n' "$background" >&2
    exit 1
  fi
done
if [[ ! "$threads" =~ ^[1-9][0-9]*$ ]]; then
  printf 'RRNA_BAIT_THREADS must be a positive integer\n' >&2
  exit 1
fi

for output in "$output_baits" "$output_manifest" "$decision_tsv" "$near_hits_tsv" "$index"; do
  mkdir -p "$(dirname "$output")"
done

temporary="$(mktemp -d)"
trap 'rm -rf -- "$temporary"' EXIT
candidate_baits="$temporary/candidate_baits.fasta"
candidate_manifest="$temporary/candidate_manifest.tsv"
candidate_decisions="$temporary/candidate_decisions.tsv"
candidate_near_hits="$temporary/candidate_near_hits.tsv"
candidate_index="$temporary/candidate.idx"

export PYTHONPATH="$(pwd)/src${PYTHONPATH:+:$PYTHONPATH}"
python -m rrna_bait.core_nt \
  --baits "$input_baits" \
  --blast-tsv "$blast_tsv" \
  --manifest-in "$input_manifest" \
  --manifest-out "$candidate_manifest" \
  --baits-out "$candidate_baits" \
  --decisions-out "$candidate_decisions" \
  --near-hits-out "$candidate_near_hits"

retained_count="$(grep -c '^>' "$candidate_baits" || true)"
if (( retained_count == 0 )); then
  printf 'Core-nt validation rejected every bait; refusing to emit an empty index.\n' >&2
  exit 1
fi

deacon index build \
  -k 31 -w 1 -e 0 -t "$threads" \
  "$candidate_baits" \
  -o "$candidate_index"

roundtrip="$temporary/bait_roundtrip.fasta"
deacon filter -a 1 -r 0 \
  "$candidate_index" \
  "$candidate_baits" \
  -o "$roundtrip"
cmp \
  <(seqkit seq -s "$candidate_baits" | sort) \
  <(seqkit seq -s "$roundtrip" | sort)

for background in "${backgrounds[@]}"; do
  background_roundtrip="$temporary/$(basename "$background").roundtrip.fasta"
  deacon filter -a 1 -r 0 \
    "$candidate_index" \
    "$background" \
    -o "$background_roundtrip"
  retained_background="$(
    awk '/^>/{count++} END{print count+0}' "$background_roundtrip"
  )"
  if (( retained_background != 0 )); then
    printf 'Static background retained %s record(s): %s\n' \
      "$retained_background" "$background" >&2
    exit 1
  fi
done

mv -f -- "$candidate_baits" "$output_baits"
mv -f -- "$candidate_manifest" "$output_manifest"
mv -f -- "$candidate_decisions" "$decision_tsv"
mv -f -- "$candidate_near_hits" "$near_hits_tsv"
mv -f -- "$candidate_index" "$index"
printf 'Published %s core-nt-validated bait(s).\n' "$retained_count"
