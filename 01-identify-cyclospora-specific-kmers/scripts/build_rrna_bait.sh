#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
project_root="$(pwd)"

threads="${RRNA_BAIT_THREADS:-8}"
output_root="${RRNA_BAIT_OUTPUT_ROOT:-$(pwd)}"
fixture_dir="${RRNA_BAIT_FIXTURE_DIR:-}"
offline="${RRNA_BAIT_OFFLINE:-0}"

mkdir -p "$output_root"/{curated,kmers,provenance,reports,work,logs}
mkdir -p "$output_root/background"

export PYTHONPATH="$(pwd)/src${PYTHONPATH:+:$PYTHONPATH}"

if [[ ! "$threads" =~ ^[1-9][0-9]*$ ]]; then
  echo "RRNA_BAIT_THREADS must be a positive integer" >&2
  exit 1
fi
if [[ "$offline" != "0" && "$offline" != "1" ]]; then
  echo "RRNA_BAIT_OFFLINE must be 0 or 1" >&2
  exit 1
fi

retry() {
  local attempt
  for attempt in 1 2 3 4 5; do
    "$@" && return 0
    sleep "$attempt"
  done
  return 1
}

require_nonempty() {
  if [[ ! -s "$1" ]]; then
    echo "Required input is absent or empty: $1" >&2
    return 1
  fi
}

provenance_file="$output_root/provenance/input_manifest.tsv"

cache_is_verified() {
  local name="$1"
  local release="$2"
  local source="$3"
  local local_path="$4"
  local expected_bytes="$5"
  local expected_sha256="$6"
  python -m rrna_bait.sources verify-provenance \
    "$provenance_file" \
    --name "$name" \
    --release "$release" \
    --url "$source" \
    --local-path "$local_path" \
    --project-root "$project_root" \
    --expected-bytes "$expected_bytes" \
    --expected-sha256 "$expected_sha256"
}

record_cache_provenance() {
  local name="$1"
  local release="$2"
  local source="$3"
  local local_path="$4"
  python -m rrna_bait.sources record-provenance \
    "$provenance_file" \
    --name "$name" \
    --release "$release" \
    --url "$source" \
    --local-path "$local_path" \
    --project-root "$project_root"
}

download_once() {
  local url="$1"
  local destination="$2"
  local expected_bytes="$3"
  local expected_sha256="$4"
  local temporary="${destination}.tmp.$$"
  rm -f -- "$temporary"
  mkdir -p "$(dirname "$destination")"
  if curl --fail --location --silent --show-error "$url" -o "$temporary" \
    && [[ -s "$temporary" ]] \
    && python -m rrna_bait.sources verify-identity "$temporary" \
      --expected-bytes "$expected_bytes" \
      --expected-sha256 "$expected_sha256"; then
    mv -f -- "$temporary" "$destination"
    return 0
  fi
  rm -f -- "$temporary"
  return 1
}

fetch_accessions_once() {
  local accession_file="$1"
  local destination="$2"
  local rules="$3"
  local expected_bytes="$4"
  local expected_sha256="$5"
  local temporary="${destination}.tmp.$$"
  local validated="${temporary}.validated"
  local accessions
  rm -f -- "$temporary"
  rm -f -- "$validated"
  accessions="$(paste -sd, "$accession_file")"
  if curl --fail --location --silent --show-error --http1.1 \
      --request POST \
      --data-urlencode "db=nuccore" \
      --data-urlencode "id=$accessions" \
      --data-urlencode "rettype=fasta" \
      --data-urlencode "retmode=text" \
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi" \
      -o "$temporary" \
    && [[ -s "$temporary" ]] \
    && python -m rrna_bait.targets prepare-queries \
      --downloaded-fasta "$temporary" \
      --rules "$rules" \
      --output-fasta "$validated" \
    && python -m rrna_bait.sources verify-identity "$temporary" \
      --expected-bytes "$expected_bytes" \
      --expected-sha256 "$expected_sha256"; then
    mv -f -- "$temporary" "$destination"
    rm -f -- "$validated"
    return 0
  fi
  rm -f -- "$temporary" "$validated"
  return 1
}

fetch_versioned_fasta_once() {
  local accession_file="$1"
  local destination="$2"
  local expected_bytes="$3"
  local expected_sha256="$4"
  local temporary="${destination}.tmp.$$"
  local normalized="${temporary}.normalized"
  local accessions
  rm -f -- "$temporary"
  rm -f -- "$normalized"
  accessions="$(paste -sd, "$accession_file")"
  if curl --fail --location --silent --show-error --http1.1 \
      --request POST \
      --data-urlencode "db=nuccore" \
      --data-urlencode "id=$accessions" \
      --data-urlencode "rettype=fasta" \
      --data-urlencode "retmode=text" \
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi" \
      -o "$temporary" \
    && [[ -s "$temporary" ]] \
    && awk 'NF' "$temporary" > "$normalized" \
    && python -m rrna_bait.sources verify-identity "$normalized" \
      --expected-bytes "$expected_bytes" \
      --expected-sha256 "$expected_sha256"; then
    mv -f -- "$normalized" "$destination"
    rm -f -- "$temporary"
    return 0
  fi
  rm -f -- "$temporary" "$normalized"
  return 1
}

prepare_fixture_loci() {
  local fasta="$1"
  local loci="$2"
  awk '
    BEGIN {
      OFS = "\t"
      print "target_record_id", "rrna_class", "query_accession", \
        "subject_accession", "start", "end", "strand"
    }
    /^>/ {
      if (record_id != "") {
        print record_id, rrna_class, record_id, record_id, 1, length(sequence), "+"
      }
      record_id = substr($0, 2)
      split(record_id, parts, "|")
      rrna_class = parts[1]
      sequence = ""
      next
    }
    {
      gsub(/[[:space:]]/, "")
      sequence = sequence toupper($0)
    }
    END {
      if (record_id != "") {
        print record_id, rrna_class, record_id, record_id, 1, length(sequence), "+"
      }
    }
  ' "$fasta" > "$loci"
}

combine_gzip_fastas() {
  local destination="$1"
  shift
  local temporary="${destination}.tmp.$$"
  gzip -cd -- "$@" | gzip -n -c > "$temporary"
  mv -f -- "$temporary" "$destination"
}

target_fasta="$output_root/curated/target_rrna.fasta"
target_loci="$output_root/curated/target_loci.tsv"
silva_background="$output_root/background/silva_non_target.fasta.gz"
rfam_background="$output_root/background/rfam_non_target.fasta.gz"
other_background="$output_root/background/other_cyclospora_rrna.fasta"

if [[ -n "$fixture_dir" ]]; then
  require_nonempty "$fixture_dir/build_curated_targets.fasta"
  require_nonempty "$fixture_dir/background_silva.fasta"
  require_nonempty "$fixture_dir/background_rfam.fasta"
  require_nonempty "$fixture_dir/background_related.fasta"

  cp "$fixture_dir/build_curated_targets.fasta" "$target_fasta"
  prepare_fixture_loci "$target_fasta" "$target_loci"
  cp "$fixture_dir/background_silva.fasta" "$output_root/background/silva_non_target.fasta"
  cp "$fixture_dir/background_rfam.fasta" "$output_root/background/rfam_non_target.fasta"
  cp "$fixture_dir/background_related.fasta" "$other_background"
  printf \
    'name\trelease\turl\tlocal_path\tretrieved_at_utc\tbytes\tsha256\n' \
    > "$provenance_file"
  silva_background="$output_root/background/silva_non_target.fasta"
  rfam_background="$output_root/background/rfam_non_target.fasta"
else
  declare -A source_release source_url source_local source_bytes source_sha256
  while IFS=$'\t' read -r name release url local_path expected_bytes expected_sha256; do
    [[ "$name" == "name" || -z "$name" ]] && continue
    source_release["$name"]="$release"
    source_url["$name"]="$url"
    source_local["$name"]="$local_path"
    source_bytes["$name"]="$expected_bytes"
    source_sha256["$name"]="$expected_sha256"
  done < config/sources.tsv

  config_accessions=()
  while IFS=$'\t' read -r accession _; do
    [[ "$accession" == "accession" || -z "$accession" ]] && continue
    config_accessions+=("$accession")
  done < config/target_queries.tsv

  silva_sources=()
  rfam_sources=()
  refseq_sources=()
  for name in \
    silva_ssu silva_lsu rfam_5s rfam_5_8s \
    refseq_ccayref3 refseq_asm76915v2; do
    release="${source_release[$name]}"
    url="${source_url[$name]}"
    local_path="${source_local[$name]}"
    expected_bytes="${source_bytes[$name]}"
    expected_sha256="${source_sha256[$name]}"
    source_path="$(pwd)/$local_path"
    if ! cache_is_verified \
      "$name" "$release" "$url" "$source_path" \
      "$expected_bytes" "$expected_sha256"; then
      if [[ "$offline" == "1" ]]; then
        echo "Offline mode requires a verified cached source: $source_path" >&2
        exit 1
      fi
      retry download_once \
        "$url" "$source_path" "$expected_bytes" "$expected_sha256"
      record_cache_provenance "$name" "$release" "$url" "$source_path"
    fi
    case "$name" in
      silva_*) silva_sources+=("$source_path") ;;
      rfam_*) rfam_sources+=("$source_path") ;;
      refseq_*) refseq_sources+=("$source_path") ;;
    esac
  done

  if (( ${#silva_sources[@]} == 0 || ${#rfam_sources[@]} == 0 || ${#refseq_sources[@]} != 2 )); then
    echo "config/sources.tsv must define SILVA, Rfam, and two RefSeq sources" >&2
    exit 1
  fi

  refseq_subject="$project_root/${source_local[refseq_nuclear_combined]}"
  python -m rrna_bait.sources reconstruct-gzip-concat \
    "$refseq_subject" "${refseq_sources[@]}" \
    --expected-bytes "${source_bytes[refseq_nuclear_combined]}" \
    --expected-sha256 "${source_sha256[refseq_nuclear_combined]}"

  query_accession_file="$output_root/work/target_query_accessions.txt"
  printf '%s\n' "${config_accessions[@]}" > "$query_accession_file"
  downloaded_queries="$output_root/work/target_queries_downloaded.fasta"
  target_query_source="${source_url[ncbi_eutils_target_queries]}"
  target_query_release="${source_release[ncbi_eutils_target_queries]}"
  target_query_bytes="${source_bytes[ncbi_eutils_target_queries]}"
  target_query_sha256="${source_sha256[ncbi_eutils_target_queries]}"
  query_cache_valid=0
  if cache_is_verified \
      "ncbi_eutils_target_queries" "$target_query_release" \
      "$target_query_source" "$downloaded_queries" \
      "$target_query_bytes" "$target_query_sha256" \
    && python -m rrna_bait.targets prepare-queries \
      --downloaded-fasta "$downloaded_queries" \
      --rules config/target_queries.tsv \
      --output-fasta "$output_root/work/target_queries.fasta"; then
    query_cache_valid=1
  fi
  if [[ "$query_cache_valid" == "0" ]]; then
    if [[ "$offline" == "1" ]]; then
      echo "Offline mode requires complete, verified cached target queries: $downloaded_queries" >&2
      exit 1
    fi
    retry fetch_accessions_once \
      "$query_accession_file" "$downloaded_queries" config/target_queries.tsv \
      "$target_query_bytes" "$target_query_sha256"
    record_cache_provenance \
      "ncbi_eutils_target_queries" "$target_query_release" \
      "$target_query_source" "$downloaded_queries"
  fi
  python -m rrna_bait.targets prepare-queries \
    --downloaded-fasta "$downloaded_queries" \
    --rules config/target_queries.tsv \
    --output-fasta "$output_root/work/target_queries.fasta"

  blastn -task megablast \
    -query "$output_root/work/target_queries.fasta" \
    -subject "$refseq_subject" \
    -outfmt '6 qseqid qlen qstart qend sseqid sstart send length pident' \
    -out "$output_root/work/target_hits.tsv"

  python -m rrna_bait.targets \
    --queries "$output_root/work/target_queries.fasta" \
    --subject "$refseq_subject" \
    --rules config/target_queries.tsv \
    --hits "$output_root/work/target_hits.tsv" \
    --fasta-out "$target_fasta" \
    --loci-out "$target_loci"

  target_rrna_accessions="$project_root/${source_local[target_rrna_accession_snapshot]}"
  python -m rrna_bait.sources verify-identity "$target_rrna_accessions" \
    --expected-bytes "${source_bytes[target_rrna_accession_snapshot]}" \
    --expected-sha256 "${source_sha256[target_rrna_accession_snapshot]}"

  silva_exclusions=()
  for accession in "${config_accessions[@]}"; do
    silva_exclusions+=(--excluded-accession "$accession")
  done
  silva_parts=()
  for index in "${!silva_sources[@]}"; do
    part="$output_root/work/silva_non_target_${index}.fasta.gz"
    python -m rrna_bait.sources filter-silva \
      "${silva_sources[$index]}" "$part" "${silva_exclusions[@]}"
    silva_parts+=("$part")
  done
  combine_gzip_fastas "$silva_background" "${silva_parts[@]}"

  rfam_exclusions=()
  for accession in "${config_accessions[@]}"; do
    rfam_exclusions+=(--excluded-accession "$accession")
  done
  while IFS= read -r accession; do
    [[ -z "$accession" ]] && continue
    rfam_exclusions+=(--excluded-accession "$accession")
  done < "$target_rrna_accessions"
  rfam_parts=()
  for index in "${!rfam_sources[@]}"; do
    part="$output_root/work/rfam_non_target_${index}.fasta.gz"
    python -m rrna_bait.sources filter-rfam \
      "${rfam_sources[$index]}" "$part" "${rfam_exclusions[@]}"
    rfam_parts+=("$part")
  done
  combine_gzip_fastas "$rfam_background" "${rfam_parts[@]}"

  other_accessions="$project_root/${source_local[other_cyclospora_accession_snapshot]}"
  python -m rrna_bait.sources verify-identity "$other_accessions" \
    --expected-bytes "${source_bytes[other_cyclospora_accession_snapshot]}" \
    --expected-sha256 "${source_sha256[other_cyclospora_accession_snapshot]}"
  other_rrna_source="${source_url[ncbi_eutils_other_cyclospora_rrna]}"
  other_rrna_release="${source_release[ncbi_eutils_other_cyclospora_rrna]}"
  other_rrna_bytes="${source_bytes[ncbi_eutils_other_cyclospora_rrna]}"
  other_rrna_sha256="${source_sha256[ncbi_eutils_other_cyclospora_rrna]}"
  if ! cache_is_verified \
    "ncbi_eutils_other_cyclospora_rrna" "$other_rrna_release" \
    "$other_rrna_source" "$other_background" \
    "$other_rrna_bytes" "$other_rrna_sha256"; then
    if [[ "$offline" == "1" ]]; then
      echo "Offline mode requires verified cached other-Cyclospora rRNAs: $other_background" >&2
      exit 1
    fi
    retry fetch_versioned_fasta_once \
      "$other_accessions" "$other_background" \
      "$other_rrna_bytes" "$other_rrna_sha256"
    record_cache_provenance \
      "ncbi_eutils_other_cyclospora_rrna" "$other_rrna_release" \
      "$other_rrna_source" "$other_background"
  fi
fi

require_nonempty "$target_fasta"
require_nonempty "$target_loci"
require_nonempty "$silva_background"
require_nonempty "$rfam_background"
require_nonempty "$other_background"

meryl_log="$output_root/logs/meryl.log"
: > "$meryl_log"
for database in \
  target silva rfam other_cyclospora static_non_target genus_compatible \
  all_non_target exact_specific \
  target_silva target_rfam target_other_cyclospora; do
  rm -rf -- "$output_root/work/${database}.meryl"
done

meryl count k=31 threads="$threads" "$target_fasta" \
  output "$output_root/work/target.meryl" >> "$meryl_log" 2>&1
meryl count k=31 threads="$threads" "$silva_background" \
  output "$output_root/work/silva.meryl" >> "$meryl_log" 2>&1
meryl count k=31 threads="$threads" "$rfam_background" \
  output "$output_root/work/rfam.meryl" >> "$meryl_log" 2>&1
meryl count k=31 threads="$threads" "$other_background" \
  output "$output_root/work/other_cyclospora.meryl" >> "$meryl_log" 2>&1
meryl union \
  "$output_root/work/silva.meryl" \
  "$output_root/work/rfam.meryl" \
  output "$output_root/work/static_non_target.meryl" >> "$meryl_log" 2>&1
meryl difference \
  "$output_root/work/target.meryl" \
  "$output_root/work/static_non_target.meryl" \
  output "$output_root/work/genus_compatible.meryl" >> "$meryl_log" 2>&1
meryl union \
  "$output_root/work/silva.meryl" \
  "$output_root/work/rfam.meryl" \
  "$output_root/work/other_cyclospora.meryl" \
  output "$output_root/work/all_non_target.meryl" >> "$meryl_log" 2>&1
meryl difference \
  "$output_root/work/target.meryl" \
  "$output_root/work/all_non_target.meryl" \
  output "$output_root/work/exact_specific.meryl" >> "$meryl_log" 2>&1

meryl intersect "$output_root/work/target.meryl" "$output_root/work/silva.meryl" \
  output "$output_root/work/target_silva.meryl" >> "$meryl_log" 2>&1
meryl intersect "$output_root/work/target.meryl" "$output_root/work/rfam.meryl" \
  output "$output_root/work/target_rfam.meryl" >> "$meryl_log" 2>&1
meryl intersect \
  "$output_root/work/target.meryl" \
  "$output_root/work/other_cyclospora.meryl" \
  output "$output_root/work/target_other_cyclospora.meryl" >> "$meryl_log" 2>&1

exact_print="$output_root/work/exact_specific.kmers.tsv"
genus_print="$output_root/work/genus_compatible_pre_entropy.kmers.tsv"
silva_print="$output_root/work/target_silva.kmers.tsv"
rfam_print="$output_root/work/target_rfam.kmers.tsv"
other_print="$output_root/work/target_other_cyclospora.kmers.tsv"
meryl print "$output_root/work/exact_specific.meryl" \
  > "$exact_print" 2>> "$meryl_log"
meryl print "$output_root/work/genus_compatible.meryl" \
  > "$genus_print" 2>> "$meryl_log"
meryl print "$output_root/work/target_silva.meryl" \
  > "$silva_print" 2>> "$meryl_log"
meryl print "$output_root/work/target_rfam.meryl" \
  > "$rfam_print" 2>> "$meryl_log"
meryl print "$output_root/work/target_other_cyclospora.meryl" \
  > "$other_print" 2>> "$meryl_log"

raw_manifest="$output_root/kmers/cyclospora_cayetanensis_rrna_specific.raw.kmers.tsv"
raw_baits="$output_root/work/raw_specific_baits.fasta"
final_manifest="$output_root/kmers/cyclospora_cayetanensis_rrna_specific.kmers.tsv"
final_baits="$output_root/kmers/cyclospora_cayetanensis_rrna_baits.fasta"
python -m rrna_bait.manifest raw \
  --target-fasta "$target_fasta" \
  --loci-tsv "$target_loci" \
  --silva-print "$silva_print" \
  --rfam-print "$rfam_print" \
  --other-cyclospora-print "$other_print" \
  --exact-difference-print "$exact_print" \
  --raw-manifest "$raw_manifest" \
  --raw-baits "$raw_baits" \
  -k 31

entropy_index="$output_root/work/entropy_filtered.idx"
entropy_pass="$output_root/work/entropy_pass.fasta"
rm -f -- "$entropy_index" "$entropy_pass"
raw_count="$(awk '/^>/{count++} END{print count+0}' "$raw_baits")"
if (( raw_count > 0 )); then
  deacon index build \
    -k 31 -w 1 -e 0.6 -t "$threads" \
    "$raw_baits" \
    -o "$entropy_index"
  deacon index dump "$entropy_index" -o "$entropy_pass"
else
  : > "$entropy_pass"
fi

python -m rrna_bait.manifest finalize \
  --raw-manifest "$raw_manifest" \
  --entropy-pass-fasta "$entropy_pass" \
  --final-manifest "$final_manifest" \
  --final-baits "$final_baits"

final_count="$(awk '/^>/{count++} END{print count+0}' "$final_baits")"
genus_count="$(awk 'NF{count++} END{print count+0}' "$genus_print")"
if (( genus_count == 0 )); then
  genus_level_status="ZERO_GENUS_CANDIDATES"
elif (( final_count > 0 )); then
  genus_level_status="SPECIES_SPECIFIC_CANDIDATES"
else
  genus_level_status="GENUS_ONLY_CANDIDATES"
fi
if (( final_count > 0 )); then
  build_status="FEASIBLE"
  conclusion="available mature-rRNA references support a C. cayetanensis-specific k=31 assay"
else
  build_status="INFEASIBLE"
  conclusion="available mature-rRNA references do not support a C. cayetanensis-specific k=31 assay"
fi
summary_report="$output_root/reports/build_summary.tsv"
{
  printf 'metric\tvalue\n'
  printf 'build_status\t%s\n' "$build_status"
  printf 'conclusion\t%s\n' "$conclusion"
  printf 'raw_exact_candidate_count\t%s\n' "$raw_count"
  printf 'final_bait_count\t%s\n' "$final_count"
  printf 'genus_compatible_pre_entropy_count\t%s\n' "$genus_count"
  printf 'genus_level_status\t%s\n' "$genus_level_status"
} > "$summary_report"

final_index="$output_root/cyclospora_cayetanensis_rrna_k31w1.idx"
verification_summary="$output_root/reports/cyclospora_cayetanensis_rrna_index_summary.tsv"
verification_report="$output_root/reports/cyclospora_cayetanensis_rrna_index_report.md"
rm -f -- "$final_index"
if (( final_count == 0 )); then
  verify_args=(
    --target "$target_fasta"
    --manifest "$final_manifest"
    --baits "$final_baits"
    --exact-difference "$exact_print"
    --silva-intersection "$silva_print"
    --rfam-intersection "$rfam_print"
    --other-cyclospora-intersection "$other_print"
    --genus-compatible "$genus_print"
    --provenance "$provenance_file"
    --loci "$target_loci"
    --summary "$verification_summary"
    --report "$verification_report"
  )
  if [[ -n "$fixture_dir" ]]; then
    verify_args+=(--allow-empty-provenance)
  fi
  python -m rrna_bait.verify "${verify_args[@]}"
  echo "No rRNA bait k-mers survived filtering; see $summary_report" >&2
  exit 2
fi

deacon index build \
  -k 31 -w 1 -e 0 -t "$threads" \
  "$final_baits" \
  -o "$final_index"

final_index_info="$output_root/work/final_index_info.txt"
deacon index info "$final_index" > "$final_index_info" 2>&1

bait_roundtrip="$output_root/work/bait_roundtrip.fasta"
rm -f -- "$bait_roundtrip"
deacon filter -a 1 -r 0 \
  "$final_index" \
  "$final_baits" \
  -o "$bait_roundtrip"

cmp \
  <(seqkit seq -s "$final_baits" | sort) \
  <(seqkit seq -s "$bait_roundtrip" | sort)

background_names=(silva rfam other_cyclospora)
background_fastas=("$silva_background" "$rfam_background" "$other_background")
for index in "${!background_names[@]}"; do
  background_name="${background_names[$index]}"
  background_roundtrip="$output_root/work/${background_name}_background_roundtrip.fasta"
  rm -f -- "$background_roundtrip"
  deacon filter -a 1 -r 0 \
    "$final_index" \
    "${background_fastas[$index]}" \
    -o "$background_roundtrip"
  retained_count="$(
    awk '/^>/{count++} END{print count+0}' "$background_roundtrip"
  )"
  if (( retained_count != 0 )); then
    echo \
      "Static background round trip retained $retained_count ${background_name} record(s)" \
      >&2
    exit 1
  fi
done

verify_args=(
  --target "$target_fasta"
  --manifest "$final_manifest"
  --baits "$final_baits"
  --index "$final_index"
  --exact-difference "$exact_print"
  --silva-intersection "$silva_print"
  --rfam-intersection "$rfam_print"
  --other-cyclospora-intersection "$other_print"
  --genus-compatible "$genus_print"
  --provenance "$provenance_file"
  --loci "$target_loci"
  --summary "$verification_summary"
  --report "$verification_report"
)
if [[ -n "$fixture_dir" ]]; then
  verify_args+=(--allow-empty-provenance)
fi
python -m rrna_bait.verify "${verify_args[@]}"
