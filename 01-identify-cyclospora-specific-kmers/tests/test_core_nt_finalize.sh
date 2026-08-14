#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

test_root="$(mktemp -d)"
trap 'rm -rf "$test_root"' EXIT

if bash scripts/build_core_nt_validated_index.sh >"$test_root/no_argument.log" 2>&1; then
  printf 'Finalization without explicit current-bait BLAST evidence succeeded unexpectedly.\n' >&2
  exit 1
else
  status=$?
fi
test "$status" -eq 2
grep -q '^Usage:' "$test_root/no_argument.log"

input_baits="$test_root/input.fasta"
input_manifest="$test_root/input.tsv"
blast_tsv="$test_root/blast.tsv"
output_baits="$test_root/output.fasta"
output_manifest="$test_root/output.tsv"
decisions="$test_root/decisions.tsv"
near_hits="$test_root/near_hits.tsv"
index="$test_root/output.idx"
background_silva="$test_root/silva.fasta"
background_rfam="$test_root/rfam.fasta"
background_other="$test_root/other.fasta"

cat > "$input_baits" <<'EOF'
>bait_reject|18S
AACCGGTTACGATCGTAGCTAGGCTAACGTA
>bait_target|28S
AGTCGATGCTACCGTTAAGGCTACGATTCGA
>bait_near|5S
ATGCCGTAGCATCGATGGCATACGTTAGCTA
EOF
printf '>silva\n%s\n' 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT' > "$background_silva"
printf '>rfam\n%s\n' 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC' > "$background_rfam"
printf '>other\n%s\n' 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG' > "$background_other"

cat > "$input_manifest" <<'EOF'
kmer	rrna_classes	target_records	target_starts	target_occurrences	target_copy_count	in_silva	in_rfam	in_other_cyclospora	in_exact_difference	status	rejection_reason
AACCGGTTACGATCGTAGCTAGGCTAACGTA	18S	18S|fixture	0	18S|fixture:0	1	0	0	0	1	PASS	none
AGTCGATGCTACCGTTAAGGCTACGATTCGA	28S	28S|fixture	0	28S|fixture:0	1	0	0	0	1	PASS	none
ATGCCGTAGCATCGATGGCATACGTTAGCTA	5S	5S|fixture	0	5S|fixture:0	1	0	0	0	1	PASS	none
EOF

cat > "$blast_tsv" <<'EOF'
bait_reject|18S	MK946294.1	342097	uncultured soil eukaryote	100.000	31	0	0	1	31	54	24	soil SSU
bait_target|28S	X.1	88456	Cyclospora cayetanensis	100.000	31	0	0	1	31	1	31	target
bait_near|5S	Y.1	1	other	96.774	31	1	0	1	31	1	31	near
EOF

fake_bin="$test_root/bin"
mkdir -p "$fake_bin"
cat > "$fake_bin/deacon" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
if [[ "$1 $2" == "index build" ]]; then
  output=
  input=
  for ((i=3; i <= $#; i++)); do
    value="${!i}"
    if [[ "$value" == "-o" ]]; then
      ((i++))
      output="${!i}"
    elif [[ "$value" != -* && "$value" != "31" && "$value" != "1" && "$value" != "0" && "$value" != "2" ]]; then
      input="$value"
    fi
  done
  cp "$input" "$output"
elif [[ "$1" == "filter" ]]; then
  output=
  input=
  for ((i=2; i <= $#; i++)); do
    value="${!i}"
    if [[ "$value" == "-o" ]]; then
      ((i++))
      output="${!i}"
    elif [[ "$value" != -* && "$value" != "1" && "$value" != "0" ]]; then
      if [[ -n "$input" ]]; then
        input="$value"
      else
        input=INDEX
      fi
    fi
  done
  if [[ "$input" == *output.fasta || "$input" == */candidate_baits.fasta ]]; then
    cp "$input" "$output"
  else
    : > "$output"
  fi
else
  printf 'unexpected fake deacon invocation: %s\n' "$*" >&2
  exit 1
fi
EOF
chmod +x "$fake_bin/deacon"

PATH="$fake_bin:$PATH" PYTHONPATH=src \
bash scripts/finalize_core_nt_validation.sh \
  --blast-tsv "$blast_tsv" \
  --input-baits "$input_baits" \
  --input-manifest "$input_manifest" \
  --output-baits "$output_baits" \
  --output-manifest "$output_manifest" \
  --decision-tsv "$decisions" \
  --near-hits-tsv "$near_hits" \
  --index "$index" \
  --background "$background_silva" \
  --background "$background_rfam" \
  --background "$background_other"

test "$(grep -c '^>' "$output_baits")" -eq 2
test "$(sed -n '1p' "$output_baits")" = ">bait_target|28S"
test "$(sed -n '3p' "$output_baits")" = ">bait_near|5S"
! grep -Fq '>bait_reject|18S' "$output_baits"

awk -F '\t' '
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      if ($i == "bait_id") bait = i
      if ($i == "status") status = i
    }
    next
  }
  {
    rows++
    seen[$bait]++
    if ($bait == "bait_reject|18S" &&
        $status == "REJECT_CORE_NT_EXACT_NON_TARGET") rejected = 1
  }
  END {
    if (rows != 3 || !rejected) exit 1
    for (bait_id in seen) if (seen[bait_id] != 1) exit 1
  }
' "$decisions"

grep -Fq $'bait_reject|18S\tREJECT_CORE_NT_EXACT_NON_TARGET' "$decisions"
grep -Fq $'bait_near|5S\tPASS_CORE_NT\t0\t0\t1' "$decisions"
grep -Fq $'bait_near|5S\tY.1\t1\tother' "$near_hits"

awk -F '\t' '
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      if ($i == "core_nt_status") core_status = i
      if ($i == "core_nt_rejection_reason") core_reason = i
      if ($i == "status") status = i
    }
    if (!core_status || !core_reason) exit 1
    next
  }
  $1 == "AACCGGTTACGATCGTAGCTAGGCTAACGTA" {
    if ($status != "REJECT_CORE_NT_EXACT_NON_TARGET" ||
        $core_status != "REJECT_CORE_NT_EXACT_NON_TARGET" ||
        $core_reason != "exact_non_target_31_of_31") exit 1
    rejected = 1
  }
  END {
    if (NR != 4 || !rejected) exit 1
  }
' "$output_manifest"

test -s "$index"
