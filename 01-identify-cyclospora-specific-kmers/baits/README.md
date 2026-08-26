# *Cyclospora cayetanensis* ribosomal RNA bait k-mers

Two FASTA files, each record a single canonical 31-mer drawn from *C. cayetanensis* mature ribosomal
RNA. **Use the core-nt-validated set to reproduce the downstream results in this repository.**

| File | Records | Stage | Use it when |
|---|---:|---|---|
| [cyclospora_cayetanensis_rrna_core_nt_validated_baits.fasta](cyclospora_cayetanensis_rrna_core_nt_validated_baits.fasta) | 1,464 | After exact `core_nt` screening | You want to reproduce the downstream analysis. |
| [cyclospora_cayetanensis_rrna_specific_baits.fasta](cyclospora_cayetanensis_rrna_specific_baits.fasta) | 1,670 | Before exact `core_nt` screening | You want to inspect or repeat taxonomic screening. |

Composition of the validated set is 61 × 18S, 1,372 × 28S, and 31 × 5S. There are no 5.8S baits,
because all 126 5.8S candidates occur in Rfam RF00002 and were subtracted.

## Record format

```
>cc_rrna_kmer_000001|28S
AAAAACACGAACCTCTCCCTACTCTCACTCT
```

The identifier is a sequential number within the pre-screening file followed by the rRNA class the k-mer came from. It is not stable when the upstream bait set changes. Sequences
are canonical, meaning each is the lexicographically smaller of the 31-mer and its reverse
complement, so a bait matches a read on either strand. Because they are canonical rather than
genomic, the records are **not** in coordinate order along the rRNA and adjacent records do not
overlap. To recover genomic context, join on the k-mer column of
[the full k-mer manifest](../results/cyclospora_cayetanensis_rrna_kmer_manifest.tsv),
which carries the source locus, start coordinates, and copy count for every candidate.

## What these are and are not

**What they are.** Each 31-mer in the 1,464-bait screened set occurs in *C. cayetanensis* mature rRNA, occurs in no other rRNA
sequence in SILVA 138.2 or Rfam 15.1, occurs in no non-*cayetanensis Cyclospora* reference we could
obtain, and has no exact full-length match anywhere in NCBI core-nt assigned to a taxon other than
*C. cayetanensis* (taxid 88456).

**What they are not.** They are not a genus-level assay, they are not a quantitative assay, and they
do not distinguish *C. cayetanensis* from *C. ashfordi* or *C. henanensis* beyond what the available
reference sequence supports. Detection of these k-mers in a sample is evidence of
*C. cayetanensis* rRNA sequence, not of viable oocysts or of infectivity.

## Using them with Deacon

```bash
# Build the filtering index. k=31 and w=1 are not optional; see below.
deacon index build -k 31 -w 1 -e 0 \
  cyclospora_cayetanensis_rrna_core_nt_validated_baits.fasta \
  -o cyclospora_cayetanensis_rrna_core_nt_validated_k31w1.idx

# Retain reads carrying at least 24 distinct diagnostic 31-mers.
deacon filter -a 24 -r 0 \
  cyclospora_cayetanensis_rrna_core_nt_validated_k31w1.idx \
  reads_R1.fastq.gz reads_R2.fastq.gz
```

`w=1` makes every 31-mer its own minimizer, so nothing is subsampled and the `-a` threshold reads
directly as "number of distinct diagnostic 31-mers found in this read". At any larger `w` the
threshold stops meaning that. The index rebuilds byte-identically from the FASTA and has SHA-256
`0a695f65d973a067c690699e3025da50519ab0de9be73b8f43606a763fdf520d`.

> **Deacon pools k-mer hits across mates in paired mode.** Two mates carrying 13 and 12 *disjoint*
> hits are retained at `-a 24`, because the union is 25. If you need a per-read guarantee, recount
> each retained read against the bait FASTA independently rather than trusting the pair-level
> decision. This repository does exactly that, and it is the reason the reported counts are lower
> than Deacon's own retained-read counts.

## Using them without Deacon

Nothing about the bait set is Deacon-specific. It is a list of 31-mers. Any exact k-mer matcher
works, and the counting rule that matters is "how many *distinct* baits does this read contain",
counting both strands. See
[`src/rrna_bait/core.py`](../src/rrna_bait/core.py) for the canonical-form
helper used here.

## Provenance and how to rebuild

Every input is pinned by release, URL, byte size, and SHA-256 in
[`config/sources.tsv`](../config/sources.tsv) and
[`results/input_manifest.tsv`](../results/input_manifest.tsv). The build
verifies those identities before use and fails rather than proceeding on a changed input. Full
method, attrition at every stage, and rebuild instructions are in the
[main README](../../REPRODUCING.md).

## Citation and license

If you use these baits, please cite this repository and the *Cyclospora* reference sequences the
baits were derived from, which are listed in
[`config/target_rrna_accessions.txt`](../config/target_rrna_accessions.txt).
The k-mer sequences themselves are derived from public NCBI, SILVA, and Rfam records.
