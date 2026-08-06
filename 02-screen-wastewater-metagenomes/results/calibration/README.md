# Threshold calibration — why a read must carry 20 diagnostic 31-mers

This directory holds the evidence for the detection threshold: why a single
shared 31-mer is not enough to call a read *Cyclospora*, and why the screen runs
at `deacon filter -a 20 -r 0`. It rests on the public SRA runs alone, screened
permissively at `-a 1 -r 0` and then classified against `core_nt` along each
read's full length, so the calibration and [Figure 1](../figures/) are computed
from the same data.

## What was screened

205 SRA runs from BioProject PRJNA1247874 were screened successfully at `-a 1`:
**394.5 billion reads, 59.5 Tbp**. A further 88 runs failed to download from ENA
(aria2c timeouts and MD5 mismatches on multi-gigabyte FASTQs); they are listed
with their reason in [`sra_screen_a1_summary.tsv`](sra_screen_a1_summary.tsv)
rather than dropped, because a failed download is not a screening result.

160 runs returned at least one candidate sequence, 131,580 in all.

## The read population the sweep is computed over

Deacon's paired mode emits two kinds of sequence that must not reach a per-read
calibration, and removing them changes the numbers by an order of magnitude:

| | count |
|---|--:|
| sequences returned at `-a 1` | 131,580 |
| fragments | 65,790 |
| — duplicate fragments (both mates identical base for base) | −53,649 |
| — mates carrying no diagnostic 31-mer of their own | −7,857 |
| **reads carried into the sweep** | **16,425** |
| distinct sequences among them, the BLAST query | 8,421 |

**Duplicates.** A fragment that was amplified or resequenced is one observation,
not several. 82% of fragments here are duplicates. Two fragments count as one
only when *both* mates match base for base — a shared R1 with different R2s is
two molecules and stays two.

**Mates with nothing of their own.** Paired mode compares `-a` against the
*union* of distinct k-mer hits across both mates, so a read carrying no
diagnostic 31-mer at all is emitted whenever its partner matched. The `-a` rule
is per-read, so such a read can never pass any threshold; leaving it in would
understate what fraction of candidates are genuinely *Cyclospora*.

Fragments are collapsed **before** the empty mates are dropped. The other order
would make two distinct molecules that happen to share their k-mer-bearing mate
look like a single duplicate.

## What the candidates turn out to be

Each of the 8,421 distinct sequences was aligned along its full length against
`core_nt` and assigned by its best bit score anywhere in the collection:

| class | reads | |
|---|--:|---|
| `non_target` | 16,111 | best alignment is to something else |
| `target` | 282 | best alignment is *Cyclospora* |
| `top_tie` | 24 | the two tie — the evidence does not decide |
| `no_hit` | 8 | no core-nt alignment at all |

**98.1% of reads retained at `-a 1` are not *Cyclospora*.** They are dominated by
environmental rRNA — uncultured fungi (5,605 distinct sequences), the bdelloid
rotifer *Adineta vaga* (1,043), uncultured eukaryotes (333), *Paramecium*,
*Asplanchna* — and, importantly, by other apicomplexa: *Eimeria acervulina*
(380), *Voromonas pontica* (160), *Babesia microti* (66). A read can carry a
genuinely unique diagnostic 31-mer, one with no exact match outside *Cyclospora*
anywhere in core-nt, and still align better over its remaining 120 nt to a
relative. That is precisely why the threshold exists.

## The sweep

[`threshold_blast_read_counts.tsv`](threshold_blast_read_counts.tsv), reads
surviving each threshold by class:

| threshold | target | non-target | tie | no hit |
|--:|--:|--:|--:|--:|
| 1 | 282 | 16,111 | 24 | 8 |
| 2 | 275 | 2,683 | 24 | 0 |
| 3 | 263 | 491 | 14 | 0 |
| 4 | 259 | 11 | 2 | 0 |
| 6 | 248 | 3 | 2 | 0 |
| 10 | 227 | 2 | 0 | 0 |
| 12 | 221 | 2 | 0 | 0 |
| **13** | **220** | **0** | **0** | **0** |
| 20 | 187 | 0 | 0 | 0 |

**The empirical noise ceiling is 12.** The highest count on any read that is not
*Cyclospora* is 12, carried by one sequence found in two fragments of
`MO_Monett_20250610__SRR35939776`. It is not an artifact: it aligns to *Eimeria
acervulina* 28S at 97.2% identity and bit score 241, against 95.9% and 237 to
*C. cayetanensis* — a genuine near-relative apicomplexan that beats the target by
four bits. Ambiguous reads are gone by 10, and **13 is the lowest threshold at
which nothing survives that is better explained by another organism.**

The screen runs at **20**, which is 1.5× the empirical minimum. It was fixed
before this calibration and is left where it is: it costs 33 of 282 target reads
(220 → 187) and buys a margin against the next *Eimeria*-like read in data not
yet screened. Genuine *Cyclospora* reads carry a median of 34 diagnostic 31-mers
and up to 105, so the threshold sits far below the bulk of real signal.

## Target is the genus, not the species

A read is *target* when its best alignment carries a taxid in the genus
*Cyclospora* — the 26 taxids in
[`cyclospora_genus_taxids.txt`](cyclospora_genus_taxids.txt), expanded from the
NCBI taxonomy that ships with core_nt. GenBank holds genus-level deposits such as
U40261.1 "*Cyclospora* sp." whose 18S is identical to *C. cayetanensis*; under a
species-strict rule a read matching both scores as a tie and looks like a
specificity failure, when all that happened is that a depositor did not name a
species. The genus is the meaningful unit for "is this *Cyclospora*".

On this dataset the two rules agree exactly — `--target-taxids 88456` reproduces
the same 282 / 16,111 / 24 / 8 split and the same ceiling of 12 — so the choice
changes nothing here. It is stated because it did matter on an earlier, smaller
read set, where species-strict moved the apparent clean threshold from 7 to 26
without a single genuinely non-*Cyclospora* read being involved.

## Files

| File | What it holds |
|---|---|
| [`sra_screen_a1_summary.tsv`](sra_screen_a1_summary.tsv) | one row per SRA run screened at `-a 1`: depth, reads retained, wall time, and the reason for each failure |
| [`candidate_reads_a1.fasta.gz`](candidate_reads_a1.fasta.gz) | the 8,421 distinct candidate sequences — the BLAST query, so the classification can be repeated without re-screening 59.5 Tbp |
| [`read_blast_map.tsv`](read_blast_map.tsv) | every kept read, its diagnostic 31-mer count, and the distinct sequence that represents it |
| [`read_blast_deacon.tsv`](read_blast_deacon.tsv) | one row per read: its 31-mer count and its core-nt class — the evidence the sweep reads |
| [`cyclospora_genus_taxids.txt`](cyclospora_genus_taxids.txt) | the taxids that count as target |
| [`threshold_blast_read_counts.tsv`](threshold_blast_read_counts.tsv) | reads surviving each threshold, by class |
| [`absolute_threshold_curve.tsv`](absolute_threshold_curve.tsv) | runs and fragments surviving each threshold |

## Rebuilding it

The sweep regenerates from committed evidence with no database, cluster, or
network:

```bash
python3 ../../scripts/sweep_threshold.py
```

Redoing the classification itself needs the ~285 GB nucleotide database
(`update_blastdb.pl --decompress core_nt`), but not the sequencing data:

```bash
gunzip -c candidate_reads_a1.fasta.gz > query.fasta
blastn -task blastn -db core_nt -query query.fasta \
  -evalue 1e-10 -max_target_seqs 25 -dust no \
  -outfmt '6 qseqid qlen saccver staxids pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp stitle' \
  -out read_blast.tsv
python3 ../../scripts/classify_reads.py --blast read_blast.tsv
python3 ../../scripts/sweep_threshold.py
```

Starting from the SRA runs themselves — the only step that needs the sequencing
data — screen each run listed in `sra_screen_a1_summary.tsv` and rebuild the
query:

```bash
deacon filter -m -a 1 -r 0 cyclospora_k31w1.idx R1.fastq.gz R2.fastq.gz \
  --fasta -o <run>.diagnostic_pairs.fasta.gz
python3 ../../scripts/prepare_read_blast_query.py --reads <dir> \
  --query-output query.fasta --map-output read_blast_map.tsv
```

`prepare_read_blast_query.py` recounts every read against the stage 01 bait FASTA
rather than trusting the screen, so it needs no Deacon and no network. Its counts
reproduce `deacon filter --debug` exactly on all 16,425 reads.
