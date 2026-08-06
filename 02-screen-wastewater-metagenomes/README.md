# Stage 02. Screen wastewater metagenomes

Applies the bait set from [stage 01](../01-identify-cyclospora-specific-kmers/) to the public
wastewater sequencing in NCBI SRA BioProject
[PRJNA1247874](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1247874).

The detection threshold, the screening results across 30 sewersheds, and the figure built from them
are all here. Every sewershed is named by its `public_code_CASPER_SRA` — the identity it is deposited
under in SRA — so nothing here depends on an internal sample name.

## What is here

```text
scripts/
  plot_heatmap.py               builds Figure 1 and the matrix behind it
  verify_published_reads.py     recounts every published read against the baits
  sweep_threshold.py            re-derives the two threshold sweep tables
  recount_diagnostic_reads.py   recounts each retained read against the bait FASTA
  prepare_read_blast_query.py   dedups candidate reads to the unique BLAST query
  classify_reads.py             assigns each read target/non-target from core-nt
results/
  sra_sample_summary.tsv               per-run results, all 2,292 screened SRA runs
  casper_sites.tsv                     the SRA sewershed codes, with coordinates and role
  site_fortnight_matrix_per_billion.tsv  the matrix plotted in Figure 1
  reads/                               one FASTA per SRA run with diagnostic reads (73 runs)
  figures/cyclospora_heatmap.svg       static figure, embedded in the main README
  figures/cyclospora_heatmap.vl.json   Vega-Lite spec behind the interactive figure
  figures/cyclospora_heatmap.html      vega-embed wrapper around that spec
  calibration/                         threshold evidence — the -a 1 SRA screen and its core-nt classification
pixi.toml, pixi.lock            Deacon and BLAST+, only needed to screen your own reads
```

`python3 scripts/plot_heatmap.py` and `python3 scripts/verify_published_reads.py` use the standard
library only and run on a bare checkout with no environment, no network, and no database.

## The short version

A read counts as *Cyclospora* when it carries **at least 20 diagnostic 31-mers of its own**. In
Deacon that is `-a 20 -r 0`. Requiring fewer admits reads that belong to other organisms: calibrated
against the SRA runs alone, the highest count on any read that is not *Cyclospora* is 12, so 13 is
the lowest fully specific threshold and 20 is a conservative 1.5× margin above it. See
[the threshold section](#the-threshold-evidence) below.

## Screening results

Applying `deacon filter -a 20 -r 0` to every run in the BioProject, then recounting each retained
read on its own, gives the results in
[`results/sra_sample_summary.tsv`](results/sra_sample_summary.tsv) — 2,292 runs collected between
2023-12-26 and 2026-06-30. Restricting to the 30 sewersheds sampled at 10 or more timepoints leaves
2,287 runs and 4.20 trillion reads, of which **351 distinct diagnostic reads** (1,156 before
collapsing duplicates) met the threshold across 73 positive runs. The signal is seasonal and recurs:
detection rises through the summers of both 2025 and 2026.

```bash
python3 scripts/plot_heatmap.py     # rebuilds the figure and the matrix
```

The value plotted is **distinct diagnostic reads per billion reads sequenced**. Distinct here is
Deacon's own count of unique read sequences per run, so PCR and optical copies are collapsed before
anything is pooled; the raw retained count is in the summary alongside. A cell covering more than one
run pools summed reads over summed depth, never a mean of per-run rates. Pass `--raw` to plot the
undeduplicated counts and `--min-timepoints N` to vary the inclusion rule.

### The reads for Figure 1

Every read counted in the summary and plotted in Figure 1 is committed under
[`results/reads/`](results/reads/) as **one gzipped FASTA per run**, named
`<CODE>_<YYYYMMDD>__<accession>` for the sewershed's SRA code and the run it came from — 1,156 reads
across 73 runs ([index](results/reads/README.md)), each header carrying its diagnostic 31-mer count.
Runs that were screened and yielded nothing have no file; that a run was screened and came back clean
is recorded in [`sra_sample_summary.tsv`](results/sra_sample_summary.tsv), which covers all 2,292.

```bash
python3 scripts/verify_published_reads.py    # recount them all against the bait set
```

This confirms that every published read reaches 20 diagnostic 31-mers on its own, that each header
count matches an independent recount, and that per-run totals agree with the summary. It matters
because Deacon pools k-mer hits across mates in paired mode, so a pair can clear `-a 20` without
either read reaching 20 alone. The minimum across all 1,156 is exactly 20.

The static figure is [`results/figures/cyclospora_heatmap.svg`](results/figures/cyclospora_heatmap.svg)
and the interactive one is a
[Vega-Lite specification](results/figures/cyclospora_heatmap.vl.json) with a
[wrapper page](results/figures/cyclospora_heatmap.html) that embeds it, reporting per-cell read
counts, depth, and contributing runs on hover. Open the wrapper from a local clone or serve it from
Pages; GitHub serves committed HTML as source rather than rendering it. The specification is the
portable artifact: it also renders in VS Code, JupyterLab, Observable, and the Vega editor, and it
carries its data inline, so nothing else needs to be fetched.

Columns are ordered west to east from the coordinates in
[`results/casper_sites.tsv`](results/casper_sites.tsv), which are approximate city centroids used only
to order the figure, not survey coordinates.

## The threshold evidence

The calibration rests on the public SRA runs alone: 205 of them screened at `-a 1 -r 0` — 394.5
billion reads, 59.5 Tbp — and every candidate read then aligned along its full length against
`core_nt`. The evidence, the sweep tables, and the candidate reads themselves are in
[`results/calibration/`](results/calibration/), which explains the analysis in full.

At `-a 1`, **98.1% of retained reads are not *Cyclospora***: uncultured fungi, bdelloid rotifers,
and other apicomplexa including *Eimeria* and *Babesia*. The highest diagnostic 31-mer count carried
by any read that is not *Cyclospora* is **12** — an *Eimeria acervulina* 28S read that beats
*C. cayetanensis* by four bits along its length — so **13 is the lowest fully specific threshold**.
The screen runs at **20**, 1.5× that minimum. Of the 282 target reads found at `-a 1`, 220 survive
the fully specific threshold of 13 and 187 survive 20, so the extra margin costs 33 of them and buys
room against the next such read in data not yet screened.

Both sweep tables regenerate from committed evidence, with no database, cluster, or network:

```bash
python3 scripts/sweep_threshold.py     # -> threshold_blast_read_counts.tsv, absolute_threshold_curve.tsv
```

The narrative that interprets the sweep is in the
[Results](../README.md#setting-a-calibration-threshold-of-twenty-diagnostic-31-mers-before-a-read-counts-as-cyclospora) of
the main README.

The reads carried into the sweep are not simply what Deacon returned. Paired mode pools k-mer hits
across mates and emits both, so a read carrying no diagnostic 31-mer of its own rides out on its
partner, and amplified fragments are returned many times over. Collapsing duplicate fragments — a
duplicate only when *both* mates match base for base — and then dropping the mates that carry
nothing takes 131,580 returned sequences to **16,425 reads**, of which 8,421 are distinct. Both
corrections are applied by `scripts/prepare_read_blast_query.py`, which recounts each read against
the stage 01 baits rather than trusting the screen.

### Repeating the whole-read alignment

The classification behind the sweep aligns each distinct candidate read along its full length against
`core_nt` and takes its global best bit score. The 8,421 distinct reads are committed as
[`results/calibration/candidate_reads_a1.fasta.gz`](results/calibration/candidate_reads_a1.fasta.gz),
so this can be repeated without re-screening 59.5 Tbp. It needs the database, roughly 285 GB,
retrieved with `update_blastdb.pl --decompress core_nt`.

```bash
gunzip -c results/calibration/candidate_reads_a1.fasta.gz > query.fasta
pixi run blastn -task blastn -db core_nt -query query.fasta \
  -evalue 1e-10 -max_target_seqs 25 -dust no \
  -outfmt '6 qseqid qlen saccver staxids pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp stitle' \
  -out read_blast.tsv
python3 scripts/classify_reads.py --blast read_blast.tsv    # -> read_blast_deacon.tsv
```

A read is **target** when its single highest-scoring alignment anywhere in the collection is to the
genus *Cyclospora*, **non-target** when it is to anything else, and a **tie** when the top bit score
is shared between the two. Ties are held apart rather than assigned, because a tie is exactly the
case where the evidence does not decide. Target is the genus rather than *C. cayetanensis* alone
because GenBank carries genus-level deposits with identical 18S; the taxids are committed in
[`results/calibration/cyclospora_genus_taxids.txt`](results/calibration/cyclospora_genus_taxids.txt),
and `--target-taxids 88456` reproduces the species-strict variant.

## Screening your own reads

Build the index once from the stage 01 baits, then filter.

```bash
deacon index build -k 31 -w 1 -e 0 \
  ../01-identify-cyclospora-specific-kmers/baits/cyclospora_cayetanensis_rrna_core_nt_validated_baits.fasta \
  -o cyclospora_k31w1.idx

deacon filter -m -a 20 -r 0 cyclospora_k31w1.idx reads_R1.fastq.gz reads_R2.fastq.gz
```

Two things will bite you if you skip them.

**Deacon pools k-mer hits across mates in paired mode.** A pair whose mates carry 13 and 10 *disjoint*
hits passes `-a 20`, because the union is 23, even though neither read reaches 20 on its own. Every
read reported in this work was therefore recounted against the bait FASTA after filtering, and only
reads reaching 20 by themselves are counted or published. This is why our numbers are lower than
Deacon's retained-read counts.

**The threshold is specific to 151 nt reads.** A 151 nt read has 121 possible 31-mer positions, so
demanding 20 is demanding roughly a sixth of them. On 100 nt reads the same number is a far harsher
demand, and on 250 nt reads a far softer one. Repeat the calibration if your read lengths differ.