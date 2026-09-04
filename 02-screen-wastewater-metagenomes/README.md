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
  import_nvd_screen.py          freezes primary public outputs from a returned NVD screen
  summarize_nvd_screen.py       recounts the frozen screen and writes results
  plot_heatmap.py               builds Figure 1 and the matrix behind it
  sweep_threshold.py            re-derives the threshold table
  summarize_bait_calibration.py derives the bait and calibration summary
  prepare_read_blast_query.py   dedups candidate reads to the unique BLAST query
  classify_reads.py             assigns each read target/non-target from core-nt
results/
  screen-source/                       lossless primary public reports and retained FASTQs
  sra_sample_summary.tsv               per-run results, all 2,333 screened SRA runs
  casper_sites.tsv                     the SRA sewershed codes, with coordinates and role
  site_fortnight_matrix_per_billion.tsv  the matrix plotted in Figure 1
  diagnostic_reads.tar.gz              one FASTA member per positive SRA run (81 runs)
  figures/cyclospora_heatmap.svg       static figure, embedded in the main README
  figures/cyclospora_heatmap.vl.json   Vega-Lite spec behind the interactive figure
  figures/cyclospora_heatmap.html      vega-embed wrapper around that spec
  calibration/                         threshold evidence — the -a 1 calibration and its core-nt classification
    bait_calibration_summary.json      derived bait topology and calibration measurements
pixi.toml, pixi.lock            Deacon and BLAST+, only needed to screen your own reads
```

The import, summary, and plotting scripts use the standard library only. Summarization and plotting
run from committed inputs without the returned NVD directory, a network connection, or a database.

## The short version

A read counts as *Cyclospora* when it carries **at least 24 diagnostic 31-mers of its own**. In
Deacon that is `-a 24 -r 0`. Non-target or tied reads remained through a bait count of 23; 24 is
the lowest threshold at which only target-classified reads remain. See
[the threshold section](#the-threshold-evidence) below.

## Screening results

NVD 3.5.0 (revision `bd421147bba2a8ba22aafffb5dbf61c47d241feb`, clean) screened 2,333 public runs with Deacon 0.16.0 at `-a 24 -r 0`. The committed source bundle
preserves every public Deacon report and retained FASTQ. The summary then recounts each read on its
own, leaving 1,244 diagnostic reads across 81 positive runs. Within runs, those collapse to 379
distinct sequences when reverse complements are treated as identical. Restricting to the 30
sewersheds sampled at 10 or more timepoints leaves 2,328 runs and 4.32 trillion reads. The signal is
seasonal and recurs through the summers of both 2025 and 2026.

```bash
pixi run summarize-nvd-screen  # rebuilds the summary and published reads
pixi run plot-heatmap           # rebuilds the figure and matrix from that summary
pixi run rebuild-screen-artifacts  # performs both steps in order
```

The one-time import from the returned NVD directory is:

```bash
pixi run import-nvd-screen --run-dir /path/to/returned/run --replace
```

The import boundary is NVD's primary per-sample reports and retained FASTQs. The returned directory's
`summaries/postmerge/` recounts are downstream audit work and are not imported.

The value plotted is **distinct diagnostic reads per billion reads sequenced**. Distinct here is the
summary's count of unique read sequences per run, treating reverse complements as identical, so PCR
and optical copies are collapsed before anything is pooled; the raw retained count is in the summary
alongside. A cell covering more than one run pools summed reads over summed depth, never a mean of
per-run rates. Pass `--raw` to plot the undeduplicated counts and `--min-timepoints N` to vary the
inclusion rule.

### The reads for Figure 1

Every read counted in the summary and plotted in Figure 1 is committed in
[`results/diagnostic_reads.tar.gz`](results/diagnostic_reads.tar.gz) as **one FASTA member per run**,
named `<CODE>_<YYYYMMDD>__<accession>` for the sewershed's SRA code and the run it came from — 1,244
reads across 81 runs, each header carrying its diagnostic 31-mer count.
Runs that were screened and yielded nothing have no archive member; that a run was screened and came
back clean is recorded in [`sra_sample_summary.tsv`](results/sra_sample_summary.tsv), which covers all
2,333.

```bash
pixi run summarize-nvd-screen
```

This verifies the frozen archive against its manifest, recounts every retained read, and writes both
the read archive and per-run summary. Deacon pools k-mer hits across mates in paired mode, so the
summary counts only reads that reach 24 on their own.

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

The calibration contains 108,474 reads carrying at least one of the 1,464 validated baits. The
10,894 distinct sequences are committed as
[`calibration_read_blast_queries.fasta.gz`](results/calibration/calibration_read_blast_queries.fasta.gz)
and were submitted in full to BLASTN against `core_nt`. Classification compared their best local-alignment bit scores within and outside the declared *Cyclospora* target scope. The per-read evidence and threshold table are in
[`results/calibration/`](results/calibration/), which explains the analysis in full.

At `-a 1`, the calibration contains 1,284 target reads, 107,040 non-target reads, 123 ties, and 27
reads with no `core_nt` hit. At 23 baits, 855 target reads and 16 ties remain. At **24**, 853 target
reads remain and no non-target, tied, or no-hit reads remain, making 24 the lowest threshold
supported by the calibration rule.

The threshold table regenerates from committed evidence, with no database, cluster, or network:

```bash
python3 scripts/sweep_threshold.py     # -> threshold_read_counts.tsv
```

The bait distribution, spatial runs, calibration-read and bait-utilization distributions, target-alignment support, and sensitivity to baits admitted without an exact `core_nt` match regenerate from the same committed evidence:

```bash
pixi run summarize-bait-calibration
# -> results/calibration/bait_calibration_summary.json
```

The task independently recounts each calibration sequence against the bait FASTA and fails if that count disagrees with the committed read evidence. It uses the archived BLAST output rather than repeating the search, so it requires no database or network access. The resulting score margins compare the best reported target and non-target alignments; the archived search used `-max_target_seqs 100`, and the summary records how many threshold-passing queries reached that reporting limit. This task summarizes calibration evidence only—it does not reproduce the BioProject-wide production screen.

The narrative that interprets the sweep is in the
[Results](../README.md#setting-a-calibration-threshold-of-24-diagnostic-31-mers-before-a-read-counts-as-cyclospora) of
the main README.

The threshold sweep uses each read's own count of distinct diagnostic 31-mers. This avoids treating
the pooled count across paired mates as if both individual reads reached the threshold.

### Repeating the read classification

The committed raw alignments reproduce the classification without the database:

```bash
gzip -dc results/calibration/read_blast_hits.tsv.gz > read_blast.tsv
python3 scripts/classify_reads.py --blast read_blast.tsv
```

Repeating the BLASTN search needs the database, roughly 285 GB, retrieved with
`update_blastdb.pl --decompress core_nt`. Each complete read sequence is submitted as a query, but
BLASTN reports local alignments and no minimum query coverage is imposed.

```bash
gzip -dc results/calibration/calibration_read_blast_queries.fasta.gz > query.fasta
pixi run blastn -task blastn -db core_nt -query query.fasta \
  -evalue 1e-10 -max_target_seqs 100 -dust no \
  -outfmt '6 qseqid qlen saccver staxids pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp stitle' \
  -out read_blast.tsv
```

A read is **target** when its highest local-alignment bit score is to a taxon in the declared 26-taxid
*Cyclospora* scope, **non-target** when the highest score is outside that scope, and a **tie** when the
best scores within and outside the scope differ by less than 0.1 bits. Ties are held apart rather than
assigned because the alignment evidence does not decide.

## Screening your own reads

Build the index once from the stage 01 baits, then filter.

```bash
deacon index build -k 31 -w 1 -e 0 \
  ../01-identify-cyclospora-specific-kmers/baits/cyclospora_cayetanensis_rrna_core_nt_validated_baits.fasta \
  -o cyclospora_k31w1.idx

deacon filter -a 24 -r 0 cyclospora_k31w1.idx reads_R1.fastq.gz reads_R2.fastq.gz
```

Two things will bite you if you skip them.

**Deacon pools k-mer hits across mates in paired mode.** A pair whose mates carry 13 and 12 *disjoint*
hits passes `-a 24`, because the union is 25, even though neither read reaches 24 on its own. Each
retained read must therefore be recounted against the bait FASTA after filtering, and only reads
reaching 24 by themselves should be counted.

**The threshold was calibrated on 150–151 nt reads.** Those reads have 120–121 possible 31-mer
positions, so demanding 24 is demanding roughly a fifth of them. On 100 nt reads the same number is
a far harsher demand, and on 250 nt reads a far softer one. Repeat the calibration if your read
lengths differ substantially.
