# Reproducing this analysis

This document is the methodological companion to the [main README](README.md), which reports the findings.

## Repository organization

Each analysis stage is a numbered, self-contained directory with its own environment, code, inputs, committed results, and README.

```text
01-identify-cyclospora-specific-kmers/   design and validate the bait set
  config/     pinned inputs, accession lists, alignment thresholds
  scripts/    the three entry points
  src/        the Python behind them
  tests/      offline unit and fixture tests
  curated/    mature rRNA loci located in the assemblies
  results/    committed outputs of our run
  baits/      the published bait FASTAs
02-screen-wastewater-metagenomes/        apply the baits to wastewater libraries
  scripts/    threshold sweep and figure generation
  results/    calibration evidence, per-sample screen results, figures
```

The stage 01 README, [`01-identify-cyclospora-specific-kmers/README.md`](01-identify-cyclospora-specific-kmers/README.md), describes the contents of that directory path by path. The steps below are run from within it.

## Prerequisites

[Pixi](https://pixi.sh) is the only prerequisite. It installs every tool at the version used here from [`pixi.toml`](01-identify-cyclospora-specific-kmers/pixi.toml) and the accompanying `pixi.lock`: Python 3.12, BLAST+ 2.15, meryl 1.4.1, Deacon 0.15, seqkit 2.13, barrnap 0.9, and entrez-direct 25.3. The lock is solved for `linux-64`; on another platform, add it to `platforms` in `pixi.toml` and re-solve.

Requirements are approximately 330 MB of network downloads and a few GB of working space for intermediate k-mer databases.

```bash
git clone <this repository>
cd cyclospora-in-wastewater-metagenomics/01-identify-cyclospora-specific-kmers
pixi install
pixi run test    # 89 unit tests and 3 fixture builds; no network, no prior build
```

## Step 1. Input data

The build retrieves each input automatically.

| Input | Pinned version | Role |
|---|---|---|
| [SILVA SSURef_NR99](https://ftp.arb-silva.de/release_138.2/Exports/) and LSURef_NR99 | 138.2 | rRNA background to subtract |
| [Rfam RF00001](https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/fasta_files/) (5S) and RF00002 (5.8S) | 15.1 | rRNA background to subtract |
| *C. cayetanensis* nuclear assemblies [GCF_002999335.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002999335.1/) (CcayRef3) and [GCF_000769155.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000769155.1/) (ASM76915v2) | RefSeq | the target the k-mers come from |
| Nine reference transcripts: `AF111183.1` (18S), `MPGL01000046.1` (28S), and seven `XR_` records (5S, 5.8S) | GenBank, retrieved 2026-07-23 | queries used to locate the mature loci |
| 94 non-*cayetanensis Cyclospora* records | GenBank, versioned accessions, retrieved 2026-07-23 | background from the parasite's close relatives |

Every input is pinned by URL, byte count, and SHA-256 in [`config/sources.tsv`](01-identify-cyclospora-specific-kmers/config/sources.tsv), and the accession lists are in [`config/target_queries.tsv`](01-identify-cyclospora-specific-kmers/config/target_queries.tsv), [`config/target_rrna_accessions.txt`](01-identify-cyclospora-specific-kmers/config/target_rrna_accessions.txt), and [`config/other_cyclospora_accessions.txt`](01-identify-cyclospora-specific-kmers/config/other_cyclospora_accessions.txt). Retrieval and checksum verification are handled by [`src/rrna_bait/sources.py`](01-identify-cyclospora-specific-kmers/src/rrna_bait/sources.py), which verifies each identity before use and halts on a changed input rather than proceeding. The files downloaded during our own run, and their timestamps, are recorded in [`results/input_manifest.tsv`](01-identify-cyclospora-specific-kmers/results/input_manifest.tsv). Downloads are cached and reused; setting `RRNA_BAIT_OFFLINE=1` disables network access and requires a verified cache.

## Step 2. Build the bait set

```bash
pixi run build          # scripts/build_rrna_bait.sh; RRNA_BAIT_THREADS=8 by default
```

This command runs [`scripts/build_rrna_bait.sh`](01-identify-cyclospora-specific-kmers/scripts/build_rrna_bait.sh), which executes the design end to end in a single pass: it retrieves and checksum-verifies the Step 1 inputs, removes the target's own accessions from the SILVA and Rfam backgrounds so the target cannot subtract itself, and then carries out the five stages below. It writes the located loci to `curated/`, the bait FASTA and per-k-mer manifest to `kmers/`, the index `cyclospora_cayetanensis_rrna_k31w1.idx` to the stage directory, the verification summary and report to `reports/`, the record of what was downloaded to `provenance/`, and intermediates to `work/` and `logs/`. Any failed check aborts the run with a nonzero exit status.

The parameters at each stage are as follows.

1. **Locate the mature loci.** The nine reference transcripts are aligned to the two assemblies with `blastn -task megablast`, and a hit is kept only at ≥98% identity and ≥90% query coverage. The per-query thresholds and the 18S and 28S trim coordinates are in [`config/target_queries.tsv`](01-identify-cyclospora-specific-kmers/config/target_queries.tsv), and the selection logic in [`src/rrna_bait/targets.py`](01-identify-cyclospora-specific-kmers/src/rrna_bait/targets.py). Internal transcribed spacers, organelle rRNA, and genomic flanks are excluded at this stage, so only mature rRNA reaches the k-mer step. Output: [`curated/target_loci.tsv`](01-identify-cyclospora-specific-kmers/curated/target_loci.tsv) and [`curated/target_rrna.fasta`](01-identify-cyclospora-specific-kmers/curated/target_rrna.fasta), 83 alignment rows collapsing to 7 distinct sequences in our run.
2. **Enumerate 31-mers.** `meryl count k=31` over the target loci and over each background separately.
3. **Subtract.** `meryl union` of the backgrounds, then `meryl difference` against the target. Subtracting SILVA and Rfam alone yields the 1,839-member genus-compatible set; adding the other-*Cyclospora* records yields the species-specific set. Each other-*Cyclospora* record is subtracted whole rather than trimmed to annotated rRNA, which can only discard baits, never admit them. Per-k-mer bookkeeping is assembled into the manifest by [`src/rrna_bait/manifest.py`](01-identify-cyclospora-specific-kmers/src/rrna_bait/manifest.py).
4. **Drop low-complexity k-mers.** `deacon index build -k 31 -w 1 -e 0.6`, an entropy threshold of 0.6.
5. **Build and check the index.** `deacon index build -k 31 -w 1 -e 0`, followed by two round trips that must both hold: every bait is recovered when the bait FASTA is filtered against its own index, and zero records are retained when each background FASTA is filtered against it. The full invariant list is in [`src/rrna_bait/verify.py`](01-identify-cyclospora-specific-kmers/src/rrna_bait/verify.py).

Expected output is 1,388 baits in `kmers/cyclospora_cayetanensis_rrna_baits.fasta` and the attrition of Table 1 in `reports/cyclospora_cayetanensis_rrna_index_summary.tsv`, matching the committed [`results/cyclospora_cayetanensis_rrna_index_summary.tsv`](01-identify-cyclospora-specific-kmers/results/cyclospora_cayetanensis_rrna_index_summary.tsv).

## Step 3. Validate against the complete nucleotide collection

This step reduces the 1,388 baits from Step 2 to the 1,184 published baits. Each bait is searched against NCBI `core_nt` for exact, full-length 31-of-31 matches, and a bait is discarded if any such match is assigned to a taxon other than 88456 (*C. cayetanensis*); a match carrying no taxid counts as non-target. Of the 1,388 baits, 204 have at least one non-target match and are discarded, leaving 1,184.

The search results are committed, so the reduction reproduces in a few minutes without the database.

```bash
pixi run validate       # applies the committed BLAST evidence to the baits from Step 2
```

The evidence is [`results/core_nt_bait_exact_match_blast.tsv`](01-identify-cyclospora-specific-kmers/results/core_nt_bait_exact_match_blast.tsv), which holds one row per bait–subject match: 5,012 rows, because a bait present in many database records produces a row for each. Those rows cover the 1,280 baits with at least one exact match anywhere in `core_nt`; the other 108 match nothing and are retained.

The rule is implemented in [`src/rrna_bait/core_nt.py`](01-identify-cyclospora-specific-kmers/src/rrna_bait/core_nt.py), and the per-bait verdicts, with target and non-target hit counts, are written to [`results/cyclospora_cayetanensis_core_nt_validation.tsv`](01-identify-cyclospora-specific-kmers/results/cyclospora_cayetanensis_core_nt_validation.tsv). [`scripts/finalize_core_nt_validation.sh`](01-identify-cyclospora-specific-kmers/scripts/finalize_core_nt_validation.sh) re-runs the index and background round trips from Step 2 before publishing any output.

### Repeating the core-nt search

This requires the `core_nt` BLAST database, approximately 285 GB across 89 volumes, retrieved with `update_blastdb.pl --decompress core_nt` from the BLAST+ installation Pixi provides. The search is a single command over all 1,388 baits.

```bash
blastn -task blastn -word_size 31 -ungapped \
       -perc_identity 100 -qcov_hsp_perc 100 -dust no \
       -db core_nt \
       -query kmers/cyclospora_cayetanensis_rrna_baits.fasta \
       -outfmt '6 qseqid saccver staxids sscinames pident length mismatch gapopen qstart qend sstart send stitle' \
       -out core_nt_blast.tsv

bash scripts/build_core_nt_validated_index.sh core_nt_blast.tsv
```

We ran this query as one job per database volume on an HTCondor pool, staging each volume to node-local scratch, then concatenated the per-volume TSVs into the committed evidence file. The submit files are no longer carried in the repository. The result is determined by the query, its parameters, and the decision rule above, which are the same on a single machine.

## Step 4. Check the result

These three commands, run from the stage directory after Step 3, confirm that the rebuild matches ours. They use the build's own file names; the same two bait sets are published under `baits/` with longer names.

```bash
grep -c '^>' kmers/cyclospora_cayetanensis_rrna_baits.fasta                     # 1388, after Step 2
grep -c '^>' kmers/cyclospora_cayetanensis_rrna_core_nt_validated_baits.fasta   # 1184, after Step 3
shasum -a 256 cyclospora_cayetanensis_rrna_core_nt_validated_k31w1.idx
# 4bd2ee592ab7dfff30b56bfebd8346f7b2b91e903d1f8a88639d8e19b0d8e248
```

The checksum is the decisive one: it matches only if the 1,184 baits are identical to ours. `git diff curated/` should also stay empty, since `curated/` is committed and overwritten in place, which checks the locus-finding step. Everything else a rerun writes is gitignored and can be compared against the committed copies in [`results/`](01-identify-cyclospora-specific-kmers/results/).

## Step 5. Screen reads with the bait set

Steps 1 to 4 may be skipped if only the baits are needed: both FASTA files are committed under [`01-identify-cyclospora-specific-kmers/baits/`](01-identify-cyclospora-specific-kmers/baits/), and [`baits/README.md`](01-identify-cyclospora-specific-kmers/baits/README.md) documents the record format and the counting rule.

```bash
deacon index build -k 31 -w 1 -e 0 \
  baits/cyclospora_cayetanensis_rrna_core_nt_validated_baits.fasta \
  -o cyclospora_cayetanensis_rrna_core_nt_validated_k31w1.idx

deacon filter -m -a 20 -r 0 \
  cyclospora_cayetanensis_rrna_core_nt_validated_k31w1.idx \
  reads_R1.fastq.gz reads_R2.fastq.gz
```

`w=1` is required: it makes every 31-mer its own minimizer, so `-a` corresponds directly to the number of distinct diagnostic 31-mers found. In paired mode Deacon pools hits across mates, so a pair passes `-a 20` on 13 and 10 disjoint hits between the two reads. Each retained read is therefore recounted against the bait FASTA individually, which is why our reported counts are lower than Deacon's retained-read counts.

## Step 6. Reproduce the wastewater screen

Stage 02 applies the bait set to the libraries and calibrates the detection threshold. Its README, [`02-screen-wastewater-metagenomes/README.md`](02-screen-wastewater-metagenomes/README.md), is the entry point.

Two parts reproduce without the sequencing data. Every read counted in the analysis is committed under [`02-screen-wastewater-metagenomes/results/reads/`](02-screen-wastewater-metagenomes/results/reads/), and `python scripts/verify_published_reads.py` recounts all 1,156 of them against the bait set to confirm none falls below the threshold of 20. The heatmap in Figure 1, the Vega-Lite specification behind its interactive version, and the site-by-fortnight matrix regenerate with `python scripts/plot_heatmap.py` from [`sra_sample_summary.tsv`](02-screen-wastewater-metagenomes/results/sra_sample_summary.tsv). The threshold sweep behind Table 2 regenerates with `pixi run sweep` (or `python scripts/sweep_threshold.py`), which requires only Python 3, with no database or network access; it reads the per-read classifications in [`results/calibration/`](02-screen-wastewater-metagenomes/results/calibration/), derived from the public SRA `-a 1` screen — 205 runs, 394.5 billion reads — with every candidate read classified against `core_nt` along its full length.

Screening the libraries themselves requires the reads. Each pair is filtered with `deacon filter -m -a 20 -r 0` against the 1,184-bait index, and every retained read is then recounted against the bait FASTA individually, because Deacon pools k-mer hits across mates while the reported counts are per-read. More than 2,200 publicly available wastewater metagenomics datasets are available as of August 2026 in [BioProject PRJNA1247874](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1247874).
