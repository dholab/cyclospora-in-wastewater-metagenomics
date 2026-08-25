# Threshold calibration — why a read must carry 24 diagnostic 31-mers

This directory holds the evidence for the detection threshold. The corrected 1,464-bait index was applied permissively at `deacon filter -a 1 -r 0`, and every retained read was classified by comparing its best local-alignment bit scores within and outside the declared *Cyclospora* target scope. The threshold sweep then asks how many target, non-target, tied, and no-hit reads remain at each minimum bait count.

## Calibration population

The calibration contains 108,474 reads carrying at least one bait. Duplicate sequences share a representative query, leaving 10,894 distinct sequences for the `core_nt` search. Those exact queries are committed as [`calibration_read_blast_queries.fasta.gz`](calibration_read_blast_queries.fasta.gz).

| class | reads | rule |
|---|---:|---|
| `target` | 1,284 | the best target bit score exceeds the best non-target score |
| `non_target` | 107,040 | the best non-target bit score exceeds the best target score |
| `top_tie` | 123 | the best target and non-target scores tie within 0.1 bits |
| `no_hit` | 27 | no usable `core_nt` alignment |

The target scope comprises 26 declared *Cyclospora* taxids, including *C. cayetanensis* taxid 88456. A subject whose usable taxids fall entirely outside that scope, or which has no usable taxid, is non-target evidence. Ties remain separate because the alignment evidence does not decide between target and non-target.

## Threshold sweep

[`threshold_read_counts.tsv`](threshold_read_counts.tsv) records all thresholds. Selected values are:

| threshold | target | non-target | tie | no hit |
|---:|---:|---:|---:|---:|
| 1 | 1,284 | 107,040 | 123 | 27 |
| 5 | 1,121 | 246 | 38 | 0 |
| 10 | 1,071 | 18 | 33 | 0 |
| 15 | 988 | 17 | 32 | 0 |
| 20 | 941 | 0 | 32 | 0 |
| 23 | 855 | 0 | 16 | 0 |
| **24** | **853** | **0** | **0** | **0** |
| 25 | 809 | 0 | 0 | 0 |

The highest bait count on a non-target read is 15. Tied reads remain through a bait count of 23. Therefore **24 is the lowest threshold at which every retained calibration read is target-classified**.

## Files

| File | What it holds |
|---|---|
| [`read_blast_hits.tsv.gz`](read_blast_hits.tsv.gz) | compressed raw `core_nt` BLAST alignments used to classify the distinct calibration sequences |
| [`read_blast_map.tsv`](read_blast_map.tsv) | one row per calibration read: source, unique read ID, sequence-derived representative ID, bait count, and read length |
| [`read_blast_deacon.tsv`](read_blast_deacon.tsv) | one row per calibration read: bait count and independently calculated `core_nt` class |
| [`calibration_read_blast_queries.fasta.gz`](calibration_read_blast_queries.fasta.gz) | the 10,894 distinct 150–151 nt sequences submitted to `core_nt` |
| [`threshold_read_counts.tsv`](threshold_read_counts.tsv) | reads surviving each threshold, by class, and runs with at least one surviving read |

The raw alignment archive was created deterministically with `gzip -9 -n`. Its SHA-256 is `2f97a75ab265bae5fc5de0eb33de10f607d4090fc18477ab59d680d12459f4d2`; the decompressed TSV SHA-256 is `d66285d7c2b4b647789470222b9cb28cbe1afeae8f64a5528c72e1d6efa4cc8c`. The search used BLAST 2.16.0 and BLASTDB version 5. The `core_nt` database reported 117,146,392 sequences, 894,401,043,626 bases, and a build date of 28 July 2025 at 23:24.

## Rebuilding the sweep

The threshold table regenerates from the committed per-read evidence without a database, cluster, or network:

```bash
pixi run sweep
```

Repeating the classifications requires the committed raw `core_nt` BLAST output, but not the database itself:

```bash
gzip -dc results/calibration/read_blast_hits.tsv.gz > read_blast_hits.tsv
python scripts/classify_reads.py \
  --blast read_blast_hits.tsv
pixi run sweep
```

Repeating the BLAST search requires the approximately 285 GB `core_nt` database and the distinct calibration sequences:

```bash
gzip -dc results/calibration/calibration_read_blast_queries.fasta.gz > query.fasta
blastn -task blastn -db core_nt -query query.fasta \
  -evalue 1e-10 -max_target_seqs 100 -dust no \
  -outfmt '6 qseqid qlen saccver staxids pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp stitle' \
  -out read_blast_hits.tsv
```
