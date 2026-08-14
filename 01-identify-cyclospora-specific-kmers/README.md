# Stage 01. Identify *Cyclospora cayetanensis*-specific rRNA k-mers

This stage produced the historical bait set used to screen wastewater libraries in [stage 02](../02-screen-wastewater-metagenomes/): 1,184 canonical 31-mers that occur in *C. cayetanensis* mature ribosomal RNA, occur in no other ribosomal RNA in SILVA or Rfam, occur in no non-*cayetanensis Cyclospora* reference, and have no exact full-length match anywhere in NCBI `core_nt` assigned to another taxon. The corrected pre-screening set contains 1,670 baits; 282 restored baits still require taxonomic exact-match screening.

The narrative version, with the reasoning behind each threshold, is in the [main README](../REPRODUCING.md). This file is the map of what is here.

## Running it

```bash
pixi install
pixi run test                   # 90 unit tests and 3 fixture builds
pixi run build                  # design from pinned references -> 1,670 baits
pixi run validate /path/to/core_nt  # rebuild, search, and screen the current baits
```

## Layout

| Path | What it is |
|---|---|
| [`config/sources.tsv`](config/sources.tsv) | Every input pinned by URL, release, byte count, and SHA-256. The build refuses to run on an input whose identity has changed. |
| [`config/target_queries.tsv`](config/target_queries.tsv) | The nine reference transcripts used to find the mature loci, with the ≥98% identity and ≥90% query coverage thresholds and the 18S and 28S trim coordinates. |
| [`config/target_rrna_accessions.txt`](config/target_rrna_accessions.txt) | *C. cayetanensis* rRNA accessions excluded from the Rfam background so the target cannot subtract itself. |
| [`config/other_cyclospora_accessions.txt`](config/other_cyclospora_accessions.txt) | The 94 records from other *Cyclospora* species, subtracted so the baits cannot match a close relative. |
| [`scripts/build_rrna_bait.sh`](scripts/build_rrna_bait.sh) | The whole design: fetch, verify, locate loci, count 31-mers, subtract, filter low complexity, build and round-trip the index. |
| [`scripts/validate_core_nt.sh`](scripts/validate_core_nt.sh) | Rebuilds the current bait set, searches it against an operator-supplied `core_nt` database, and applies the screening decisions. |
| [`scripts/finalize_core_nt_validation.sh`](scripts/finalize_core_nt_validation.sh) | Applies core-nt BLAST evidence to the baits, rebuilds the index, and refuses to publish unless the background round trips are clean. |
| [`scripts/build_core_nt_validated_index.sh`](scripts/build_core_nt_validated_index.sh) | Thin wrapper that calls the finalizer with this repository's paths. |
| [`src/rrna_bait/`](src/rrna_bait/) | Python modules behind those scripts: `sources` (fetch and verify), `targets` (locus finding), `manifest` (k-mer bookkeeping), `core_nt` (validation decisions), `verify` (invariant checks), `core` (canonical k-mer helpers). |
| [`tests/`](tests/) | Unit tests plus end-to-end fixture builds that exercise the pipeline offline. |
| [`curated/`](curated/) | Where the located loci land. [`target_loci.tsv`](curated/target_loci.tsv) is one row per alignment site, 83 of them across 18 distinct genomic intervals, with assembly, coordinates, and strand. [`target_rrna.fasta`](curated/target_rrna.fasta) is the 7 distinct sequences those collapse to (2 × 18S, 1 × 28S, 1 × 5.8S, 3 × 5S), deduplicated by class and sequence in [`src/rrna_bait/targets.py`](src/rrna_bait/targets.py); column 1 of the TSV names the record each row fed. |
| [`results/`](results/) | The committed outputs of our run. See below. |
| [`baits/`](baits/) | The published bait FASTAs and how to use them: [`baits/README.md`](baits/README.md). |

## What is in `results/`

These are the historical manuscript outputs from the original 1,388-bait branch.

| File | Rows | What it records |
|---|---:|---|
| [`cyclospora_cayetanensis_rrna_kmer_manifest.tsv`](results/cyclospora_cayetanensis_rrna_kmer_manifest.tsv) | 5,561 | Every candidate 31-mer: source locus and offset, copy count, which backgrounds it matched, why it was kept or dropped, and its core-nt decision. The one file to read if you want to know what happened to a particular k-mer. |
| [`cyclospora_cayetanensis_rrna_index_summary.tsv`](results/cyclospora_cayetanensis_rrna_index_summary.tsv) | 92 keys | Attrition at every stage, overall and per rRNA class. The per-class survival rates quoted in the main README come from here. Note that `attrition.silva_shared_count`, `attrition.rfam_shared_count`, and `specificity.other_cyclospora_shared_count` are absolute memberships, not sequential removals: they overlap, so subtracting them in series double-counts. Table 1 of the main README reports the incremental removals, counted from the manifest. |
| [`cyclospora_cayetanensis_rrna_index_report.md`](results/cyclospora_cayetanensis_rrna_index_report.md) | — | The same verification in prose, with a SHA-256 for each intermediate artifact. |
| [`core_nt_bait_exact_match_blast.tsv`](results/core_nt_bait_exact_match_blast.tsv) | 5,012 | Raw exact-match hits from the historical `core_nt` search, with subject accession, taxid, and title. These hits do not apply to the corrected 1,670-bait build. |
| [`cyclospora_cayetanensis_core_nt_validation.tsv`](results/cyclospora_cayetanensis_core_nt_validation.tsv) | 1,388 | Per-bait verdict: 1,184 `PASS_CORE_NT`, 204 `REJECT_CORE_NT_EXACT_NON_TARGET`, with target and non-target hit counts. |
| [`cyclospora_cayetanensis_core_nt_near_hits.tsv`](results/cyclospora_cayetanensis_core_nt_near_hits.tsv) | 0 | BLAST rows matching a bait over 30 of its 31 bases. Empty by construction: the search required exact, full-length matches, so no partial match could be reported. It would populate only if the search were repeated at a permissive identity setting. |
| [`input_manifest.tsv`](results/input_manifest.tsv) | 10 | What was actually downloaded, when, and with what checksum. |

## Regenerated versus committed

`pixi run build` writes its working files to `kmers/`, `reports/`, `provenance/`, `work/`, `logs/`, and the downloads to `silva/`, `pos/`, and `background/`. All of those are gitignored, so a rerun never masquerades as a change to the record. It does overwrite `curated/`, which is committed on purpose: if your rebuild is faithful, `git diff curated/` stays empty. The corrected attrition reports differ from the committed historical reports at the low-complexity step.
