# *Cyclospora cayetanensis* mature-rRNA bait index verification

## Conclusion

PASS: A species-specific mature-rRNA k=31 bait set was produced.

## Provenance

| artifact | file | bytes | SHA-256 |
|---|---|---:|---|
| curated target FASTA | target_rrna.fasta | 8024 | 4ba5b870b08744cf429c91992d3af2a3cc6d9c8f3d0cb337696eeea901cd55c1 |
| final manifest | cyclospora_cayetanensis_rrna_specific.kmers.tsv | 1106930 | 85641010a126260804dfa0d2639a42d880233ad04f50f0bdf77a4e646a182c82 |
| final bait FASTA | cyclospora_cayetanensis_rrna_baits.fasta | 95159 | 7f219cb0b4107c69facaf40049ec914e667dc59c3d03150be5a6de78d7a17213 |
| Meryl exact difference | exact_specific.kmers.tsv | 56780 | a63357dd70dad282716cac588b03355dc12904fabab0b87c2b81fb9dc12dd895 |
| Meryl target-SILVA intersection | target_silva.kmers.tsv | 118082 | d253f98b041cd10d6f9d60600767938624d57c8298d6a4b1ba3aec1bb2b64045 |
| Meryl target-Rfam intersection | target_rfam.kmers.tsv | 9282 | 380e8563de19fcbaa9db98869cb2e04844cfb5041738864848606e2e45be96ba |
| Meryl target-other-Cyclospora intersection | target_other_cyclospora.kmers.tsv | 62764 | 6c2e1d57d92189a2643b82446ddcde53ea34cd52617c76fc6730e9b714724419 |
| Meryl genus-compatible difference | genus_compatible_pre_entropy.kmers.tsv | 62526 | 9d26753c7b5b3018a06d0a9adc77a762ebc80a745592640d4723721c6fcba84e |
| source provenance manifest | input_manifest.tsv | 2215 | d149dbaa3d9ca625ee6c07830ae3a910fd05687b0224bec01f1fefe6209d7953 |
| locus metadata | target_loci.tsv | 7275 | a96785de56d7ec0bf6dc357d0e22bdbbcdb81d108fd05a269f69c161d80388ec |

### Pinned input sources

| source | release | URL | retrieved UTC | bytes | SHA-256 |
|---|---|---|---|---:|---|
| ncbi_eutils_other_cyclospora_rrna | reviewed-versioned-accessions-2026-07-23 | ncbi-eutils:efetch-versioned:config/other_cyclospora_accessions.txt | 2026-08-25T12:21:41Z | 62811 | 2fbc89b02a27523fc236b474ed91b9dd08410ab94cb17c429681345f63a68891 |
| ncbi_eutils_target_queries | 2026-07-23 | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&ids=AF111183.1,MPGL01000046.1,XR_003297357.1,XR_003297348.1,XR_003297351.1,XR_003297352.1,XR_003297353.1,XR_003297354.1,XR_003297364.1 | 2026-08-25T12:15:07Z | 7122 | 47e3ced57f0977bc4a1043b03e71d237414d626ec28d7f134239424435682874 |
| refseq_asm76915v2 | GCF_000769155.1 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/769/155/GCF_000769155.1_ASM76915v2/GCF_000769155.1_ASM76915v2_genomic.fna.gz | 2026-08-25T12:15:06Z | 14184212 | 84cb25992e1a0c8fc7f247141e50204a0937b39cde9230ba67801841dd432199 |
| refseq_ccayref3 | GCF_002999335.1 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/999/335/GCF_002999335.1_CcayRef3/GCF_002999335.1_CcayRef3_genomic.fna.gz | 2026-08-25T12:15:02Z | 14149982 | 3d37b65a32abd7879470dbf5688e3a5d3d1b194d8a62ce09ada3b1c8d3193b53 |
| rfam_5_8s | 15.1 | https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/fasta_files/RF00002.fa.gz | 2026-08-25T12:14:57Z | 794307 | bf9dc7aac7cc1b52a7ff117c6174ec2f13bd6a09c03a65c57915d03360257e71 |
| rfam_5s | 15.1 | https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/fasta_files/RF00001.fa.gz | 2026-08-25T12:14:56Z | 25482430 | c51060970ca8182891881787890eb706cda4c0ea499c04f874af488451e2ac0d |
| silva_lsu | 138.2 | https://ftp.arb-silva.de/release_138.2/Exports/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz | 2026-08-25T12:14:34Z | 70566330 | 25abefa760384984874f4f27afce25f44fbb1aa9a0779b3a5457e064e3422248 |
| silva_ssu | 138.2 | https://ftp.arb-silva.de/release_138.2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz | 2026-08-25T12:14:17Z | 201107829 | c779f51f1c605377f23d240005cbdd96068a77fdf85117f15d9f598b138a2072 |

## Stage attrition

| rRNA class | target candidates | shared with SILVA | shared with Rfam | shared with other *Cyclospora* | exact species-specific | low-complexity rejected | passing baits |
|---|---:|---:|---:|---:|---:|---:|---:|
| 18S | 1824 | 1566 | 24 | 1720 | 89 | 0 | 89 |
| 5.8S | 126 | 0 | 126 | 126 | 0 | 0 | 0 |
| 28S | 3457 | 1907 | 0 | 0 | 1550 | 0 | 1550 |
| 5S | 154 | 0 | 123 | 0 | 31 | 0 | 31 |
| **all** | **5561** | **3473** | **273** | **1846** | **1670** | **0** | **1670** |

## Per-locus coverage

Coverage is descriptive. A covered 150 bp read start is a valid full-length start whose read contains an entire bait occurrence.

| target locus | class | bases | candidates | passing baits | bases covered | 150 bp read starts covered |
|---|---|---:|---:|---:|---:|---:|
| 18S\|AF111183.1\|NW_019209939.1\|1-1780\|+ | 18S | 1780 | 1750 | 63 | 153 | 286/1631 |
| 18S\|AF111183.1\|NW_020312507.1\|24-1815\|+ | 18S | 1792 | 1762 | 67 | 157 | 373/1643 |
| 28S\|MPGL01000046.1\|NW_019209236.1\|46-3532\|- | 28S | 3487 | 3457 | 1550 | 2173 | 2781/3338 |
| 5.8S\|XR_003297357.1\|NW_019209236.1\|4247-4402\|- | 5.8S | 156 | 126 | 0 | 0 | 0/7 |
| 5S\|XR_003297348.1\|NW_019209216.1\|779-900\|+ | 5S | 122 | 92 | 0 | 0 | 0/0 |
| 5S\|XR_003297348.1\|NW_019210658.1\|27886-28007\|+ | 5S | 122 | 92 | 0 | 0 | 0/0 |
| 5S\|XR_003297348.1\|NW_020312400.1\|6413-6534\|+ | 5S | 122 | 92 | 31 | 61 | 0/0 |

## Other-*Cyclospora* sharing

1846 target candidate(s) were shared with the available non-*cayetanensis Cyclospora* references and were excluded from the species-specific bait set.

Subtraction is intentionally conservative: complete response records are used, so any ITS or other non-rRNA sequence carried in a mixed record is also subtracted.

| rRNA class | shared candidates |
|---|---:|
| 18S | 1720 |
| 5.8S | 126 |
| 28S | 0 |
| 5S | 0 |

## Genus-level feasibility

Status: `SPECIES_SPECIFIC_CANDIDATES`.

1839 target candidate(s) remain after SILVA and Rfam subtraction before other-*Cyclospora* and entropy filtering.

## Final index properties

| property | value |
|---|---:|
| index emitted | yes |
| k | 31 |
| w | 1 |
| bait records | 1670 |
| distinct indexed minimizers | 1670 |

## Scope limitation

Wastewater performance has not been evaluated. Static correctness does not establish wastewater sensitivity, specificity, or a validated sample-calling threshold.
