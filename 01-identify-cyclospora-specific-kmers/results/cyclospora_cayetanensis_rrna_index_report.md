# *Cyclospora cayetanensis* mature-rRNA bait index verification

## Conclusion

PASS: A species-specific mature-rRNA k=31 bait set was produced.

## Provenance

| artifact | file | bytes | SHA-256 |
|---|---|---:|---|
| curated target FASTA | target_rrna.fasta | 8024 | 4ba5b870b08744cf429c91992d3af2a3cc6d9c8f3d0cb337696eeea901cd55c1 |
| final manifest | cyclospora_cayetanensis_rrna_specific.kmers.tsv | 1114544 | 689f7b5995473261ea2579ac3adbba0e0d30cb50fe19a71ccef29b3016777213 |
| final bait FASTA | cyclospora_cayetanensis_rrna_baits.fasta | 79093 | 21106d46dbecafb3cd784050b969b3a664cd21d454363583c1b1ff9f200c3200 |
| Meryl exact difference | exact_specific.kmers.tsv | 56780 | 1dcab0af3f5bc505ae70970568524ab4c4697e9b88b87a90ea9d7087ce3eefad |
| Meryl target-SILVA intersection | target_silva.kmers.tsv | 118082 | eb54f0d4b9523398687d7f460450ad6cb1e0f0db83d895568499d5c6d22a07c9 |
| Meryl target-Rfam intersection | target_rfam.kmers.tsv | 9282 | aa8fb97fbf56b4a4d111734a82e4dd41c92a75a6055e778ad0a899dcef08e1f3 |
| Meryl target-other-Cyclospora intersection | target_other_cyclospora.kmers.tsv | 62764 | 0e7dfcad69692b1925df439196e01f65ceb545cda28e810a575bf4eb7edad7b0 |
| Meryl genus-compatible difference | genus_compatible_pre_entropy.kmers.tsv | 62526 | e84565baf3e5c5cfd5c776bb8e690c7835f2147e47a50bcfd0b21037d3cd4be6 |
| source provenance manifest | input_manifest.tsv | 2666 | 84551b61ecd489ecc7f47a6a94091a3f14013d2f6856d83758a112f93452cd3b |
| locus metadata | target_loci.tsv | 7275 | a96785de56d7ec0bf6dc357d0e22bdbbcdb81d108fd05a269f69c161d80388ec |

### Pinned input sources

| source | release | URL | retrieved UTC | bytes | SHA-256 |
|---|---|---|---|---:|---|
| ncbi_eutils_other_cyclospora_rrna | reviewed-versioned-accessions-2026-07-23 | ncbi-eutils:efetch-versioned:config/other_cyclospora_accessions.txt | 2026-07-23T23:29:54Z | 62811 | 2fbc89b02a27523fc236b474ed91b9dd08410ab94cb17c429681345f63a68891 |
| ncbi_eutils_target_queries | 2026-07-23 | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&ids=AF111183.1,MPGL01000046.1,XR_003297357.1,XR_003297348.1,XR_003297351.1,XR_003297352.1,XR_003297353.1,XR_003297354.1,XR_003297364.1 | 2026-07-23T23:16:57Z | 7122 | 47e3ced57f0977bc4a1043b03e71d237414d626ec28d7f134239424435682874 |
| other_cyclospora_accession_snapshot | reviewed-2026-07-23 | tracked:config/other_cyclospora_accessions.txt | 2026-07-23T23:29:54Z | 1032 | b2dde1698b7ffc78fc77aac6f076e27823bdab76bf2a7060706c1105ae1a51ad |
| refseq_asm76915v2 | GCF_000769155.1 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/769/155/GCF_000769155.1_ASM76915v2/GCF_000769155.1_ASM76915v2_genomic.fna.gz | 2026-07-24T00:44:58Z | 14184212 | 84cb25992e1a0c8fc7f247141e50204a0937b39cde9230ba67801841dd432199 |
| refseq_ccayref3 | GCF_002999335.1 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/999/335/GCF_002999335.1_CcayRef3/GCF_002999335.1_CcayRef3_genomic.fna.gz | 2026-07-24T00:44:57Z | 14149982 | 3d37b65a32abd7879470dbf5688e3a5d3d1b194d8a62ce09ada3b1c8d3193b53 |
| rfam_5_8s | 15.1 | https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/fasta_files/RF00002.fa.gz | 2026-07-23T23:12:40Z | 794307 | bf9dc7aac7cc1b52a7ff117c6174ec2f13bd6a09c03a65c57915d03360257e71 |
| rfam_5s | 15.1 | https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/fasta_files/RF00001.fa.gz | 2026-07-23T23:12:39Z | 25482430 | c51060970ca8182891881787890eb706cda4c0ea499c04f874af488451e2ac0d |
| silva_lsu | 138.2 | https://ftp.arb-silva.de/release_138.2/Exports/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz | 2026-07-23T23:12:19Z | 70566330 | 25abefa760384984874f4f27afce25f44fbb1aa9a0779b3a5457e064e3422248 |
| silva_ssu | 138.2 | https://ftp.arb-silva.de/release_138.2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz | 2026-07-23T23:12:15Z | 201107829 | c779f51f1c605377f23d240005cbdd96068a77fdf85117f15d9f598b138a2072 |
| target_rrna_accession_snapshot | reviewed-2026-07-23 | tracked:config/target_rrna_accessions.txt | 2026-07-23T23:15:44Z | 3711 | f038c855ec1c77eeda405228cefb21d4729423166d20f001ee94d70141d62ae6 |

## Stage attrition

| rRNA class | target candidates | shared with SILVA | shared with Rfam | shared with other *Cyclospora* | exact species-specific | low-complexity rejected | passing baits |
|---|---:|---:|---:|---:|---:|---:|---:|
| 18S | 1824 | 1566 | 24 | 1720 | 89 | 17 | 72 |
| 5.8S | 126 | 0 | 126 | 126 | 0 | 0 | 0 |
| 28S | 3457 | 1907 | 0 | 0 | 1550 | 257 | 1293 |
| 5S | 154 | 0 | 123 | 0 | 31 | 8 | 23 |
| **all** | **5561** | **3473** | **273** | **1846** | **1670** | **282** | **1388** |

## Per-locus coverage

Coverage is descriptive. A covered 150 bp read start is a valid full-length start whose read contains an entire bait occurrence.

| target locus | class | bases | candidates | passing baits | bases covered | 150 bp read starts covered |
|---|---|---:|---:|---:|---:|---:|
| 18S\|AF111183.1\|NW_019209939.1\|1-1780\|+ | 18S | 1780 | 1750 | 52 | 152 | 285/1631 |
| 18S\|AF111183.1\|NW_020312507.1\|24-1815\|+ | 18S | 1792 | 1762 | 53 | 157 | 373/1643 |
| 28S\|MPGL01000046.1\|NW_019209236.1\|46-3532\|- | 28S | 3487 | 3457 | 1293 | 2169 | 2778/3338 |
| 5.8S\|XR_003297357.1\|NW_019209236.1\|4247-4402\|- | 5.8S | 156 | 126 | 0 | 0 | 0/7 |
| 5S\|XR_003297348.1\|NW_019209216.1\|779-900\|+ | 5S | 122 | 92 | 0 | 0 | 0/0 |
| 5S\|XR_003297348.1\|NW_019210658.1\|27886-28007\|+ | 5S | 122 | 92 | 0 | 0 | 0/0 |
| 5S\|XR_003297348.1\|NW_020312400.1\|6413-6534\|+ | 5S | 122 | 92 | 23 | 61 | 0/0 |

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
| bait records | 1388 |
| distinct indexed minimizers | 1388 |

## Scope limitation

Wastewater performance has not been evaluated. Static correctness does not establish wastewater sensitivity, specificity, or a validated sample-calling threshold.
