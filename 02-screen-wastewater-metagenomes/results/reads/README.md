# Diagnostic reads, one FASTA per SRA run

Each file holds the reads from one public SRA run that carry at least 20
diagnostic *Cyclospora cayetanensis* 31-mers of their own, at the calibrated
`deacon filter -a 20 -r 0`. Files are named `<CODE>_<YYYYMMDD>__<accession>`,
where `<CODE>` is the sewershed's `public_code_CASPER_SRA` and `<accession>` is
the SRA run. Only the 73 runs with at least one diagnostic read appear here;
every screened run, positive or not, is in
[`../sra_sample_summary.tsv`](../sra_sample_summary.tsv).

`scripts/verify_published_reads.py` recounts every read below against the bait
set and confirms none falls under the threshold.

| Run | Sewershed | Collection date | Diagnostic reads | Distinct |
|---|---|---|--:|--:|
| [`SRR36790976`](CA_Ontario_20250903__SRR36790976.diagnostic_reads.fasta.gz) | `CA_Ontario` | 2025-09-03 | 2 | 2 |
| [`SRR36790971`](CA_Ontario_20250922__SRR36790971.diagnostic_reads.fasta.gz) | `CA_Ontario` | 2025-09-22 | 16 | 9 |
| [`SRR39695241`](CA_Ontario_20260413__SRR39695241.diagnostic_reads.fasta.gz) | `CA_Ontario` | 2026-04-13 | 62 | 24 |
| [`SRR39695379`](CA_PaloAlto_20260622__SRR39695379.diagnostic_reads.fasta.gz) | `CA_PaloAlto` | 2026-06-22 | 2 | 1 |
| [`SRR35931370`](CA_Riverside_20250529__SRR35931370.diagnostic_reads.fasta.gz) | `CA_Riverside` | 2025-05-29 | 1 | 1 |
| [`SRR35931375`](CA_Riverside_20250616__SRR35931375.diagnostic_reads.fasta.gz) | `CA_Riverside` | 2025-06-16 | 20 | 4 |
| [`SRR35931376`](CA_Riverside_20250619__SRR35931376.diagnostic_reads.fasta.gz) | `CA_Riverside` | 2025-06-19 | 53 | 15 |
| [`SRR35931346`](CA_Riverside_20250721__SRR35931346.diagnostic_reads.fasta.gz) | `CA_Riverside` | 2025-07-21 | 9 | 1 |
| [`SRR35931343`](CA_Riverside_20250804__SRR35931343.diagnostic_reads.fasta.gz) | `CA_Riverside` | 2025-08-04 | 10 | 6 |
| [`SRR35931364`](CA_Riverside_20250807__SRR35931364.diagnostic_reads.fasta.gz) | `CA_Riverside` | 2025-08-07 | 1 | 1 |
| [`SRR36814892`](CA_Sacramento_20251125__SRR36814892.diagnostic_reads.fasta.gz) | `CA_Sacramento` | 2025-11-25 | 2 | 2 |
| [`SRR35987577`](CHI-A_20250420__SRR35987577.diagnostic_reads.fasta.gz) | `CHI-A` | 2025-04-20 | 2 | 2 |
| [`SRR35987589`](CHI-A_20250622__SRR35987589.diagnostic_reads.fasta.gz) | `CHI-A` | 2025-06-22 | 8 | 2 |
| [`SRR35987646`](CHI-A_20250810__SRR35987646.diagnostic_reads.fasta.gz) | `CHI-A` | 2025-08-10 | 6 | 2 |
| [`SRR39695338`](CHI-A_20260524__SRR39695338.diagnostic_reads.fasta.gz) | `CHI-A` | 2026-05-24 | 8 | 3 |
| [`SRR39695329`](CHI-A_20260628__SRR39695329.diagnostic_reads.fasta.gz) | `CHI-A` | 2026-06-28 | 9 | 5 |
| [`SRR35987601`](CHI-B_20250810__SRR35987601.diagnostic_reads.fasta.gz) | `CHI-B` | 2025-08-10 | 26 | 6 |
| [`SRR39695057`](CHI-B_20260628__SRR39695057.diagnostic_reads.fasta.gz) | `CHI-B` | 2026-06-28 | 24 | 5 |
| [`SRR35987579`](CHI-C_20250602__SRR35987579.diagnostic_reads.fasta.gz) | `CHI-C` | 2025-06-02 | 21 | 9 |
| [`SRR35987635`](CHI-C_20250728__SRR35987635.diagnostic_reads.fasta.gz) | `CHI-C` | 2025-07-28 | 14 | 2 |
| [`SRR39695110`](CHI-C_20260615__SRR39695110.diagnostic_reads.fasta.gz) | `CHI-C` | 2026-06-15 | 9 | 6 |
| [`SRR35987603`](CHI-D2_20250622__SRR35987603.diagnostic_reads.fasta.gz) | `CHI-D2` | 2025-06-22 | 23 | 3 |
| [`SRR36995842`](FL_Miami_20250513__SRR36995842.diagnostic_reads.fasta.gz) | `FL_Miami` | 2025-05-13 | 1 | 1 |
| [`SRR39695226`](FL_Miami_20260526__SRR39695226.diagnostic_reads.fasta.gz) | `FL_Miami` | 2026-05-26 | 2 | 2 |
| [`SRR39695224`](FL_Miami_20260608__SRR39695224.diagnostic_reads.fasta.gz) | `FL_Miami` | 2026-06-08 | 9 | 4 |
| [`SRR39695220`](FL_Miami_20260629__SRR39695220.diagnostic_reads.fasta.gz) | `FL_Miami` | 2026-06-29 | 6 | 4 |
| [`SRR38294946`](IA_Ottumwa_20260121__SRR38294946.diagnostic_reads.fasta.gz) | `IA_Ottumwa` | 2026-01-21 | 9 | 2 |
| [`SRR36876975`](ID_Boise_20250630__SRR36876975.diagnostic_reads.fasta.gz) | `ID_Boise` | 2025-06-30 | 9 | 1 |
| [`SRR36876973`](ID_Boise_20250707__SRR36876973.diagnostic_reads.fasta.gz) | `ID_Boise` | 2025-07-07 | 14 | 3 |
| [`SRR39695079`](ID_Boise_20260629__SRR39695079.diagnostic_reads.fasta.gz) | `ID_Boise` | 2026-06-29 | 10 | 5 |
| [`SRR37076038`](MA_Boston_DITPN_20250528__SRR37076038.diagnostic_reads.fasta.gz) | `MA_Boston_DITPN` | 2025-05-28 | 1 | 1 |
| [`SRR37076056`](MA_Boston_DITPN_20250702__SRR37076056.diagnostic_reads.fasta.gz) | `MA_Boston_DITPN` | 2025-07-02 | 10 | 2 |
| [`SRR37059701`](MA_Boston_DITPN_20250709__SRR37059701.diagnostic_reads.fasta.gz) | `MA_Boston_DITPN` | 2025-07-09 | 4 | 3 |
| [`SRR37367552`](MA_Boston_DITPN_20250709__SRR37367552.diagnostic_reads.fasta.gz) | `MA_Boston_DITPN` | 2025-07-09 | 4 | 3 |
| [`SRR37367644`](MA_Boston_DITPN_20250716__SRR37367644.diagnostic_reads.fasta.gz) | `MA_Boston_DITPN` | 2025-07-16 | 20 | 4 |
| [`SRR37059724`](MA_Boston_DITPN_20250813__SRR37059724.diagnostic_reads.fasta.gz) | `MA_Boston_DITPN` | 2025-08-13 | 2 | 2 |
| [`SRR37367633`](MA_Boston_DITPN_20250813__SRR37367633.diagnostic_reads.fasta.gz) | `MA_Boston_DITPN` | 2025-08-13 | 6 | 3 |
| [`SRR37367631`](MA_Boston_DITPN_20250820__SRR37367631.diagnostic_reads.fasta.gz) | `MA_Boston_DITPN` | 2025-08-20 | 8 | 4 |
| [`SRR37367603`](MA_Boston_DITPN_20250903__SRR37367603.diagnostic_reads.fasta.gz) | `MA_Boston_DITPN` | 2025-09-03 | 2 | 2 |
| [`SRR39695202`](MA_Boston_DITPN_20260422__SRR39695202.diagnostic_reads.fasta.gz) | `MA_Boston_DITPN` | 2026-04-22 | 12 | 8 |
| [`SRR39695135`](MA_Boston_DITPN_20260617__SRR39695135.diagnostic_reads.fasta.gz) | `MA_Boston_DITPN` | 2026-06-17 | 16 | 2 |
| [`SRR37076043`](MA_Boston_DITPS_20250611__SRR37076043.diagnostic_reads.fasta.gz) | `MA_Boston_DITPS` | 2025-06-11 | 10 | 3 |
| [`SRR37076041`](MA_Boston_DITPS_20250625__SRR37076041.diagnostic_reads.fasta.gz) | `MA_Boston_DITPS` | 2025-06-25 | 12 | 2 |
| [`SRR37367654`](MA_Boston_DITPS_20250910__SRR37367654.diagnostic_reads.fasta.gz) | `MA_Boston_DITPS` | 2025-09-10 | 4 | 4 |
| [`SRR35904059`](MO_Columbia_20240723__SRR35904059.diagnostic_reads.fasta.gz) | `MO_Columbia` | 2024-07-23 | 1 | 1 |
| [`SRR35904060`](MO_Columbia_20240723__SRR35904060.diagnostic_reads.fasta.gz) | `MO_Columbia` | 2024-07-23 | 2 | 2 |
| [`SRR35904050`](MO_Columbia_20240917__SRR35904050.diagnostic_reads.fasta.gz) | `MO_Columbia` | 2024-09-17 | 2 | 2 |
| [`SRR37092442`](MO_Columbia_20250729__SRR37092442.diagnostic_reads.fasta.gz) | `MO_Columbia` | 2025-07-29 | 22 | 4 |
| [`SRR37092444`](MO_Columbia_20250812__SRR37092444.diagnostic_reads.fasta.gz) | `MO_Columbia` | 2025-08-12 | 8 | 2 |
| [`SRR39695205`](MO_Columbia_20260623__SRR39695205.diagnostic_reads.fasta.gz) | `MO_Columbia` | 2026-06-23 | 1 | 1 |
| [`SRR39695198`](MO_KC_BlueRiver_20260618__SRR39695198.diagnostic_reads.fasta.gz) | `MO_KC_BlueRiver` | 2026-06-18 | 3 | 3 |
| [`SRR39695196`](MO_KC_BlueRiver_20260628__SRR39695196.diagnostic_reads.fasta.gz) | `MO_KC_BlueRiver` | 2026-06-28 | 10 | 6 |
| [`SRR39695291`](MO_KC_Westside_20260618__SRR39695291.diagnostic_reads.fasta.gz) | `MO_KC_Westside` | 2026-06-18 | 5 | 1 |
| [`SRR35939760`](MO_Monett_20250728__SRR35939760.diagnostic_reads.fasta.gz) | `MO_Monett` | 2025-07-28 | 12 | 8 |
| [`SRR36876989`](MO_Monett_20251117__SRR36876989.diagnostic_reads.fasta.gz) | `MO_Monett` | 2025-11-17 | 16 | 9 |
| [`SRR39695014`](MO_Monett_20260622__SRR39695014.diagnostic_reads.fasta.gz) | `MO_Monett` | 2026-06-22 | 7 | 5 |
| [`SRR38294111`](MO_STL_Coldwater_20250707__SRR38294111.diagnostic_reads.fasta.gz) | `MO_STL_Coldwater` | 2025-07-07 | 20 | 8 |
| [`SRR38294083`](MO_STL_Coldwater_20250806__SRR38294083.diagnostic_reads.fasta.gz) | `MO_STL_Coldwater` | 2025-08-06 | 5 | 4 |
| [`SRR38294061`](MO_STL_Coldwater_20251027__SRR38294061.diagnostic_reads.fasta.gz) | `MO_STL_Coldwater` | 2025-10-27 | 1 | 1 |
| [`SRR38294146`](MO_STL_Lemay_20250806__SRR38294146.diagnostic_reads.fasta.gz) | `MO_STL_Lemay` | 2025-08-06 | 2 | 2 |
| [`SRR39695279`](MO_STL_Lemay_20260414__SRR39695279.diagnostic_reads.fasta.gz) | `MO_STL_Lemay` | 2026-04-14 | 8 | 2 |
| [`SRR39694993`](MO_STL_Lemay_20260630__SRR39694993.diagnostic_reads.fasta.gz) | `MO_STL_Lemay` | 2026-06-30 | 12 | 5 |
| [`SRR38294184`](MO_STL_MORiver_20250620__SRR38294184.diagnostic_reads.fasta.gz) | `MO_STL_MORiver` | 2025-06-20 | 5 | 1 |
| [`SRR38294180`](MO_STL_MORiver_20250710__SRR38294180.diagnostic_reads.fasta.gz) | `MO_STL_MORiver` | 2025-07-10 | 6 | 1 |
| [`SRR38294178`](MO_STL_MORiver_20250714__SRR38294178.diagnostic_reads.fasta.gz) | `MO_STL_MORiver` | 2025-07-14 | 4 | 3 |
| [`SRR38294173`](MO_STL_MORiver_20250811__SRR38294173.diagnostic_reads.fasta.gz) | `MO_STL_MORiver` | 2025-08-11 | 11 | 5 |
| [`SRR38294165`](MO_STL_MORiver_20251103__SRR38294165.diagnostic_reads.fasta.gz) | `MO_STL_MORiver` | 2025-11-03 | 4 | 2 |
| [`SRR39695272`](MO_STL_MORiver_20260617__SRR39695272.diagnostic_reads.fasta.gz) | `MO_STL_MORiver` | 2026-06-17 | 4 | 4 |
| [`SRR37006679`](NYC-Hospital-B_20250612__SRR37006679.diagnostic_reads.fasta.gz) | `NYC-Hospital-B` | 2025-06-12 | 447 | 76 |
| [`SRR39695152`](OK-A_20260602__SRR39695152.diagnostic_reads.fasta.gz) | `OK-A` | 2026-06-02 | 10 | 7 |
| [`SRR36861556`](OK-B_20250826__SRR36861556.diagnostic_reads.fasta.gz) | `OK-B` | 2025-08-26 | 23 | 5 |
| [`SRR36876823`](Southern_California_20250706__SRR36876823.diagnostic_reads.fasta.gz) | `Southern_California` | 2025-07-06 | 2 | 1 |
| [`SRR36876846`](Southern_California_20250806__SRR36876846.diagnostic_reads.fasta.gz) | `Southern_California` | 2025-08-06 | 4 | 4 |
