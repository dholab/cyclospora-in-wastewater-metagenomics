import csv
from pathlib import Path

from rrna_bait.core import canonical_kmer, reverse_complement
from rrna_bait.manifest import (
    finalize_manifest,
    main,
    read_meryl_print,
    target_kmer_metadata,
    write_raw_manifest,
)


PASS_KMER = "ACGTTGCAACGTTGCAACGTTGCAACGTTGC"


def _write_loci(path: Path, rows: list[tuple[str, str, str, int, int, str]]) -> None:
    path.write_text(
        "target_record_id\trrna_class\tquery_accession\t"
        "subject_accession\tstart\tend\tstrand\n"
        + "".join(
            f"{record_id}\t{rrna_class}\tquery\t{subject}\t{start}\t{end}\t{strand}\n"
            for record_id, rrna_class, subject, start, end, strand in rows
        )
    )


def _write_meryl_print(path: Path, *kmers: str) -> None:
    path.write_text("".join(f"{kmer}\t1\n" for kmer in kmers))


def test_target_metadata_canonicalizes_strands_and_retains_all_loci(tmp_path: Path):
    reverse = reverse_complement(PASS_KMER)
    target_fasta = tmp_path / "targets.fasta"
    target_fasta.write_text(
        f">18S|representative_a\n{PASS_KMER}\n"
        f">28S|representative_b\n{reverse}\n"
    )
    loci_tsv = tmp_path / "loci.tsv"
    _write_loci(
        loci_tsv,
        [
            ("18S|representative_a", "18S", "REF_A", 101, 131, "+"),
            ("18S|representative_a", "18S", "REF_B", 201, 231, "-"),
            ("28S|representative_b", "28S", "REF_C", 301, 331, "+"),
        ],
    )

    metadata = target_kmer_metadata(target_fasta, loci_tsv)

    assert list(metadata) == [PASS_KMER]
    assert metadata[PASS_KMER].rrna_classes == {"18S", "28S"}
    assert metadata[PASS_KMER].target_records == {
        "18S|representative_a",
        "28S|representative_b",
    }
    assert metadata[PASS_KMER].zero_based_starts == {0}
    assert metadata[PASS_KMER].occurrences == {
        ("18S|representative_a", 0),
        ("28S|representative_b", 0),
    }
    assert metadata[PASS_KMER].copy_count == 2


def test_target_metadata_falls_back_to_unmapped_fasta_record(tmp_path: Path):
    target_fasta = tmp_path / "targets.fasta"
    target_fasta.write_text(f">18S|synthetic\n{PASS_KMER}\n")
    loci_tsv = tmp_path / "loci.tsv"
    _write_loci(loci_tsv, [])

    metadata = target_kmer_metadata(target_fasta, loci_tsv)

    assert metadata[PASS_KMER].rrna_classes == {"18S"}
    assert metadata[PASS_KMER].target_records == {"18S|synthetic"}
    assert metadata[PASS_KMER].copy_count == 1


def test_read_meryl_print_uses_first_field_as_authoritative_set(tmp_path: Path):
    path = tmp_path / "printed.txt"
    path.write_text(
        f"{PASS_KMER}\t4\n\n{'A' * 31} 2\n"
        f"{reverse_complement(PASS_KMER)}\t9\n"
    )

    assert read_meryl_print(path) == {PASS_KMER, "A" * 31}


def test_raw_manifest_distinguishes_each_static_background(tmp_path: Path):
    sequences = [
        "A" * 31,
        "AC" * 15 + "A",
        "AG" * 15 + "A",
    ]
    target_fasta = tmp_path / "targets.fasta"
    target_fasta.write_text(
        "".join(f">18S|record_{index}\n{sequence}\n" for index, sequence in enumerate(sequences))
    )
    loci_tsv = tmp_path / "loci.tsv"
    _write_loci(loci_tsv, [])
    canonical = [canonical_kmer(sequence) for sequence in sequences]
    silva = tmp_path / "silva.txt"
    rfam = tmp_path / "rfam.txt"
    related = tmp_path / "related.txt"
    exact = tmp_path / "exact.txt"
    _write_meryl_print(silva, canonical[1])
    _write_meryl_print(rfam, canonical[2])
    _write_meryl_print(related, canonical[0])
    _write_meryl_print(exact)
    raw_manifest = tmp_path / "raw.tsv"

    count = write_raw_manifest(
        target_fasta,
        loci_tsv,
        silva,
        rfam,
        related,
        exact,
        raw_manifest,
        tmp_path / "raw.fasta",
    )

    rows = {
        row["kmer"]: row
        for row in csv.DictReader(raw_manifest.open(), delimiter="\t")
    }
    assert count == 0
    assert rows[canonical[0]]["rejection_reason"] == "other_cyclospora"
    assert rows[canonical[1]]["rejection_reason"] == "silva"
    assert rows[canonical[2]]["rejection_reason"] == "rfam"
    assert rows[canonical[0]]["in_other_cyclospora"] == "1"
    assert rows[canonical[0]]["in_silva"] == "0"
    assert rows[canonical[0]]["in_rfam"] == "0"


def test_finalize_uses_entropy_dump_and_writes_stable_sorted_baits(tmp_path: Path):
    static_kmer = "A" * 31
    low_complexity_kmer = "AC" * 15 + "A"
    target_fasta = tmp_path / "targets.fasta"
    target_fasta.write_text(
        f">18S|static\n{static_kmer}\n"
        f">5S|low\n{low_complexity_kmer}\n"
        f">18S|pass\n{PASS_KMER}\n"
    )
    loci_tsv = tmp_path / "loci.tsv"
    _write_loci(loci_tsv, [])
    canonical_static = canonical_kmer(static_kmer)
    canonical_low = canonical_kmer(low_complexity_kmer)
    silva = tmp_path / "silva.txt"
    rfam = tmp_path / "rfam.txt"
    related = tmp_path / "related.txt"
    exact = tmp_path / "exact.txt"
    _write_meryl_print(silva, canonical_static)
    _write_meryl_print(rfam, canonical_static)
    _write_meryl_print(related, canonical_static)
    _write_meryl_print(exact, canonical_low, PASS_KMER)
    raw_manifest = tmp_path / "raw.tsv"
    raw_baits = tmp_path / "raw.fasta"

    assert (
        write_raw_manifest(
            target_fasta,
            loci_tsv,
            silva,
            rfam,
            related,
            exact,
            raw_manifest,
            raw_baits,
        )
        == 2
    )
    raw_rows = list(csv.DictReader(raw_manifest.open(), delimiter="\t"))
    assert raw_rows[0]["rejection_reason"] == "other_cyclospora;silva;rfam"
    entropy_pass = tmp_path / "entropy_pass.fasta"
    entropy_pass.write_text(
        f">dumped_by_deacon\n{reverse_complement(PASS_KMER)}\n"
    )
    final_manifest = tmp_path / "final.tsv"
    final_baits = tmp_path / "final.fasta"

    assert (
        finalize_manifest(
            raw_manifest,
            entropy_pass,
            final_manifest,
            final_baits,
        )
        == 1
    )

    rows = list(csv.DictReader(final_manifest.open(), delimiter="\t"))
    assert [row["status"] for row in rows] == [
        "REJECT_STATIC_BACKGROUND",
        "REJECT_LOW_COMPLEXITY",
        "PASS",
    ]
    assert rows[-1]["rejection_reason"] == "none"
    assert all(
        not line.endswith("\t")
        for line in final_manifest.read_text().splitlines()
    )
    assert final_baits.read_text() == (
        ">cc_rrna_kmer_000001|18S\n"
        "ACGTTGCAACGTTGCAACGTTGCAACGTTGC\n"
    )


def test_final_bait_ids_follow_lexicographic_kmer_order(tmp_path: Path):
    later_kmer = "C" * 31
    raw_manifest = tmp_path / "raw.tsv"
    raw_manifest.write_text(
        "kmer\trrna_classes\ttarget_records\ttarget_starts\t"
        "target_copy_count\tin_silva\tin_rfam\tin_other_cyclospora\t"
        "in_exact_difference\tstatus\trejection_reason\n"
        f"{later_kmer}\t28S\tlater\t0\t1\t0\t0\t0\t1\tCANDIDATE\t\n"
        f"{PASS_KMER}\t18S\tearlier\t0\t1\t0\t0\t0\t1\tCANDIDATE\t\n"
    )
    entropy_pass = tmp_path / "entropy.fasta"
    entropy_pass.write_text(
        f">later\n{later_kmer}\n>earlier\n{PASS_KMER}\n"
    )
    final_baits = tmp_path / "final.fasta"

    assert finalize_manifest(
        raw_manifest,
        entropy_pass,
        tmp_path / "final.tsv",
        final_baits,
    ) == 2
    assert final_baits.read_text() == (
        f">cc_rrna_kmer_000001|18S\n{PASS_KMER}\n"
        f">cc_rrna_kmer_000002|28S\n{later_kmer}\n"
    )


def test_cli_exposes_raw_and_finalize_subcommands(tmp_path: Path):
    target_fasta = tmp_path / "targets.fasta"
    target_fasta.write_text(f">18S|synthetic\n{PASS_KMER}\n")
    loci_tsv = tmp_path / "loci.tsv"
    _write_loci(loci_tsv, [])
    empty = tmp_path / "empty.txt"
    _write_meryl_print(empty)
    exact = tmp_path / "exact.txt"
    _write_meryl_print(exact, PASS_KMER)
    raw_manifest = tmp_path / "raw.tsv"
    raw_baits = tmp_path / "raw.fasta"

    main(
        [
            "raw",
            "--target-fasta",
            str(target_fasta),
            "--loci-tsv",
            str(loci_tsv),
            "--silva-print",
            str(empty),
            "--rfam-print",
            str(empty),
            "--other-cyclospora-print",
            str(empty),
            "--exact-difference-print",
            str(exact),
            "--raw-manifest",
            str(raw_manifest),
            "--raw-baits",
            str(raw_baits),
        ]
    )
    entropy_pass = tmp_path / "entropy.fasta"
    entropy_pass.write_text(f">accepted\n{PASS_KMER}\n")
    final_manifest = tmp_path / "final.tsv"
    final_baits = tmp_path / "final.fasta"
    main(
        [
            "finalize",
            "--raw-manifest",
            str(raw_manifest),
            "--entropy-pass-fasta",
            str(entropy_pass),
            "--final-manifest",
            str(final_manifest),
            "--final-baits",
            str(final_baits),
        ]
    )

    assert "status" in final_manifest.read_text().splitlines()[0]
    assert final_baits.read_text().startswith(">cc_rrna_kmer_000001|18S\n")
