from pathlib import Path

import pytest
from Bio import SeqIO

from rrna_bait.targets import (
    BlastHit,
    TargetRule,
    accept_hit,
    curate_targets,
    prepare_queries,
)


def _write_fixture_rules(path: Path) -> None:
    path.write_text(
        "accession\trrna_class\tquery_start\tquery_end\tmin_identity\tmin_query_coverage\n"
        "AF111183.1\t18S\t1\t12\t98.0\t90.0\n"
        "XR_003297348.1\t5S\t1\t0\t98.0\t90.0\n"
    )


def test_accept_hit_enforces_identity_and_query_coverage():
    rule = TargetRule("Q1.1", "18S", 1, 0, 98.0, 90.0)
    accepted = BlastHit("Q1.1", 100, 1, 100, "REF1", 11, 110, 100, 99.0)
    short = BlastHit("Q1.1", 100, 1, 80, "REF1", 11, 90, 80, 100.0)
    assert accept_hit(accepted, rule)
    assert not accept_hit(short, rule)


def test_curate_targets_extracts_reverse_strand_and_deduplicates(tmp_path: Path):
    fasta_out = tmp_path / "targets.fasta"
    loci_out = tmp_path / "loci.tsv"
    rules_path = tmp_path / "rules.tsv"
    _write_fixture_rules(rules_path)
    curate_targets(
        queries_path=Path("tests/fixtures/target_queries.fasta"),
        subject_path=Path("tests/fixtures/refseq_subject.fasta"),
        rules_path=rules_path,
        hits_path=Path("tests/fixtures/target_hits.tsv"),
        fasta_out=fasta_out,
        loci_out=loci_out,
        expected_classes=None,
    )
    records = list(SeqIO.parse(fasta_out, "fasta"))
    assert [(record.id.split("|")[0], str(record.seq)) for record in records] == [
        ("18S", "ACGTACGTACGT"),
        ("5S", "AACCGGTT"),
    ]
    locus_rows = loci_out.read_text().splitlines()
    assert locus_rows[0].startswith("target_record_id\trrna_class\t")
    representative_18s = next(record.id for record in records if record.id.startswith("18S|"))
    supporting_18s = [row.split("\t") for row in locus_rows[1:] if "\t18S\t" in row]
    assert {row[0] for row in supporting_18s} == {representative_18s}
    assert {row[3] for row in supporting_18s} == {"REF_DUPLICATE", "REF_FORWARD"}
    assert any("REF_REVERSE\t5\t12\t-" in row for row in locus_rows)
    assert sum("18S\tAF111183.1" in row for row in locus_rows) == 2


def test_curate_targets_rejects_untrimmed_bounded_query(tmp_path: Path):
    rules_path = tmp_path / "rules.tsv"
    rules_path.write_text(
        "accession\trrna_class\tquery_start\tquery_end\tmin_identity\tmin_query_coverage\n"
        "MPGL01000046.1\t28S\t11\t3500\t98.0\t90.0\n"
    )
    queries_path = tmp_path / "queries.fasta"
    queries_path.write_text(">MPGL01000046.1\n" + "A" * 3500 + "\n")
    subject_path = tmp_path / "subjects.fasta"
    subject_path.write_text(">REF\n" + "A" * 3490 + "\n")
    hits_path = tmp_path / "hits.tsv"
    hits_path.write_text(
        "MPGL01000046.1\t3500\t1\t3490\tREF\t1\t3490\t3490\t100.0\n"
    )

    with pytest.raises(ValueError, match="MPGL01000046.1.*expected 3490.*got 3500"):
        curate_targets(
            queries_path=queries_path,
            subject_path=subject_path,
            rules_path=rules_path,
            hits_path=hits_path,
            fasta_out=tmp_path / "targets.fasta",
            loci_out=tmp_path / "loci.tsv",
            expected_classes=None,
        )


def test_curate_targets_rejects_blast_query_length_mismatch(tmp_path: Path):
    rules_path = tmp_path / "rules.tsv"
    _write_fixture_rules(rules_path)
    hits_path = tmp_path / "hits.tsv"
    hits_path.write_text(
        "AF111183.1\t11\t1\t11\tREF_FORWARD\t1\t11\t11\t100.0\n"
    )

    with pytest.raises(ValueError, match="BLAST query length.*AF111183.1.*expected 12.*got 11"):
        curate_targets(
            queries_path=Path("tests/fixtures/target_queries.fasta"),
            subject_path=Path("tests/fixtures/refseq_subject.fasta"),
            rules_path=rules_path,
            hits_path=hits_path,
            fasta_out=tmp_path / "targets.fasta",
            loci_out=tmp_path / "loci.tsv",
            expected_classes=None,
        )


def test_extract_loci_rejects_ambiguous_subject_interval():
    rule = TargetRule("Q1", "18S", 1, 0, 98.0, 90.0)
    subjects = list(SeqIO.parse(Path("tests/fixtures/refseq_subject.fasta"), "fasta"))
    subjects[0].seq = subjects[0].seq[:4] + "N" + subjects[0].seq[5:]
    hits = [BlastHit("Q1", 12, 1, 12, "REF_FORWARD", 1, 12, 12, 100.0)]

    from rrna_bait.targets import extract_loci

    assert extract_loci(subjects, hits, [rule]) == []


def test_curate_targets_output_is_stable_when_hits_are_reordered(tmp_path: Path):
    rules_path = tmp_path / "rules.tsv"
    _write_fixture_rules(rules_path)
    original_hits = Path("tests/fixtures/target_hits.tsv").read_text().splitlines()
    reordered_hits = tmp_path / "reordered_hits.tsv"
    reordered_hits.write_text("\n".join(reversed(original_hits)) + "\n")

    outputs: list[tuple[str, str]] = []
    for label, hits_path in (("original", Path("tests/fixtures/target_hits.tsv")), ("reordered", reordered_hits)):
        fasta_out = tmp_path / f"{label}.fasta"
        loci_out = tmp_path / f"{label}.tsv"
        curate_targets(
            queries_path=Path("tests/fixtures/target_queries.fasta"),
            subject_path=Path("tests/fixtures/refseq_subject.fasta"),
            rules_path=rules_path,
            hits_path=hits_path,
            fasta_out=fasta_out,
            loci_out=loci_out,
            expected_classes=None,
        )
        outputs.append((fasta_out.read_text(), loci_out.read_text()))

    assert outputs[0] == outputs[1]


def test_curate_targets_skips_unknown_query_and_missing_subject_hits(tmp_path: Path):
    rules_path = tmp_path / "rules.tsv"
    _write_fixture_rules(rules_path)
    queries_path = tmp_path / "queries.fasta"
    queries_path.write_text(
        Path("tests/fixtures/target_queries.fasta").read_text() + ">EXTRA.1\nACGTACGTACGT\n"
    )
    hits_path = tmp_path / "hits.tsv"
    hits_path.write_text(
        Path("tests/fixtures/target_hits.tsv").read_text()
        + "UNKNOWN.1\t12\t1\t12\tREF_FORWARD\t1\t12\t12\t100.0\n"
        + "EXTRA.1\t12\t1\t12\tREF_FORWARD\t1\t12\t12\t100.0\n"
        + "AF111183.1\t12\t1\t12\tMISSING_SUBJECT\t1\t12\t12\t100.0\n"
    )
    fasta_out = tmp_path / "targets.fasta"
    loci_out = tmp_path / "loci.tsv"

    curate_targets(
        queries_path=queries_path,
        subject_path=Path("tests/fixtures/refseq_subject.fasta"),
        rules_path=rules_path,
        hits_path=hits_path,
        fasta_out=fasta_out,
        loci_out=loci_out,
        expected_classes=None,
    )

    assert [(record.id.split("|")[0], str(record.seq)) for record in SeqIO.parse(fasta_out, "fasta")] == [
        ("18S", "ACGTACGTACGT"),
        ("5S", "AACCGGTT"),
    ]


def test_prepare_queries_trims_intervals_and_uses_configured_accession(
    tmp_path: Path,
):
    output_fasta = tmp_path / "queries.fasta"
    prepare_queries(
        Path("tests/fixtures/target_queries.fasta"),
        [
            TargetRule("AF111183.1", "18S", 2, 11, 98.0, 90.0),
            TargetRule("XR_003297348.1", "5S", 1, 0, 98.0, 90.0),
        ],
        output_fasta,
    )
    records = list(SeqIO.parse(output_fasta, "fasta"))
    assert [(record.id, str(record.seq)) for record in records] == [
        ("AF111183.1", "CGTACGTACG"),
        ("XR_003297348.1", "AACCGGTT"),
    ]


def test_prepare_queries_rejects_out_of_range_interval(tmp_path: Path):
    with pytest.raises(ValueError, match="out of range"):
        prepare_queries(
            Path("tests/fixtures/target_queries.fasta"),
            [TargetRule("AF111183.1", "18S", 1, 13, 98.0, 90.0)],
            tmp_path / "queries.fasta",
        )


def test_curate_targets_requires_expected_rrna_classes(tmp_path: Path):
    rules_path = tmp_path / "rules.tsv"
    _write_fixture_rules(rules_path)
    with pytest.raises(ValueError, match="Missing expected rRNA classes"):
        curate_targets(
            queries_path=Path("tests/fixtures/target_queries.fasta"),
            subject_path=Path("tests/fixtures/refseq_subject.fasta"),
            rules_path=rules_path,
            hits_path=Path("tests/fixtures/target_hits.tsv"),
            fasta_out=tmp_path / "targets.fasta",
            loci_out=tmp_path / "loci.tsv",
            expected_classes={"18S", "5.8S", "28S", "5S"},
        )
