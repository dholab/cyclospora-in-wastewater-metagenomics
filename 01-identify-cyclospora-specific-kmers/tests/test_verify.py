from __future__ import annotations

import csv
from pathlib import Path

import pytest

from rrna_bait.core import (
    canonical_kmer,
    iter_canonical_kmers,
    reverse_complement,
)
from rrna_bait.verify import main, parse_index_info, verify_static


MANIFEST_FIELDS = (
    "kmer",
    "rrna_classes",
    "target_records",
    "target_starts",
    "target_occurrences",
    "target_copy_count",
    "in_silva",
    "in_rfam",
    "in_other_cyclospora",
    "in_exact_difference",
    "status",
    "rejection_reason",
    "core_nt_bait_id",
    "core_nt_status",
    "core_nt_rejection_reason",
)
WORD_18S = "AACCGGTTACGATCGTAGCTAGGCTAACGTA"
WORD_28S = "AGTCGATGCTACCGTTAAGGCTACGATTCGA"
BAIT_18S = canonical_kmer(WORD_18S)
BAIT_28S = canonical_kmer(WORD_28S)


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.write_text(
        "".join(f">{record_id}\n{sequence}\n" for record_id, sequence in records)
    )


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=MANIFEST_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_meryl_print(path: Path, kmers: set[str]) -> None:
    path.write_text("".join(f"{kmer}\t1\n" for kmer in sorted(kmers)))


def manifest_row(
    kmer: str,
    rrna_class: str,
    *,
    status: str = "PASS",
    in_silva: str = "0",
    in_rfam: str = "0",
    in_other_cyclospora: str = "0",
    in_exact_difference: str = "1",
    rejection_reason: str = "none",
    core_nt_bait_id: str | None = None,
    core_nt_status: str | None = None,
    core_nt_rejection_reason: str | None = None,
) -> dict[str, str]:
    if core_nt_status is None:
        core_nt_status = "PASS_CORE_NT" if status == "PASS" else "NOT_APPLICABLE"
    if core_nt_bait_id is None:
        core_nt_bait_id = f"bait_{rrna_class}" if status == "PASS" else ""
    if core_nt_rejection_reason is None:
        core_nt_rejection_reason = (
            "none" if status == "PASS" else "not_screened_static_reject"
        )
    return {
        "kmer": kmer,
        "rrna_classes": rrna_class,
        "target_records": f"{rrna_class}|fixture_{rrna_class}",
        "target_starts": "20",
        "target_occurrences": f"{rrna_class}|fixture_{rrna_class}:20",
        "target_copy_count": "1",
        "in_silva": in_silva,
        "in_rfam": in_rfam,
        "in_other_cyclospora": in_other_cyclospora,
        "in_exact_difference": in_exact_difference,
        "status": status,
        "rejection_reason": rejection_reason,
        "core_nt_bait_id": core_nt_bait_id,
        "core_nt_status": core_nt_status,
        "core_nt_rejection_reason": core_nt_rejection_reason,
    }


@pytest.fixture
def passing_inputs(tmp_path: Path) -> dict[str, Path | str]:
    target = tmp_path / "target.fasta"
    manifest = tmp_path / "manifest.tsv"
    baits = tmp_path / "baits.fasta"
    exact_difference = tmp_path / "exact.tsv"
    silva_intersection = tmp_path / "silva.tsv"
    rfam_intersection = tmp_path / "rfam.tsv"
    other_cyclospora_intersection = tmp_path / "other.tsv"
    genus_compatible = tmp_path / "genus.tsv"
    loci = tmp_path / "loci.tsv"
    provenance = tmp_path / "sources.tsv"
    summary = tmp_path / "summary.tsv"
    report = tmp_path / "report.md"
    core_nt_decisions = tmp_path / "core_nt_decisions.tsv"

    target_records = [
        ("18S|fixture_18S", "T" * 20 + WORD_18S + "C" * 129),
        ("28S|fixture_28S", "G" * 20 + WORD_28S + "T" * 129),
    ]
    write_fasta(target, target_records)
    classes_by_kmer: dict[str, set[str]] = {}
    occurrences_by_kmer: dict[str, set[tuple[str, int]]] = {}
    for record_id, sequence in target_records:
        rrna_class = record_id.split("|", 1)[0]
        for start, kmer in iter_canonical_kmers(sequence):
            classes_by_kmer.setdefault(kmer, set()).add(rrna_class)
            occurrences_by_kmer.setdefault(kmer, set()).add((record_id, start))
    passing_kmers = {BAIT_18S, BAIT_28S}
    silva_kmers = set(classes_by_kmer) - passing_kmers
    write_manifest(
        manifest,
        [
            {
                **manifest_row(
                    kmer,
                    ",".join(sorted(classes_by_kmer[kmer])),
                    status=(
                        "PASS"
                        if kmer in passing_kmers
                        else "REJECT_STATIC_BACKGROUND"
                    ),
                    in_silva=str(int(kmer in silva_kmers)),
                    in_exact_difference=str(int(kmer in passing_kmers)),
                    rejection_reason=(
                        "none" if kmer in passing_kmers else "silva"
                    ),
                ),
                "target_records": ",".join(
                    sorted({
                        record_id
                        for record_id, _ in occurrences_by_kmer[kmer]
                    })
                ),
                "target_starts": ",".join(
                    str(start)
                    for start in sorted({
                        start for _, start in occurrences_by_kmer[kmer]
                    })
                ),
                "target_occurrences": ";".join(
                    f"{record_id}:{start}"
                    for record_id, start in sorted(occurrences_by_kmer[kmer])
                ),
                "target_copy_count": str(len(occurrences_by_kmer[kmer])),
            }
            for kmer in sorted(classes_by_kmer)
        ],
    )
    write_fasta(
        baits,
        [
            ("cc_rrna_kmer_000001|18S", BAIT_18S),
            ("cc_rrna_kmer_000002|28S", BAIT_28S),
        ],
    )
    exact_difference.write_text(
        f"{reverse_complement(BAIT_28S)}\t1\n"
        f"{reverse_complement(BAIT_18S)}\t1\n"
    )
    write_meryl_print(silva_intersection, silva_kmers)
    write_meryl_print(rfam_intersection, set())
    write_meryl_print(other_cyclospora_intersection, set())
    write_meryl_print(genus_compatible, passing_kmers)
    loci.write_text(
        "target_record_id\trrna_class\tquery_accession\tsubject_accession"
        "\tstart\tend\tstrand\n"
        "18S|fixture_18S\t18S\t18S\tfixture_18S\t1\t180\t+\n"
        "28S|fixture_28S\t28S\t28S\tfixture_28S\t1\t180\t+\n"
    )
    provenance.write_text(
        "name\trelease\turl\tlocal_path\tretrieved_at_utc\tbytes\tsha256\n"
        "fixture\t1\thttps://example.test/source\tfixture.fasta"
        "\t2026-07-23T12:00:00+00:00\t31\tabc123\n"
    )
    core_nt_decisions.write_text(
        "bait_id\tstatus\texact_target_count\texact_non_target_count"
        "\tnear_match_count\n"
        "bait_18S\tPASS_CORE_NT\t1\t0\t0\n"
        "bait_28S\tPASS_CORE_NT\t1\t0\t0\n"
    )
    return {
        "target_fasta": target,
        "manifest": manifest,
        "baits": baits,
        "exact_difference": exact_difference,
        "silva_intersection": silva_intersection,
        "rfam_intersection": rfam_intersection,
        "other_cyclospora_intersection": other_cyclospora_intersection,
        "genus_compatible": genus_compatible,
        "loci_tsv": loci,
        "provenance": provenance,
        "core_nt_decisions": core_nt_decisions,
        "summary_path": summary,
        "report_path": report,
        "index_info": (
            "K-mer length: 31\n"
            "Window size: 1\n"
            "Distinct minimizer count: 2\n"
        ),
    }


def test_verify_rejects_missing_core_nt_decision(
    passing_inputs: dict[str, Path | str],
) -> None:
    rows = read_manifest(Path(passing_inputs["manifest"]))
    next(row for row in rows if row["kmer"] == BAIT_18S)[
        "core_nt_status"
    ] = ""
    write_manifest(Path(passing_inputs["manifest"]), rows)

    with pytest.raises(ValueError, match="core-nt status"):
        verify_static(**passing_inputs)


def test_verify_rejects_non_target_exact_pass_bait(
    passing_inputs: dict[str, Path | str],
) -> None:
    rows = read_manifest(Path(passing_inputs["manifest"]))
    row = next(row for row in rows if row["kmer"] == BAIT_18S)
    row["core_nt_status"] = "REJECT_CORE_NT_EXACT_NON_TARGET"
    row["core_nt_rejection_reason"] = "exact_non_target_31_of_31"
    write_manifest(Path(passing_inputs["manifest"]), rows)
    Path(passing_inputs["core_nt_decisions"]).write_text(
        "bait_id\tstatus\texact_target_count\texact_non_target_count"
        "\tnear_match_count\n"
        "bait_18S\tREJECT_CORE_NT_EXACT_NON_TARGET\t0\t1\t0\n"
        "bait_28S\tPASS_CORE_NT\t1\t0\t0\n"
    )

    with pytest.raises(ValueError, match="non-target exact"):
        verify_static(**passing_inputs)


def test_verify_rejects_incomplete_core_nt_decision_table(
    passing_inputs: dict[str, Path | str],
) -> None:
    Path(passing_inputs["core_nt_decisions"]).write_text(
        "bait_id\tstatus\texact_target_count\texact_non_target_count"
        "\tnear_match_count\n"
        "bait_18S\tPASS_CORE_NT\t1\t0\t0\n"
    )

    with pytest.raises(ValueError, match="decision table.*complete"):
        verify_static(**passing_inputs)


@pytest.mark.parametrize("length", [30, 32])
def test_verify_rejects_non_31_base_bait(
    passing_inputs: dict[str, Path | str], length: int
) -> None:
    write_fasta(Path(passing_inputs["baits"]), [("bad_length", "A" * length)])

    with pytest.raises(ValueError, match="exactly 31"):
        verify_static(**passing_inputs)


def test_verify_rejects_ambiguous_bait(
    passing_inputs: dict[str, Path | str],
) -> None:
    write_fasta(
        Path(passing_inputs["baits"]),
        [("ambiguous", "A" * 15 + "N" + "A" * 15)],
    )

    with pytest.raises(ValueError, match="ACGT"):
        verify_static(**passing_inputs)


def test_verify_rejects_bait_absent_from_target(
    passing_inputs: dict[str, Path | str],
) -> None:
    absent = "A" * 31
    write_fasta(
        Path(passing_inputs["target_fasta"]),
        [("18S|fixture_18S", WORD_18S), ("28S|fixture_28S", WORD_28S)],
    )
    write_fasta(Path(passing_inputs["baits"]), [("absent", absent)])
    write_manifest(
        Path(passing_inputs["manifest"]),
        [manifest_row(absent, "18S")],
    )
    Path(passing_inputs["exact_difference"]).write_text(f"{absent}\t1\n")
    passing_inputs["index_info"] = (
        "K-mer length: 31\nWindow size: 1\nDistinct minimizers: 1\n"
    )

    with pytest.raises(ValueError, match="eligible target"):
        verify_static(**passing_inputs)


def test_verify_rejects_bait_marked_rejected_in_manifest(
    passing_inputs: dict[str, Path | str],
) -> None:
    write_manifest(
        Path(passing_inputs["manifest"]),
        [
            manifest_row(
                BAIT_18S,
                "18S",
                status="REJECT_LOW_COMPLEXITY",
                rejection_reason="low_complexity",
            ),
            manifest_row(BAIT_28S, "28S"),
        ],
    )

    with pytest.raises(ValueError, match="status=PASS"):
        verify_static(**passing_inputs)


def test_verify_rejects_duplicate_bait_records(
    passing_inputs: dict[str, Path | str],
) -> None:
    write_fasta(
        Path(passing_inputs["baits"]),
        [
            ("first", BAIT_18S),
            ("duplicate", BAIT_18S),
            ("second", BAIT_28S),
        ],
    )
    passing_inputs["index_info"] = (
        "K-mer length: 31\nWindow size: 1\nDistinct minimizers: 3\n"
    )

    with pytest.raises(ValueError, match="duplicate"):
        verify_static(**passing_inputs)


def test_verify_rejects_indexed_count_mismatch(
    passing_inputs: dict[str, Path | str],
) -> None:
    passing_inputs["index_info"] = (
        "K-mer length: 31\n"
        "Window size: 1\n"
        "Distinct minimizer count: 3\n"
    )

    with pytest.raises(ValueError, match="indexed minimizer count"):
        verify_static(**passing_inputs)


def test_verify_rejects_exact_difference_disagreement(
    passing_inputs: dict[str, Path | str],
) -> None:
    Path(passing_inputs["exact_difference"]).write_text(f"{BAIT_18S}\t1\n")

    with pytest.raises(ValueError, match="exact-difference"):
        verify_static(**passing_inputs)


def test_verify_requires_authoritative_exact_difference(
    passing_inputs: dict[str, Path | str],
) -> None:
    incomplete = dict(passing_inputs)
    del incomplete["exact_difference"]

    with pytest.raises(TypeError, match="exact_difference"):
        verify_static(**incomplete)


def test_cli_requires_all_authoritative_evidence(
    passing_inputs: dict[str, Path | str],
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as error:
        main(
            [
                "--target",
                str(passing_inputs["target_fasta"]),
                "--manifest",
                str(passing_inputs["manifest"]),
                "--baits",
                str(passing_inputs["baits"]),
                "--summary",
                str(passing_inputs["summary_path"]),
                "--report",
                str(passing_inputs["report_path"]),
            ]
        )

    assert error.value.code == 2
    cli_error = capsys.readouterr().err
    assert "--exact-difference" in cli_error
    assert "--silva-intersection" in cli_error
    assert "--rfam-intersection" in cli_error
    assert "--other-cyclospora-intersection" in cli_error
    assert "--genus-compatible" in cli_error
    assert "--loci" in cli_error
    assert "--provenance" in cli_error


def test_verify_requires_locus_metadata(
    passing_inputs: dict[str, Path | str],
) -> None:
    incomplete = dict(passing_inputs)
    del incomplete["loci_tsv"]

    with pytest.raises(TypeError, match="loci_tsv"):
        verify_static(**incomplete)


@pytest.mark.parametrize(
    ("replacement", "message"),
    [
        (
            "18S|unknown\t18S\t18S\tunknown\t1\t180\t+\n",
            "unknown target_record_id",
        ),
        (
            "18S|fixture_18S\t28S\t18S\tfixture_18S\t1\t180\t+\n"
            "28S|fixture_28S\t28S\t28S\tfixture_28S\t1\t180\t+\n",
            "rrna_class",
        ),
        (
            "18S|fixture_18S\t18S\t18S\tfixture_18S\t1\t180\t+\n"
            "18S|fixture_18S\t18S\t18S\tfixture_18S\t1\t180\t+\n"
            "28S|fixture_28S\t28S\t28S\tfixture_28S\t1\t180\t+\n",
            "duplicate locus",
        ),
    ],
)
def test_verify_rejects_invalid_locus_metadata(
    passing_inputs: dict[str, Path | str],
    replacement: str,
    message: str,
) -> None:
    Path(passing_inputs["loci_tsv"]).write_text(
        "target_record_id\trrna_class\tquery_accession\tsubject_accession"
        "\tstart\tend\tstrand\n"
        + replacement
    )

    with pytest.raises(ValueError, match=message):
        verify_static(**passing_inputs)


def test_verify_rejects_locus_metadata_missing_target_record(
    passing_inputs: dict[str, Path | str],
) -> None:
    Path(passing_inputs["loci_tsv"]).write_text(
        "target_record_id\trrna_class\tquery_accession\tsubject_accession"
        "\tstart\tend\tstrand\n"
        "18S|fixture_18S\t18S\t18S\tfixture_18S\t1\t180\t+\n"
    )

    with pytest.raises(ValueError, match="missing target_record_id"):
        verify_static(**passing_inputs)


@pytest.mark.parametrize("artifact", ["baits", "manifest"])
def test_verify_rejects_raw_u_before_normalization(
    passing_inputs: dict[str, Path | str], artifact: str
) -> None:
    raw_u = BAIT_18S.replace("T", "U", 1)
    assert raw_u != BAIT_18S
    if artifact == "baits":
        write_fasta(
            Path(passing_inputs["baits"]),
            [("raw_u", raw_u), ("other", BAIT_28S)],
        )
    else:
        rows = read_manifest(Path(passing_inputs["manifest"]))
        next(row for row in rows if row["kmer"] == BAIT_18S)["kmer"] = raw_u
        write_manifest(Path(passing_inputs["manifest"]), rows)

    with pytest.raises(ValueError, match="uppercase ACGT"):
        verify_static(**passing_inputs)


@pytest.mark.parametrize(
    ("field", "evidence_name"),
    [
        ("in_silva", "SILVA"),
        ("in_rfam", "Rfam"),
        ("in_other_cyclospora", "other-Cyclospora"),
        ("in_exact_difference", "exact-difference"),
    ],
)
def test_verify_rejects_manifest_flag_disagreement(
    passing_inputs: dict[str, Path | str],
    field: str,
    evidence_name: str,
) -> None:
    rows = read_manifest(Path(passing_inputs["manifest"]))
    row = next(row for row in rows if row["kmer"] == BAIT_18S)
    row[field] = str(1 - int(row[field]))
    write_manifest(Path(passing_inputs["manifest"]), rows)

    with pytest.raises(ValueError, match=evidence_name):
        verify_static(**passing_inputs)


def test_verify_rejects_genus_compatible_disagreement(
    passing_inputs: dict[str, Path | str],
) -> None:
    Path(passing_inputs["genus_compatible"]).write_text("")

    with pytest.raises(ValueError, match="genus-compatible"):
        verify_static(**passing_inputs)


def test_verify_rejects_rrna_class_annotation_disagreement(
    passing_inputs: dict[str, Path | str],
) -> None:
    rows = read_manifest(Path(passing_inputs["manifest"]))
    next(row for row in rows if row["kmer"] == BAIT_18S)[
        "rrna_classes"
    ] = "5S"
    write_manifest(Path(passing_inputs["manifest"]), rows)

    with pytest.raises(ValueError, match="rrna_classes"):
        verify_static(**passing_inputs)


@pytest.mark.parametrize(
    ("field", "replacement", "message"),
    [
        ("target_records", "wrong_record", "target_records"),
        ("target_starts", "999", "target_starts"),
        ("target_copy_count", "99", "target_copy_count"),
        ("target_occurrences", "18S|fixture_18S:999", "target_occurrences"),
        ("target_occurrences", "wrong_record:20", "target_occurrences"),
    ],
)
def test_verify_rejects_corrupt_target_occurrence_metadata(
    passing_inputs: dict[str, Path | str],
    field: str,
    replacement: str,
    message: str,
) -> None:
    rows = read_manifest(Path(passing_inputs["manifest"]))
    row = next(row for row in rows if row["kmer"] == BAIT_18S)
    row[field] = replacement
    write_manifest(Path(passing_inputs["manifest"]), rows)

    with pytest.raises(ValueError, match=message):
        verify_static(**passing_inputs)


@pytest.mark.parametrize(
    ("status", "reason", "message"),
    [
        ("REJECT_LOW_COMPLEXITY", "low_complexity", "status"),
        ("PASS", "silva", "rejection_reason"),
    ],
)
def test_verify_derives_status_and_reason_from_authoritative_evidence(
    passing_inputs: dict[str, Path | str],
    status: str,
    reason: str,
    message: str,
) -> None:
    rows = read_manifest(Path(passing_inputs["manifest"]))
    row = next(row for row in rows if row["kmer"] == BAIT_18S)
    row["status"] = status
    row["rejection_reason"] = reason
    write_manifest(Path(passing_inputs["manifest"]), rows)

    with pytest.raises(ValueError, match=message):
        verify_static(**passing_inputs)


def test_verify_rejects_nonbinary_manifest_flag(
    passing_inputs: dict[str, Path | str],
) -> None:
    rows = read_manifest(Path(passing_inputs["manifest"]))
    rows[0]["in_silva"] = "2"
    write_manifest(Path(passing_inputs["manifest"]), rows)

    with pytest.raises(ValueError, match="binary"):
        verify_static(**passing_inputs)


def test_verify_rejects_unknown_manifest_status(
    passing_inputs: dict[str, Path | str],
) -> None:
    rows = read_manifest(Path(passing_inputs["manifest"]))
    rows[0]["status"] = "CANDIDATE"
    write_manifest(Path(passing_inputs["manifest"]), rows)

    with pytest.raises(ValueError, match="status"):
        verify_static(**passing_inputs)


def test_verify_rejects_status_exact_difference_inconsistency(
    passing_inputs: dict[str, Path | str],
) -> None:
    rows = read_manifest(Path(passing_inputs["manifest"]))
    row = next(row for row in rows if row["status"] == "REJECT_STATIC_BACKGROUND")
    row["status"] = "REJECT_LOW_COMPLEXITY"
    row["rejection_reason"] = "low_complexity"
    write_manifest(Path(passing_inputs["manifest"]), rows)

    with pytest.raises(ValueError, match="status.*exact-difference"):
        verify_static(**passing_inputs)


@pytest.mark.parametrize(
    "retrieved_at_utc",
    [
        "",
        "not-a-date",
        "2026-07-23T12:00:00",
        "2026-07-23T07:00:00-05:00",
    ],
)
def test_verify_rejects_missing_or_malformed_provenance_date(
    passing_inputs: dict[str, Path | str],
    retrieved_at_utc: str,
) -> None:
    Path(passing_inputs["provenance"]).write_text(
        "name\trelease\turl\tlocal_path\tretrieved_at_utc\tbytes\tsha256\n"
        "fixture\t1\thttps://example.test/source\tfixture.fasta"
        f"\t{retrieved_at_utc}\t31\tabc123\n"
    )

    with pytest.raises(ValueError, match="retrieved_at_utc"):
        verify_static(**passing_inputs)


def test_verify_rejects_passing_bait_not_marked_exact_difference(
    passing_inputs: dict[str, Path | str],
) -> None:
    rows = read_manifest(Path(passing_inputs["manifest"]))
    next(row for row in rows if row["kmer"] == BAIT_18S)[
        "in_exact_difference"
    ] = "0"
    write_manifest(Path(passing_inputs["manifest"]), rows)

    with pytest.raises(ValueError, match="exact-difference"):
        verify_static(**passing_inputs)


def test_verify_rejects_incomplete_target_manifest(
    passing_inputs: dict[str, Path | str],
) -> None:
    rows = read_manifest(Path(passing_inputs["manifest"]))
    removed = next(row for row in rows if row["status"] != "PASS")
    rows.remove(removed)
    write_manifest(Path(passing_inputs["manifest"]), rows)

    with pytest.raises(ValueError, match="manifest k-mers.*complete target"):
        verify_static(**passing_inputs)


def test_parse_index_info_accepts_deacon_015_stderr_format() -> None:
    assert parse_index_info(
        "Index information:\n"
        "  Format version: 3\n"
        "  K-mer length (k): 31\n"
        "  Window size (w): 1\n"
        "  Distinct minimizer count: 27\n"
        "Loaded index info in 71.47µs\n"
    ) == (31, 1, 27)


def test_verify_writes_deterministic_success_reports_and_coverage(
    passing_inputs: dict[str, Path | str],
) -> None:
    summary = verify_static(**passing_inputs)
    first_summary = Path(passing_inputs["summary_path"]).read_bytes()
    first_report = Path(passing_inputs["report_path"]).read_bytes()

    assert summary.bait_count == 2
    assert summary.k == 31
    assert summary.w == 1
    assert summary.indexed_count == 2
    assert {item.rrna_class for item in summary.classes} == {
        "18S",
        "5.8S",
        "28S",
        "5S",
    }
    assert {locus.rrna_class for locus in summary.loci} == {"18S", "28S"}
    assert all(locus.candidate_count > 1 for locus in summary.loci)
    assert all(locus.pass_count == 1 for locus in summary.loci)
    assert all(locus.bases_covered == 31 for locus in summary.loci)
    assert all(locus.read_starts_total == 31 for locus in summary.loci)
    assert all(locus.read_starts_covered == 21 for locus in summary.loci)

    summary_text = first_summary.decode()
    report_text = first_report.decode()
    assert summary_text.startswith("key\tvalue\n")
    assert (
        "conclusion\tPASS: A species-specific mature-rRNA k=31 bait set "
        "was produced.\n"
    ) in summary_text
    assert f"attrition.silva_shared_count\t{summary.candidate_count - 2}\n" in summary_text
    assert "attrition.rfam_shared_count\t0\n" in summary_text
    assert "class.18S.silva_shared_count\t" in summary_text
    assert "genus.genus_level_status\tSPECIES_SPECIFIC_CANDIDATES\n" in summary_text
    assert "## Provenance" in report_text
    assert "## Stage attrition" in report_text
    assert "shared with SILVA" in report_text
    assert "shared with Rfam" in report_text
    assert "## Per-locus coverage" in report_text
    assert "## Other-*Cyclospora* sharing" in report_text
    assert "## Final index properties" in report_text
    assert "retrieved UTC" in report_text
    assert "2026-07-23T12:00:00+00:00" in report_text
    assert "Wastewater performance has not been evaluated." in report_text

    verify_static(**passing_inputs)
    assert Path(passing_inputs["summary_path"]).read_bytes() == first_summary
    assert Path(passing_inputs["report_path"]).read_bytes() == first_report


def test_verify_reports_infeasible_without_index(tmp_path: Path) -> None:
    target = tmp_path / "target.fasta"
    manifest = tmp_path / "manifest.tsv"
    baits = tmp_path / "baits.fasta"
    exact_difference = tmp_path / "exact.tsv"
    silva_intersection = tmp_path / "silva.tsv"
    rfam_intersection = tmp_path / "rfam.tsv"
    other_intersection = tmp_path / "other.tsv"
    genus_compatible = tmp_path / "genus.tsv"
    loci = tmp_path / "loci.tsv"
    provenance = tmp_path / "sources.tsv"
    core_nt_decisions = tmp_path / "core_nt_decisions.tsv"
    summary_path = tmp_path / "summary.tsv"
    report_path = tmp_path / "report.md"
    write_fasta(target, [("18S|fixture", WORD_18S)])
    write_manifest(
        manifest,
        [
                {
                    **manifest_row(
                        BAIT_18S,
                        "18S",
                        status="REJECT_STATIC_BACKGROUND",
                        in_other_cyclospora="1",
                        in_exact_difference="0",
                        rejection_reason="other_cyclospora",
                    ),
                    "target_records": "18S|fixture",
                    "target_starts": "0",
                    "target_occurrences": "18S|fixture:0",
                }
        ],
    )
    baits.write_text("")
    exact_difference.write_text("")
    write_meryl_print(silva_intersection, set())
    write_meryl_print(rfam_intersection, set())
    write_meryl_print(other_intersection, {BAIT_18S})
    write_meryl_print(genus_compatible, {BAIT_18S})
    loci.write_text(
        "target_record_id\trrna_class\tquery_accession\tsubject_accession"
        "\tstart\tend\tstrand\n"
        "18S|fixture\t18S\t18S\tfixture\t1\t31\t+\n"
    )
    provenance.write_text(
        "name\trelease\turl\tlocal_path\tretrieved_at_utc\tbytes\tsha256\n"
    )
    core_nt_decisions.write_text(
        "bait_id\tstatus\texact_target_count\texact_non_target_count"
        "\tnear_match_count\n"
    )

    result = verify_static(
        target_fasta=target,
        manifest=manifest,
        baits=baits,
        index_info=None,
        exact_difference=exact_difference,
        silva_intersection=silva_intersection,
        rfam_intersection=rfam_intersection,
        other_cyclospora_intersection=other_intersection,
        genus_compatible=genus_compatible,
        loci_tsv=loci,
        provenance=provenance,
        core_nt_decisions=core_nt_decisions,
        allow_empty_provenance=True,
        summary_path=summary_path,
        report_path=report_path,
    )

    assert result.bait_count == 0
    assert result.indexed_count == 0
    assert (
        "conclusion\tINFEASIBLE: No species-specific mature-rRNA k=31 bait "
        "survived; no index was emitted.\n"
    ) in summary_path.read_text()
    assert "genus.genus_compatible_candidate_count\t1\n" in summary_path.read_text()
    assert "genus.genus_level_status\tGENUS_ONLY_CANDIDATES\n" in summary_path.read_text()
    assert "## Genus-level feasibility" in report_path.read_text()
    assert "GENUS_ONLY_CANDIDATES" in report_path.read_text()
    assert "index emitted | no" in report_path.read_text()
    assert "Wastewater performance has not been evaluated." in report_path.read_text()
