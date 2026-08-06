from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from collections import defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from rrna_bait.core import canonical_kmer, iter_canonical_kmers, normalize_dna
from rrna_bait.core import sha256_file


K = 31
W = 1
READ_LENGTH = 150
PASS_CONCLUSION = (
    "PASS: A species-specific mature-rRNA k=31 bait set was produced."
)
INFEASIBLE_CONCLUSION = (
    "INFEASIBLE: No species-specific mature-rRNA k=31 bait survived; "
    "no index was emitted."
)
_EXPECTED_CLASSES = ("18S", "5.8S", "28S", "5S")
_BINARY_FIELDS = (
    "in_silva",
    "in_rfam",
    "in_other_cyclospora",
    "in_exact_difference",
)
_ALLOWED_STATUSES = {
    "PASS",
    "REJECT_STATIC_BACKGROUND",
    "REJECT_LOW_COMPLEXITY",
    "REJECT_CORE_NT_EXACT_NON_TARGET",
}
_CORE_NT_FIELDS = {
    "core_nt_bait_id",
    "core_nt_status",
    "core_nt_rejection_reason",
}
_CORE_NT_STATUSES = {
    "PASS_CORE_NT",
    "REJECT_CORE_NT_EXACT_NON_TARGET",
    "NOT_APPLICABLE",
}
_PROVENANCE_FIELDS = (
    "name",
    "release",
    "url",
    "local_path",
    "retrieved_at_utc",
    "bytes",
    "sha256",
)
_LOCUS_FIELDS = (
    "target_record_id",
    "rrna_class",
    "query_accession",
    "subject_accession",
    "start",
    "end",
    "strand",
)


@dataclass(frozen=True)
class ClassSummary:
    rrna_class: str
    candidate_count: int
    silva_shared: int
    rfam_shared: int
    exact_specific_count: int
    low_complexity_rejected: int
    pass_count: int
    other_cyclospora_shared: int


@dataclass(frozen=True)
class LocusSummary:
    target_record_id: str
    rrna_class: str
    length: int
    candidate_count: int
    pass_count: int
    bases_covered: int
    read_starts_total: int
    read_starts_covered: int


@dataclass(frozen=True)
class VerificationSummary:
    bait_count: int
    k: int
    w: int
    indexed_count: int
    candidate_count: int
    silva_shared: int
    rfam_shared: int
    exact_specific_count: int
    low_complexity_rejected: int
    other_cyclospora_shared: int
    genus_compatible_count: int
    classes: tuple[ClassSummary, ...]
    loci: tuple[LocusSummary, ...]

    @property
    def conclusion(self) -> str:
        return PASS_CONCLUSION if self.bait_count else INFEASIBLE_CONCLUSION

    @property
    def genus_level_status(self) -> str:
        if self.genus_compatible_count == 0:
            return "ZERO_GENUS_CANDIDATES"
        if self.bait_count:
            return "SPECIES_SPECIFIC_CANDIDATES"
        return "GENUS_ONLY_CANDIDATES"


def _read_target_records(path: Path) -> list[SeqRecord]:
    records = list(SeqIO.parse(Path(path), "fasta"))
    if not records:
        raise ValueError(f"target FASTA has no records: {path}")
    record_ids = [record.id for record in records]
    if len(record_ids) != len(set(record_ids)):
        raise ValueError("target FASTA contains duplicate record identifiers")
    return records


def _target_kmers(
    records: Iterable[SeqRecord],
) -> tuple[
    set[str],
    dict[str, list[tuple[int, str]]],
    dict[str, set[str]],
    dict[str, set[tuple[str, int]]],
]:
    all_kmers: set[str] = set()
    by_record: dict[str, list[tuple[int, str]]] = {}
    classes_by_kmer: dict[str, set[str]] = defaultdict(set)
    occurrences_by_kmer: dict[str, set[tuple[str, int]]] = defaultdict(set)
    for record in records:
        occurrences = list(iter_canonical_kmers(str(record.seq), K))
        by_record[record.id] = occurrences
        all_kmers.update(kmer for _, kmer in occurrences)
        rrna_class = record.id.split("|", 1)[0]
        for start, kmer in occurrences:
            classes_by_kmer[kmer].add(rrna_class)
            occurrences_by_kmer[kmer].add((record.id, start))
    return (
        all_kmers,
        by_record,
        dict(classes_by_kmer),
        dict(occurrences_by_kmer),
    )


def _validate_raw_kmer(sequence: str, label: str) -> str:
    if len(sequence) != K:
        raise ValueError(
            f"{label} must be exactly {K} bases; found {len(sequence)}"
        )
    if re.fullmatch(r"[ACGT]+", sequence) is None:
        raise ValueError(f"{label} must contain exactly uppercase ACGT bases")
    if sequence != canonical_kmer(sequence):
        raise ValueError(f"{label} is not in canonical orientation")
    return sequence


def _read_baits(path: Path) -> tuple[list[str], set[str]]:
    sequences: list[str] = []
    for record in SeqIO.parse(Path(path), "fasta"):
        sequence = _validate_raw_kmer(
            str(record.seq), f"bait {record.id!r}"
        )
        sequences.append(sequence)

    unique = set(sequences)
    if len(sequences) != len(unique):
        raise ValueError("bait FASTA contains a duplicate bait sequence")
    return sequences, unique


def _read_manifest(path: Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
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
    }
    fields = set(rows[0]) if rows else set()
    if not rows:
        with Path(path).open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = set(reader.fieldnames or ())
    missing = required - fields
    if missing:
        raise ValueError(
            "manifest is missing required columns: " + ", ".join(sorted(missing))
        )
    present_core_fields = fields & _CORE_NT_FIELDS
    if present_core_fields and present_core_fields != _CORE_NT_FIELDS:
        raise ValueError(
            "manifest is missing required core-nt columns: "
            + ", ".join(sorted(_CORE_NT_FIELDS - fields))
        )

    kmers = [
        _validate_raw_kmer(row["kmer"], "manifest k-mer")
        for row in rows
    ]
    if len(kmers) != len(set(kmers)):
        raise ValueError("manifest contains a duplicate k-mer row")
    for row, kmer in zip(rows, kmers, strict=True):
        row["kmer"] = kmer
        for field in _BINARY_FIELDS:
            if row[field] not in {"0", "1"}:
                raise ValueError(
                    f"manifest field {field} must be binary 0 or 1"
                )
        status = row["status"]
        if status not in _ALLOWED_STATUSES:
            raise ValueError(f"manifest contains unsupported status: {status!r}")
        in_exact_difference = row["in_exact_difference"] == "1"
        if (not in_exact_difference and status != "REJECT_STATIC_BACKGROUND") or (
            in_exact_difference
            and status
            not in {
                "PASS",
                "REJECT_LOW_COMPLEXITY",
                "REJECT_CORE_NT_EXACT_NON_TARGET",
            }
        ):
            raise ValueError(
                "manifest status is inconsistent with exact-difference "
                f"membership for {kmer}"
            )
    return rows


def _read_core_nt_decisions(path: Path) -> dict[str, dict[str, str]]:
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        rows = list(reader)
    required = {
        "bait_id",
        "status",
        "exact_target_count",
        "exact_non_target_count",
        "near_match_count",
    }
    missing = required - fields
    if missing:
        raise ValueError(
            "core-nt decision table is missing required columns: "
            + ", ".join(sorted(missing))
        )
    decisions: dict[str, dict[str, str]] = {}
    for row_number, row in enumerate(rows, start=2):
        bait_id = row["bait_id"]
        if not bait_id:
            raise ValueError(
                f"core-nt decision table row {row_number} has an empty bait_id"
            )
        if bait_id in decisions:
            raise ValueError(
                "core-nt decision table contains duplicate bait_id: " + bait_id
            )
        if row["status"] not in _CORE_NT_STATUSES - {"NOT_APPLICABLE"}:
            raise ValueError(
                "core-nt decision table contains unknown status: "
                + repr(row["status"])
            )
        try:
            counts = [
                int(row[field])
                for field in (
                    "exact_target_count",
                    "exact_non_target_count",
                    "near_match_count",
                )
            ]
        except ValueError as error:
            raise ValueError(
                f"core-nt decision table row {row_number} has a non-integer count"
            ) from error
        if any(value < 0 for value in counts):
            raise ValueError(
                f"core-nt decision table row {row_number} has a negative count"
            )
        if (
            row["status"] == "PASS_CORE_NT"
            and int(row["exact_non_target_count"]) != 0
        ) or (
            row["status"] == "REJECT_CORE_NT_EXACT_NON_TARGET"
            and int(row["exact_non_target_count"]) == 0
        ):
            raise ValueError(
                f"core-nt decision table status/count mismatch for {bait_id}"
            )
        decisions[bait_id] = row
    return decisions


def _validate_core_nt_annotations(
    rows: list[dict[str, str]], core_nt_decisions: Path | None
) -> None:
    has_annotations = bool(rows and _CORE_NT_FIELDS <= set(rows[0]))
    if not has_annotations:
        if core_nt_decisions is not None:
            raise ValueError(
                "core-nt decision table was supplied but manifest lacks annotations"
            )
        return
    if core_nt_decisions is None:
        raise ValueError(
            "core-nt annotated manifest requires a core-nt decision table"
        )
    decisions = _read_core_nt_decisions(Path(core_nt_decisions))
    screened: dict[str, dict[str, str]] = {}
    for row in rows:
        core_status = row["core_nt_status"]
        if core_status not in _CORE_NT_STATUSES:
            raise ValueError(
                f"manifest core-nt status is missing or unknown for {row['kmer']}: "
                f"{core_status!r}"
            )
        bait_id = row["core_nt_bait_id"]
        if core_status == "NOT_APPLICABLE":
            if bait_id:
                raise ValueError(
                    f"NOT_APPLICABLE core-nt row has bait_id for {row['kmer']}"
                )
            if row["status"] in {"PASS", "REJECT_CORE_NT_EXACT_NON_TARGET"}:
                raise ValueError(
                    f"passing/core-nt-rejected row was not screened: {row['kmer']}"
                )
            continue
        if not bait_id:
            raise ValueError(
                f"manifest core-nt decision is missing bait_id for {row['kmer']}"
            )
        if bait_id in screened:
            raise ValueError(
                "manifest contains duplicate core-nt bait_id: " + bait_id
            )
        screened[bait_id] = row
        decision = decisions.get(bait_id)
        if decision is not None and decision["status"] != core_status:
            raise ValueError(
                f"manifest core-nt status disagrees with decision for {bait_id}"
            )
        if core_status == "REJECT_CORE_NT_EXACT_NON_TARGET":
            if row["status"] == "PASS":
                raise ValueError(
                    f"PASS bait has a non-target exact core-nt hit: {bait_id}"
                )
            if row["status"] != "REJECT_CORE_NT_EXACT_NON_TARGET":
                raise ValueError(
                    f"non-target exact core-nt decision has wrong status: {bait_id}"
                )
            if row["core_nt_rejection_reason"] != "exact_non_target_31_of_31":
                raise ValueError(
                    f"core-nt rejection reason is invalid for {bait_id}"
                )
        elif row["status"] != "PASS":
            raise ValueError(
                f"PASS_CORE_NT bait does not have status=PASS: {bait_id}"
            )
        elif row["core_nt_rejection_reason"] != "none":
            raise ValueError(
                f"PASS_CORE_NT bait has a rejection reason: {bait_id}"
            )

    if set(screened) != set(decisions):
        missing = sorted(set(screened) - set(decisions))
        extra = sorted(set(decisions) - set(screened))
        detail = []
        if missing:
            detail.append(f"missing {missing[0]}")
        if extra:
            detail.append(f"unexpected {extra[0]}")
        raise ValueError(
            "core-nt decision table does not completely and exactly cover "
            "the input baits"
            + (": " + "; ".join(detail) if detail else "")
        )


def _read_meryl_print(path: Path) -> set[str]:
    kmers: set[str] = set()
    with Path(path).open() as handle:
        for line in handle:
            fields = line.split()
            if fields:
                kmers.add(canonical_kmer(fields[0]))
    return kmers


def _require_manifest_membership(
    rows: list[dict[str, str]],
    field: str,
    evidence: set[str],
    label: str,
) -> None:
    annotated = {row["kmer"] for row in rows if row[field] == "1"}
    if annotated != evidence:
        missing = sorted(evidence - annotated)
        extra = sorted(annotated - evidence)
        detail = []
        if missing:
            detail.append(f"missing annotation for {missing[0]}")
        if extra:
            detail.append(f"unexpected annotation for {extra[0]}")
        raise ValueError(
            f"manifest {field} flags disagree with authoritative {label} "
            "Meryl printout"
            + (": " + "; ".join(detail) if detail else "")
        )


def _parse_target_occurrences(value: str, kmer: str) -> set[tuple[str, int]]:
    occurrences: set[tuple[str, int]] = set()
    for item in value.split(";"):
        try:
            record_id, raw_start = item.rsplit(":", 1)
            start = int(raw_start)
        except (ValueError, TypeError) as error:
            raise ValueError(
                f"manifest target_occurrences is malformed for {kmer}: {item!r}"
            ) from error
        occurrence = (record_id, start)
        if not record_id or start < 0 or occurrence in occurrences:
            raise ValueError(
                f"manifest target_occurrences is malformed for {kmer}: {item!r}"
            )
        occurrences.add(occurrence)
    return occurrences


def _validate_manifest_annotations(
    rows: list[dict[str, str]],
    target_kmers: set[str],
    classes_by_kmer: dict[str, set[str]],
    occurrences_by_kmer: dict[str, set[tuple[str, int]]],
    silva_intersection: set[str],
    rfam_intersection: set[str],
    other_cyclospora_intersection: set[str],
    exact_difference: set[str],
    genus_compatible: set[str],
    passing: set[str],
) -> None:
    for row in rows:
        kmer = row["kmer"]
        expected_classes = ",".join(sorted(classes_by_kmer[row["kmer"]]))
        if row["rrna_classes"] != expected_classes:
            raise ValueError(
                "manifest rrna_classes disagree with target FASTA "
                f"occurrences for {row['kmer']}: "
                f"{row['rrna_classes']!r} != {expected_classes!r}"
            )
        expected_occurrences = occurrences_by_kmer[kmer]
        expected_records = ",".join(
            sorted({record_id for record_id, _ in expected_occurrences})
        )
        expected_starts = ",".join(
            str(start)
            for start in sorted({start for _, start in expected_occurrences})
        )
        expected_metadata = {
            "target_records": expected_records,
            "target_starts": expected_starts,
            "target_copy_count": str(len(expected_occurrences)),
        }
        for field, expected in expected_metadata.items():
            if row[field] != expected:
                raise ValueError(
                    f"manifest {field} disagrees with target FASTA "
                    f"occurrences for {kmer}: {row[field]!r} != {expected!r}"
                )
        observed_occurrences = _parse_target_occurrences(
            row["target_occurrences"], kmer
        )
        if observed_occurrences != expected_occurrences:
            raise ValueError(
                "manifest target_occurrences disagrees with target FASTA "
                f"occurrences for {kmer}"
            )

    _require_manifest_membership(
        rows, "in_silva", silva_intersection, "SILVA"
    )
    _require_manifest_membership(
        rows, "in_rfam", rfam_intersection, "Rfam"
    )
    _require_manifest_membership(
        rows,
        "in_other_cyclospora",
        other_cyclospora_intersection,
        "other-Cyclospora",
    )
    _require_manifest_membership(
        rows,
        "in_exact_difference",
        exact_difference,
        "exact-difference",
    )

    expected_exact = target_kmers - (
        silva_intersection
        | rfam_intersection
        | other_cyclospora_intersection
    )
    if exact_difference != expected_exact:
        raise ValueError(
            "authoritative exact-difference Meryl printout does not equal "
            "target minus all background intersections"
        )
    expected_genus = target_kmers - (
        silva_intersection | rfam_intersection
    )
    if genus_compatible != expected_genus:
        raise ValueError(
            "authoritative genus-compatible Meryl printout does not equal "
            "target minus SILVA and Rfam intersections"
        )

    for row in rows:
        kmer = row["kmer"]
        if kmer not in exact_difference:
            expected_status = "REJECT_STATIC_BACKGROUND"
            expected_reason = ";".join(
                label
                for label, evidence in (
                    ("other_cyclospora", other_cyclospora_intersection),
                    ("silva", silva_intersection),
                    ("rfam", rfam_intersection),
                )
                if kmer in evidence
            )
        elif row.get("core_nt_status") == "REJECT_CORE_NT_EXACT_NON_TARGET":
            expected_status = "REJECT_CORE_NT_EXACT_NON_TARGET"
            expected_reason = "core_nt_exact_non_target"
        elif kmer in passing:
            expected_status = "PASS"
            expected_reason = "none"
        else:
            expected_status = "REJECT_LOW_COMPLEXITY"
            expected_reason = "low_complexity"
        if row["status"] != expected_status:
            raise ValueError(
                f"manifest status disagrees with authoritative evidence for {kmer}: "
                f"{row['status']!r} != {expected_status!r}"
            )
        if row.get("rejection_reason", "") != expected_reason:
            raise ValueError(
                "manifest rejection_reason disagrees with authoritative "
                f"evidence for {kmer}: {row.get('rejection_reason', '')!r} "
                f"!= {expected_reason!r}"
            )


def parse_index_info(index_info: str) -> tuple[int, int, int]:
    patterns = {
        "k": r"\bK-mer (?:length|size)(?:\s*\(k\))?\s*:\s*(\d+)",
        "w": r"\bWindow size(?:\s*\(w\))?\s*:\s*(\d+)",
        "count": r"\bDistinct minimizer(?:s| count)?\s*:\s*(\d+)",
    }
    values: dict[str, int] = {}
    for name, pattern in patterns.items():
        match = re.search(pattern, index_info, flags=re.IGNORECASE)
        if match is None:
            raise ValueError(f"Deacon index info is missing {name!r} metadata")
        values[name] = int(match.group(1))
    return values["k"], values["w"], values["count"]


def read_index_info(index: Path) -> str:
    index = Path(index)
    if not index.is_file():
        raise ValueError(f"Deacon index is absent: {index}")
    completed = subprocess.run(
        ["deacon", "index", "info", str(index)],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout).strip()
        raise ValueError(f"deacon index info failed: {detail}")
    return completed.stdout + completed.stderr


def _locus_classes(
    records: Iterable[SeqRecord], loci_tsv: Path
) -> dict[str, str]:
    expected = {
        record.id: record.id.split("|", 1)[0]
        for record in records
    }
    path = Path(loci_tsv)
    if not path.is_file():
        raise ValueError(f"locus metadata is absent: {path}")
    with Path(loci_tsv).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = tuple(reader.fieldnames or ())
        rows = list(reader)
    missing_fields = set(_LOCUS_FIELDS) - set(fields)
    if missing_fields:
        raise ValueError(
            "locus metadata is missing required columns: "
            + ", ".join(sorted(missing_fields))
        )

    seen_rows: set[tuple[str, ...]] = set()
    seen_ids: set[str] = set()
    for row_number, row in enumerate(rows, start=2):
        values = tuple(row.get(field, "") for field in _LOCUS_FIELDS)
        if any(not value for value in values):
            raise ValueError(
                f"locus metadata row {row_number} has an empty required field"
            )
        if values in seen_rows:
            raise ValueError(
                f"locus metadata contains a duplicate locus row at {row_number}"
            )
        seen_rows.add(values)
        record_id = row["target_record_id"]
        if record_id not in expected:
            raise ValueError(
                "locus metadata contains unknown target_record_id: "
                f"{record_id}"
            )
        if row["rrna_class"] != expected[record_id]:
            raise ValueError(
                "locus metadata rrna_class disagrees with target FASTA "
                f"record {record_id}: {row['rrna_class']!r} != "
                f"{expected[record_id]!r}"
            )
        seen_ids.add(record_id)

    missing_ids = sorted(set(expected) - seen_ids)
    if missing_ids:
        raise ValueError(
            "locus metadata is missing target_record_id: " + missing_ids[0]
        )
    return expected


def _row_classes(row: dict[str, str]) -> tuple[str, ...]:
    return tuple(
        sorted(
            {
                value.strip()
                for value in row["rrna_classes"].split(",")
                if value.strip()
            }
        )
    )


def _class_summaries(rows: list[dict[str, str]]) -> tuple[ClassSummary, ...]:
    by_class: dict[str, dict[str, int]] = defaultdict(
        lambda: {
            "candidate": 0,
            "silva": 0,
            "rfam": 0,
            "exact": 0,
            "low_complexity": 0,
            "pass": 0,
            "other": 0,
        }
    )
    for row in rows:
        for rrna_class in _row_classes(row):
            counts = by_class[rrna_class]
            counts["candidate"] += 1
            counts["silva"] += row["in_silva"] == "1"
            counts["rfam"] += row["in_rfam"] == "1"
            counts["exact"] += row["in_exact_difference"] == "1"
            counts["low_complexity"] += (
                row["status"] == "REJECT_LOW_COMPLEXITY"
            )
            counts["pass"] += row["status"] == "PASS"
            counts["other"] += row["in_other_cyclospora"] == "1"

    ordered_classes = list(_EXPECTED_CLASSES)
    ordered_classes.extend(sorted(set(by_class) - set(_EXPECTED_CLASSES)))
    return tuple(
        ClassSummary(
            rrna_class=rrna_class,
            candidate_count=by_class[rrna_class]["candidate"],
            silva_shared=by_class[rrna_class]["silva"],
            rfam_shared=by_class[rrna_class]["rfam"],
            exact_specific_count=by_class[rrna_class]["exact"],
            low_complexity_rejected=by_class[rrna_class]["low_complexity"],
            pass_count=by_class[rrna_class]["pass"],
            other_cyclospora_shared=by_class[rrna_class]["other"],
        )
        for rrna_class in ordered_classes
    )


def _locus_summaries(
    records: list[SeqRecord],
    occurrences: dict[str, list[tuple[int, str]]],
    rows: list[dict[str, str]],
    passing: set[str],
    loci_tsv: Path,
) -> tuple[LocusSummary, ...]:
    candidate_kmers = {row["kmer"] for row in rows}
    classes = _locus_classes(records, loci_tsv)
    summaries: list[LocusSummary] = []

    for record in records:
        sequence_length = len(normalize_dna(str(record.seq)))
        record_occurrences = occurrences[record.id]
        record_kmers = {kmer for _, kmer in record_occurrences}
        candidate_count = len(record_kmers & candidate_kmers)
        pass_count = len(record_kmers & passing)

        covered_bases: set[int] = set()
        covered_read_starts: set[int] = set()
        read_starts_total = max(0, sequence_length - READ_LENGTH + 1)
        for start, kmer in record_occurrences:
            if kmer not in passing:
                continue
            covered_bases.update(range(start, start + K))
            first_read_start = max(0, start + K - READ_LENGTH)
            last_read_start = min(start, read_starts_total - 1)
            if first_read_start <= last_read_start:
                covered_read_starts.update(
                    range(first_read_start, last_read_start + 1)
                )

        summaries.append(
            LocusSummary(
                target_record_id=record.id,
                rrna_class=classes[record.id],
                length=sequence_length,
                candidate_count=candidate_count,
                pass_count=pass_count,
                bases_covered=len(covered_bases),
                read_starts_total=read_starts_total,
                read_starts_covered=len(covered_read_starts),
            )
        )
    return tuple(sorted(summaries, key=lambda row: row.target_record_id))


def _artifact_provenance(
    target_fasta: Path,
    manifest: Path,
    baits: Path,
    exact_difference: Path,
    silva_intersection: Path,
    rfam_intersection: Path,
    other_cyclospora_intersection: Path,
    genus_compatible: Path,
    provenance: Path,
    loci_tsv: Path,
) -> tuple[tuple[str, Path], ...]:
    artifacts: list[tuple[str, Path]] = [
        ("curated target FASTA", Path(target_fasta)),
        ("final manifest", Path(manifest)),
        ("final bait FASTA", Path(baits)),
        ("Meryl exact difference", Path(exact_difference)),
        ("Meryl target-SILVA intersection", Path(silva_intersection)),
        ("Meryl target-Rfam intersection", Path(rfam_intersection)),
        (
            "Meryl target-other-Cyclospora intersection",
            Path(other_cyclospora_intersection),
        ),
        ("Meryl genus-compatible difference", Path(genus_compatible)),
        ("source provenance manifest", Path(provenance)),
    ]
    artifacts.append(("locus metadata", Path(loci_tsv)))
    return tuple(artifacts)


def _summary_rows(summary: VerificationSummary) -> list[tuple[str, str]]:
    rows = [
        ("conclusion", summary.conclusion),
        ("build_status", "PASS" if summary.bait_count else "INFEASIBLE"),
        ("index.k", str(summary.k)),
        ("index.w", str(summary.w)),
        ("index.bait_count", str(summary.bait_count)),
        ("index.indexed_minimizer_count", str(summary.indexed_count)),
        ("attrition.target_candidate_count", str(summary.candidate_count)),
        ("attrition.silva_shared_count", str(summary.silva_shared)),
        ("attrition.rfam_shared_count", str(summary.rfam_shared)),
        (
            "attrition.exact_specific_candidate_count",
            str(summary.exact_specific_count),
        ),
        (
            "attrition.low_complexity_rejected_count",
            str(summary.low_complexity_rejected),
        ),
        ("attrition.pass_count", str(summary.bait_count)),
        (
            "specificity.other_cyclospora_shared_count",
            str(summary.other_cyclospora_shared),
        ),
        (
            "genus.genus_compatible_candidate_count",
            str(summary.genus_compatible_count),
        ),
        ("genus.genus_level_status", summary.genus_level_status),
    ]
    for item in summary.classes:
        prefix = f"class.{item.rrna_class}"
        rows.extend(
            [
                (f"{prefix}.target_candidate_count", str(item.candidate_count)),
                (f"{prefix}.silva_shared_count", str(item.silva_shared)),
                (f"{prefix}.rfam_shared_count", str(item.rfam_shared)),
                (
                    f"{prefix}.exact_specific_candidate_count",
                    str(item.exact_specific_count),
                ),
                (
                    f"{prefix}.low_complexity_rejected_count",
                    str(item.low_complexity_rejected),
                ),
                (f"{prefix}.pass_count", str(item.pass_count)),
                (
                    f"{prefix}.other_cyclospora_shared_count",
                    str(item.other_cyclospora_shared),
                ),
            ]
        )
    for item in summary.loci:
        prefix = f"locus.{item.target_record_id}"
        rows.extend(
            [
                (f"{prefix}.rrna_class", item.rrna_class),
                (f"{prefix}.length", str(item.length)),
                (f"{prefix}.candidate_count", str(item.candidate_count)),
                (f"{prefix}.pass_count", str(item.pass_count)),
                (f"{prefix}.bases_covered", str(item.bases_covered)),
                (
                    f"{prefix}.read_starts_150_total",
                    str(item.read_starts_total),
                ),
                (
                    f"{prefix}.read_starts_150_covered",
                    str(item.read_starts_covered),
                ),
            ]
        )
    return rows


def write_summary_tsv(path: Path, summary: VerificationSummary) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("key", "value"))
        writer.writerows(_summary_rows(summary))


def _markdown_cell(value: object) -> str:
    return str(value).replace("|", r"\|").replace("\n", " ")


def _source_provenance_rows(
    path: Path, allow_empty: bool
) -> list[dict[str, str]]:
    path = Path(path)
    if not path.is_file():
        raise ValueError(f"provenance manifest is absent: {path}")
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        rows = list(reader)
    missing = set(_PROVENANCE_FIELDS) - fields
    if missing:
        raise ValueError(
            "provenance manifest is missing required columns: "
            + ", ".join(sorted(missing))
        )
    if not rows and not allow_empty:
        raise ValueError("production provenance manifest must contain rows")
    for row_number, row in enumerate(rows, start=2):
        for field in _PROVENANCE_FIELDS:
            if not row.get(field):
                raise ValueError(
                    f"provenance row {row_number} has an empty {field}"
                )
        retrieved = row["retrieved_at_utc"]
        try:
            parsed = datetime.fromisoformat(retrieved.replace("Z", "+00:00"))
        except ValueError as error:
            raise ValueError(
                f"provenance row {row_number} has malformed "
                f"retrieved_at_utc: {retrieved!r}"
            ) from error
        if parsed.tzinfo is None or parsed.utcoffset() is None:
            raise ValueError(
                f"provenance row {row_number} retrieved_at_utc "
                "must include a timezone"
            )
        if parsed.utcoffset() != timedelta(0):
            raise ValueError(
                f"provenance row {row_number} retrieved_at_utc "
                "must have zero UTC offset"
            )
    return sorted(
        rows,
        key=lambda row: (
            row.get("name", ""),
            row.get("release", ""),
            row.get("url", ""),
        ),
    )


def write_markdown_report(
    path: Path,
    summary: VerificationSummary,
    artifacts: tuple[tuple[str, Path], ...],
    source_rows: list[dict[str, str]],
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# *Cyclospora cayetanensis* mature-rRNA bait index verification",
        "",
        "## Conclusion",
        "",
        summary.conclusion,
        "",
        "## Provenance",
        "",
        "| artifact | file | bytes | SHA-256 |",
        "|---|---|---:|---|",
    ]
    for label, artifact in artifacts:
        lines.append(
            "| "
            + " | ".join(
                (
                    _markdown_cell(label),
                    _markdown_cell(artifact.name),
                    str(artifact.stat().st_size),
                    sha256_file(artifact),
                )
            )
            + " |"
        )

    if source_rows:
        lines.extend(
            [
                "",
                "### Pinned input sources",
                "",
                "| source | release | URL | retrieved UTC | bytes | SHA-256 |",
                "|---|---|---|---|---:|---|",
            ]
        )
        for row in source_rows:
            lines.append(
                "| "
                + " | ".join(
                    _markdown_cell(row.get(field, ""))
                    for field in (
                        "name",
                        "release",
                        "url",
                        "retrieved_at_utc",
                        "bytes",
                        "sha256",
                    )
                )
                + " |"
            )

    lines.extend(
        [
            "",
            "## Stage attrition",
            "",
            "| rRNA class | target candidates | shared with SILVA | "
            "shared with Rfam | shared with other *Cyclospora* | "
            "exact species-specific | low-complexity rejected | passing baits |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for item in summary.classes:
        lines.append(
            f"| {_markdown_cell(item.rrna_class)} | {item.candidate_count} | "
            f"{item.silva_shared} | {item.rfam_shared} | "
            f"{item.other_cyclospora_shared} | {item.exact_specific_count} | "
            f"{item.low_complexity_rejected} | {item.pass_count} |"
        )
    lines.append(
        f"| **all** | **{summary.candidate_count}** | "
        f"**{summary.silva_shared}** | **{summary.rfam_shared}** | "
        f"**{summary.other_cyclospora_shared}** | "
        f"**{summary.exact_specific_count}** | "
        f"**{summary.low_complexity_rejected}** | "
        f"**{summary.bait_count}** |"
    )

    lines.extend(
        [
            "",
            "## Per-locus coverage",
            "",
            "Coverage is descriptive. A covered 150 bp read start is a valid "
            "full-length start whose read contains an entire bait occurrence.",
            "",
            "| target locus | class | bases | candidates | passing baits | "
            "bases covered | 150 bp read starts covered |",
            "|---|---|---:|---:|---:|---:|---:|",
        ]
    )
    for item in summary.loci:
        lines.append(
            f"| {_markdown_cell(item.target_record_id)} | "
            f"{_markdown_cell(item.rrna_class)} | {item.length} | "
            f"{item.candidate_count} | {item.pass_count} | "
            f"{item.bases_covered} | {item.read_starts_covered}/"
            f"{item.read_starts_total} |"
        )

    lines.extend(
        [
            "",
            "## Other-*Cyclospora* sharing",
            "",
            f"{summary.other_cyclospora_shared} target candidate(s) were shared "
            "with the available non-*cayetanensis Cyclospora* references and "
            "were excluded from the species-specific bait set.",
            "",
            "Subtraction is intentionally conservative: complete response "
            "records are used, so any ITS or other non-rRNA sequence carried "
            "in a mixed record is also subtracted.",
            "",
            "| rRNA class | shared candidates |",
            "|---|---:|",
        ]
    )
    for item in summary.classes:
        lines.append(
            f"| {_markdown_cell(item.rrna_class)} | "
            f"{item.other_cyclospora_shared} |"
        )

    lines.extend(
        [
            "",
            "## Genus-level feasibility",
            "",
            f"Status: `{summary.genus_level_status}`.",
            "",
            f"{summary.genus_compatible_count} target candidate(s) remain "
            "after SILVA and Rfam subtraction before other-*Cyclospora* and "
            "entropy filtering.",
        ]
    )

    lines.extend(
        [
            "",
            "## Final index properties",
            "",
            "| property | value |",
            "|---|---:|",
            f"| index emitted | {'yes' if summary.bait_count else 'no'} |",
            f"| k | {summary.k} |",
            f"| w | {summary.w} |",
            f"| bait records | {summary.bait_count} |",
            f"| distinct indexed minimizers | {summary.indexed_count} |",
            "",
            "## Scope limitation",
            "",
            "Wastewater performance has not been evaluated. Static correctness "
            "does not establish wastewater sensitivity, specificity, or a "
            "validated sample-calling threshold.",
            "",
        ]
    )
    path.write_text("\n".join(lines))


def verify_static(
    target_fasta: Path,
    manifest: Path,
    baits: Path,
    index_info: str | None,
    *,
    exact_difference: Path,
    silva_intersection: Path,
    rfam_intersection: Path,
    other_cyclospora_intersection: Path,
    genus_compatible: Path,
    loci_tsv: Path,
    provenance: Path,
    core_nt_decisions: Path | None = None,
    allow_empty_provenance: bool = False,
    summary_path: Path | None = None,
    report_path: Path | None = None,
) -> VerificationSummary:
    records = _read_target_records(Path(target_fasta))
    (
        target_kmers,
        occurrences,
        classes_by_kmer,
        occurrences_by_kmer,
    ) = _target_kmers(records)
    bait_records, bait_set = _read_baits(Path(baits))
    rows = _read_manifest(Path(manifest))
    source_rows = _source_provenance_rows(
        Path(provenance), allow_empty_provenance
    )
    printed_exact = _read_meryl_print(Path(exact_difference))
    printed_silva = _read_meryl_print(Path(silva_intersection))
    printed_rfam = _read_meryl_print(Path(rfam_intersection))
    printed_other = _read_meryl_print(Path(other_cyclospora_intersection))
    printed_genus = _read_meryl_print(Path(genus_compatible))

    absent = sorted(bait_set - target_kmers)
    if absent:
        raise ValueError(
            "bait is absent from every eligible target rRNA: " + absent[0]
        )

    passing = {row["kmer"] for row in rows if row["status"] == "PASS"}
    if passing != bait_set:
        missing = sorted(passing - bait_set)
        rejected = sorted(bait_set - passing)
        detail = []
        if missing:
            detail.append(f"missing bait for PASS row {missing[0]}")
        if rejected:
            detail.append(f"bait without status=PASS {rejected[0]}")
        raise ValueError(
            "status=PASS manifest k-mers must equal bait sequences exactly"
            + (": " + "; ".join(detail) if detail else "")
        )
    _validate_core_nt_annotations(rows, core_nt_decisions)

    manifest_kmers = {row["kmer"] for row in rows}
    if manifest_kmers != target_kmers:
        missing = sorted(target_kmers - manifest_kmers)
        extra = sorted(manifest_kmers - target_kmers)
        detail = []
        if missing:
            detail.append(f"missing {missing[0]}")
        if extra:
            detail.append(f"extra {extra[0]}")
        raise ValueError(
            "manifest k-mers must equal the complete target canonical "
            "31-mer set"
            + (": " + "; ".join(detail) if detail else "")
        )

    _validate_manifest_annotations(
        rows=rows,
        target_kmers=target_kmers,
        classes_by_kmer=classes_by_kmer,
        occurrences_by_kmer=occurrences_by_kmer,
        silva_intersection=printed_silva,
        rfam_intersection=printed_rfam,
        other_cyclospora_intersection=printed_other,
        exact_difference=printed_exact,
        genus_compatible=printed_genus,
        passing=passing,
    )

    non_exact_passing = sorted(passing - printed_exact)
    if non_exact_passing:
        raise ValueError(
            "every PASS bait must belong to the exact-difference set: "
            + non_exact_passing[0]
        )

    if bait_records:
        if not index_info:
            raise ValueError("a nonempty bait set requires Deacon index info")
        k, w, indexed_count = parse_index_info(index_info)
        if k != K:
            raise ValueError(f"Deacon index k must be {K}; found {k}")
        if w != W:
            raise ValueError(f"Deacon index w must be {W}; found {w}")
        if indexed_count != len(bait_records):
            raise ValueError(
                "Deacon indexed minimizer count differs from bait count: "
                f"{indexed_count} != {len(bait_records)}"
            )
    else:
        k, w, indexed_count = K, W, 0

    summary = VerificationSummary(
        bait_count=len(bait_records),
        k=k,
        w=w,
        indexed_count=indexed_count,
        candidate_count=len(rows),
        silva_shared=len(printed_silva),
        rfam_shared=len(printed_rfam),
        exact_specific_count=len(printed_exact),
        low_complexity_rejected=sum(
            row["status"] == "REJECT_LOW_COMPLEXITY" for row in rows
        ),
        other_cyclospora_shared=len(printed_other),
        genus_compatible_count=len(printed_genus),
        classes=_class_summaries(rows),
        loci=_locus_summaries(
            records, occurrences, rows, passing, loci_tsv
        ),
    )

    if summary_path is not None:
        write_summary_tsv(Path(summary_path), summary)
    if report_path is not None:
        write_markdown_report(
            Path(report_path),
            summary,
            _artifact_provenance(
                Path(target_fasta),
                Path(manifest),
                Path(baits),
                Path(exact_difference),
                Path(silva_intersection),
                Path(rfam_intersection),
                Path(other_cyclospora_intersection),
                Path(genus_compatible),
                Path(provenance),
                loci_tsv,
            ),
            source_rows,
        )
    return summary


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Independently verify and report a packaged rRNA bait index."
    )
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--baits", type=Path, required=True)
    parser.add_argument("--index", type=Path)
    parser.add_argument("--exact-difference", type=Path, required=True)
    parser.add_argument("--silva-intersection", type=Path, required=True)
    parser.add_argument("--rfam-intersection", type=Path, required=True)
    parser.add_argument(
        "--other-cyclospora-intersection", type=Path, required=True
    )
    parser.add_argument("--genus-compatible", type=Path, required=True)
    parser.add_argument("--loci", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    parser.add_argument("--core-nt-decisions", type=Path)
    parser.add_argument("--allow-empty-provenance", action="store_true")
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    args = _parser().parse_args(sys.argv[1:] if argv is None else argv)
    try:
        index_info = read_index_info(args.index) if args.index else None
        verify_static(
            target_fasta=args.target,
            manifest=args.manifest,
            baits=args.baits,
            index_info=index_info,
            exact_difference=args.exact_difference,
            silva_intersection=args.silva_intersection,
            rfam_intersection=args.rfam_intersection,
            other_cyclospora_intersection=(
                args.other_cyclospora_intersection
            ),
            genus_compatible=args.genus_compatible,
            loci_tsv=args.loci,
            provenance=args.provenance,
            core_nt_decisions=args.core_nt_decisions,
            allow_empty_provenance=args.allow_empty_provenance,
            summary_path=args.summary,
            report_path=args.report,
        )
    except ValueError as error:
        raise SystemExit(f"verification failed: {error}") from error


if __name__ == "__main__":
    main()
