"""Classify core-nt BLAST hits for candidate 31-mer baits."""

from __future__ import annotations

import argparse
import csv
import sys
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

from rrna_bait.core import normalize_dna

_BLAST_FIELD_COUNT = 13
_EXACT_NON_TARGET = "REJECT_CORE_NT_EXACT_NON_TARGET"
_PASS = "PASS_CORE_NT"
_BLAST_FIELDS = (
    "qseqid",
    "saccver",
    "staxids",
    "sscinames",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "stitle",
)


@dataclass(frozen=True)
class BlastHit:
    """One row from the fixed core-nt BLAST outfmt."""

    query_id: str
    accession: str
    taxids: str
    scientific_names: str
    percent_identity: float
    alignment_length: int
    mismatches: int
    gap_opens: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    title: str

    @classmethod
    def from_tsv(cls, line: str) -> BlastHit:
        """Parse one tab-delimited BLAST row in the required field order."""
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) != _BLAST_FIELD_COUNT:
            raise ValueError(
                f"expected {_BLAST_FIELD_COUNT} BLAST fields, found {len(fields)}"
            )
        try:
            return cls(
                query_id=fields[0],
                accession=fields[1],
                taxids=fields[2],
                scientific_names=fields[3],
                percent_identity=float(fields[4]),
                alignment_length=int(fields[5]),
                mismatches=int(fields[6]),
                gap_opens=int(fields[7]),
                query_start=int(fields[8]),
                query_end=int(fields[9]),
                subject_start=int(fields[10]),
                subject_end=int(fields[11]),
                title=fields[12],
            )
        except ValueError as error:
            raise ValueError("BLAST row contains an invalid numeric field") from error

    @property
    def taxid_values(self) -> frozenset[str]:
        """Return the reported subject taxids, preserving unknown values as non-targets."""
        return frozenset(value.strip() for value in self.taxids.split(";") if value.strip())

    @property
    def is_exact_candidate(self) -> bool:
        """Whether this HSP is an exact, full-length 31-mer match."""
        return (
            self.percent_identity == 100.0
            and self.alignment_length == 31
            and self.mismatches == 0
            and self.gap_opens == 0
            and self.query_start == 1
            and self.query_end == 31
        )

    @property
    def is_ungapped_30_of_31(self) -> bool:
        """Whether this HSP has one mismatch across the full 31-mer bait."""
        return (
            self.alignment_length == 31
            and self.mismatches == 1
            and self.gap_opens == 0
            and self.query_start == 1
            and self.query_end == 31
        )


@dataclass(frozen=True)
class BaitDecision:
    """Core-nt outcome and supporting hit counts for one bait."""

    bait_id: str
    status: str
    exact_target_count: int
    exact_non_target_count: int
    near_match_count: int


def classify_hits(
    bait_ids: set[str], hits: Iterable[BlastHit], target_taxid: str = "88456"
) -> dict[str, BaitDecision]:
    """Classify baits using exact matches only; retain 30/31 HSPs as context."""
    exact_target_counts = {bait_id: 0 for bait_id in bait_ids}
    exact_non_target_counts = {bait_id: 0 for bait_id in bait_ids}
    near_match_counts = {bait_id: 0 for bait_id in bait_ids}

    for hit in hits:
        if hit.query_id not in bait_ids:
            continue
        if hit.is_ungapped_30_of_31:
            near_match_counts[hit.query_id] += 1
        if not hit.is_exact_candidate:
            continue
        if hit.taxid_values == {target_taxid}:
            exact_target_counts[hit.query_id] += 1
        else:
            exact_non_target_counts[hit.query_id] += 1

    return {
        bait_id: BaitDecision(
            bait_id=bait_id,
            status=(
                _EXACT_NON_TARGET
                if exact_non_target_counts[bait_id]
                else _PASS
            ),
            exact_target_count=exact_target_counts[bait_id],
            exact_non_target_count=exact_non_target_counts[bait_id],
            near_match_count=near_match_counts[bait_id],
        )
        for bait_id in bait_ids
    }


def write_decisions(decisions: Mapping[str, BaitDecision], path: Path) -> None:
    """Write decisions as stable tab-delimited records, ordered by bait ID."""
    lines = [
        "bait_id\tstatus\texact_target_count\texact_non_target_count\tnear_match_count"
    ]
    lines.extend(
        "\t".join(
            (
                decision.bait_id,
                decision.status,
                str(decision.exact_target_count),
                str(decision.exact_non_target_count),
                str(decision.near_match_count),
            )
        )
        for _, decision in sorted(decisions.items())
    )
    path.write_text("\n".join(lines) + "\n")


def _read_hits(path: Path) -> list[BlastHit]:
    hits: list[BlastHit] = []
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                hits.append(BlastHit.from_tsv(line))
            except ValueError as error:
                raise ValueError(
                    f"invalid BLAST row {line_number}: {error}"
                ) from error
    return hits


def finalize_validation(
    baits: Path,
    blast_tsv: Path,
    manifest_in: Path,
    manifest_out: Path,
    baits_out: Path,
    decisions_out: Path,
    near_hits_out: Path,
) -> tuple[int, int]:
    """Apply exact-only core-nt decisions while preserving audit rows."""
    bait_records = list(SeqIO.parse(Path(baits), "fasta"))
    if not bait_records:
        raise ValueError(f"input bait FASTA has no records: {baits}")
    bait_ids = [record.id for record in bait_records]
    if len(bait_ids) != len(set(bait_ids)):
        raise ValueError("input bait FASTA contains duplicate identifiers")
    sequences = [normalize_dna(str(record.seq)) for record in bait_records]
    if any(len(sequence) != 31 for sequence in sequences):
        raise ValueError("every input bait must be exactly 31 bases")
    if len(sequences) != len(set(sequences)):
        raise ValueError("input bait FASTA contains duplicate sequences")

    hits = _read_hits(Path(blast_tsv))
    unknown_queries = sorted({hit.query_id for hit in hits} - set(bait_ids))
    if unknown_queries:
        raise ValueError(
            "BLAST TSV contains an unknown bait identifier: "
            + unknown_queries[0]
        )
    decisions = classify_hits(set(bait_ids), hits)

    with Path(manifest_in).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or ())
        rows = list(reader)
    required = {"kmer", "status", "rejection_reason"}
    missing = required - set(fieldnames)
    if missing:
        raise ValueError(
            "manifest is missing required columns: " + ", ".join(sorted(missing))
        )
    if len({row["kmer"] for row in rows}) != len(rows):
        raise ValueError("manifest contains duplicate k-mer rows")
    by_kmer = {row["kmer"]: row for row in rows}
    missing_baits = [
        sequence for sequence in sequences if sequence not in by_kmer
    ]
    if missing_baits:
        raise ValueError("input bait is absent from manifest: " + missing_baits[0])

    screened = set(sequences)
    for field in (
        "core_nt_bait_id",
        "core_nt_status",
        "core_nt_rejection_reason",
    ):
        if field not in fieldnames:
            fieldnames.append(field)
    for row in rows:
        if row["kmer"] not in screened:
            row["core_nt_bait_id"] = ""
            row["core_nt_status"] = "NOT_APPLICABLE"
            row["core_nt_rejection_reason"] = "not_screened_static_reject"

    retained_records = []
    for record, sequence in zip(bait_records, sequences, strict=True):
        row = by_kmer[sequence]
        if row["status"] != "PASS":
            raise ValueError(
                f"input bait manifest row is not PASS before core-nt: {record.id}"
            )
        decision = decisions[record.id]
        row["core_nt_bait_id"] = record.id
        row["core_nt_status"] = decision.status
        if decision.status == _EXACT_NON_TARGET:
            row["core_nt_rejection_reason"] = "exact_non_target_31_of_31"
            row["status"] = _EXACT_NON_TARGET
            row["rejection_reason"] = "core_nt_exact_non_target"
        else:
            row["core_nt_rejection_reason"] = "none"
            retained_records.append(record)

    with Path(manifest_out).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    SeqIO.write(retained_records, Path(baits_out), "fasta")

    ordered_decisions = {bait_id: decisions[bait_id] for bait_id in bait_ids}
    lines = [
        "bait_id\tstatus\texact_target_count\texact_non_target_count"
        "\tnear_match_count"
    ]
    lines.extend(
        "\t".join(
            (
                item.bait_id,
                item.status,
                str(item.exact_target_count),
                str(item.exact_non_target_count),
                str(item.near_match_count),
            )
        )
        for item in ordered_decisions.values()
    )
    Path(decisions_out).write_text("\n".join(lines) + "\n")

    with Path(near_hits_out).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(_BLAST_FIELDS)
        for hit in hits:
            if hit.is_ungapped_30_of_31:
                writer.writerow(
                    (
                        hit.query_id,
                        hit.accession,
                        hit.taxids,
                        hit.scientific_names,
                        f"{hit.percent_identity:g}",
                        hit.alignment_length,
                        hit.mismatches,
                        hit.gap_opens,
                        hit.query_start,
                        hit.query_end,
                        hit.subject_start,
                        hit.subject_end,
                        hit.title,
                    )
                )
    return len(bait_records), len(retained_records)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Finalize bait decisions from an all-bait core-nt BLAST TSV."
    )
    parser.add_argument("--baits", type=Path, required=True)
    parser.add_argument("--blast-tsv", type=Path, required=True)
    parser.add_argument("--manifest-in", type=Path, required=True)
    parser.add_argument("--manifest-out", type=Path, required=True)
    parser.add_argument("--baits-out", type=Path, required=True)
    parser.add_argument("--decisions-out", type=Path, required=True)
    parser.add_argument("--near-hits-out", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = _parser().parse_args(sys.argv[1:] if argv is None else argv)
    try:
        total, retained = finalize_validation(
            baits=args.baits,
            blast_tsv=args.blast_tsv,
            manifest_in=args.manifest_in,
            manifest_out=args.manifest_out,
            baits_out=args.baits_out,
            decisions_out=args.decisions_out,
            near_hits_out=args.near_hits_out,
        )
    except ValueError as error:
        raise SystemExit(f"core-nt finalization failed: {error}") from error
    print(f"Core-nt decisions: {total}; retained baits: {retained}")


if __name__ == "__main__":
    main()
