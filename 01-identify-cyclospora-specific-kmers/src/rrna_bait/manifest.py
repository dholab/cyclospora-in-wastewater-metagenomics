from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO

from rrna_bait.core import canonical_kmer, iter_canonical_kmers, normalize_dna


_MANIFEST_FIELDS = [
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
]
_REJECTION_SOURCES = (
    ("other_cyclospora", "in_other_cyclospora"),
    ("silva", "in_silva"),
    ("rfam", "in_rfam"),
)


@dataclass
class KmerMetadata:
    rrna_classes: set[str] = field(default_factory=set)
    target_records: set[str] = field(default_factory=set)
    zero_based_starts: set[int] = field(default_factory=set)
    occurrences: set[tuple[str, int]] = field(default_factory=set)
    copy_count: int = 0


def _locus_metadata(
    loci_tsv: Path,
) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    support_by_target: dict[str, set[str]] = defaultdict(set)
    classes_by_target: dict[str, set[str]] = defaultdict(set)
    with Path(loci_tsv).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            target_record_id = row.get("target_record_id", "")
            if not target_record_id:
                continue
            support = (
                f"{row['subject_accession']}:{row['start']}-{row['end']}:"
                f"{row['strand']}"
            )
            support_by_target[target_record_id].add(support)
            if row.get("rrna_class"):
                classes_by_target[target_record_id].add(row["rrna_class"])
    return dict(support_by_target), dict(classes_by_target)


def target_kmer_metadata(
    target_fasta: Path, loci_tsv: Path, k: int = 31
) -> dict[str, KmerMetadata]:
    _, classes_by_target = _locus_metadata(Path(loci_tsv))
    metadata: dict[str, KmerMetadata] = {}
    for record in SeqIO.parse(Path(target_fasta), "fasta"):
        rrna_classes = classes_by_target.get(
            record.id, {record.id.split("|", 1)[0]}
        )
        for start, kmer in iter_canonical_kmers(str(record.seq), k):
            annotation = metadata.setdefault(kmer, KmerMetadata())
            annotation.rrna_classes.update(rrna_classes)
            annotation.target_records.add(record.id)
            annotation.zero_based_starts.add(start)
            annotation.occurrences.add((record.id, start))

    for annotation in metadata.values():
        annotation.copy_count = len(annotation.occurrences)
    return {kmer: metadata[kmer] for kmer in sorted(metadata)}


def read_meryl_print(path: Path) -> set[str]:
    kmers: set[str] = set()
    with Path(path).open() as handle:
        for line in handle:
            fields = line.split()
            if fields:
                kmers.add(canonical_kmer(fields[0]))
    return kmers


def _write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=_MANIFEST_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _write_raw_baits(
    path: Path, candidates: list[tuple[str, KmerMetadata]]
) -> None:
    with Path(path).open("w") as handle:
        for index, (kmer, metadata) in enumerate(candidates, start=1):
            rrna_classes = ",".join(sorted(metadata.rrna_classes))
            handle.write(f">raw_kmer_{index:06d}|{rrna_classes}\n{kmer}\n")


def write_raw_manifest(
    target_fasta: Path,
    loci_tsv: Path,
    silva_print: Path,
    rfam_print: Path,
    other_cyclospora_print: Path,
    exact_difference_print: Path,
    raw_manifest: Path,
    raw_baits: Path,
    k: int = 31,
) -> int:
    metadata = target_kmer_metadata(Path(target_fasta), Path(loci_tsv), k)
    backgrounds = {
        "silva": read_meryl_print(Path(silva_print)),
        "rfam": read_meryl_print(Path(rfam_print)),
        "other_cyclospora": read_meryl_print(Path(other_cyclospora_print)),
    }
    exact_difference = read_meryl_print(Path(exact_difference_print))
    candidates: list[tuple[str, KmerMetadata]] = []
    rows: list[dict[str, str]] = []

    for kmer, annotation in metadata.items():
        membership = {
            source: kmer in printed_kmers
            for source, printed_kmers in backgrounds.items()
        }
        in_exact_difference = kmer in exact_difference
        rejection_reason = ";".join(
            source
            for source, _ in _REJECTION_SOURCES
            if membership[source]
        )
        if in_exact_difference:
            candidates.append((kmer, annotation))
        rows.append(
            {
                "kmer": kmer,
                "rrna_classes": ",".join(sorted(annotation.rrna_classes)),
                "target_records": ",".join(sorted(annotation.target_records)),
                "target_starts": ",".join(
                    str(start) for start in sorted(annotation.zero_based_starts)
                ),
                "target_occurrences": ";".join(
                    f"{record_id}:{start}"
                    for record_id, start in sorted(annotation.occurrences)
                ),
                "target_copy_count": str(annotation.copy_count),
                "in_silva": str(int(membership["silva"])),
                "in_rfam": str(int(membership["rfam"])),
                "in_other_cyclospora": str(
                    int(membership["other_cyclospora"])
                ),
                "in_exact_difference": str(int(in_exact_difference)),
                "status": (
                    "CANDIDATE"
                    if in_exact_difference
                    else "REJECT_STATIC_BACKGROUND"
                ),
                "rejection_reason": (
                    "none" if in_exact_difference else rejection_reason
                ),
            }
        )

    _write_manifest(Path(raw_manifest), rows)
    _write_raw_baits(Path(raw_baits), candidates)
    return len(candidates)


def finalize_manifest(
    raw_manifest: Path,
    entropy_pass_fasta: Path,
    final_manifest: Path,
    final_baits: Path,
) -> int:
    entropy_pass = {
        canonical_kmer(normalize_dna(str(record.seq)))
        for record in SeqIO.parse(Path(entropy_pass_fasta), "fasta")
    }
    with Path(raw_manifest).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rows.sort(key=lambda row: row["kmer"])

    passing_rows: list[dict[str, str]] = []
    for row in rows:
        if row["in_exact_difference"] != "1":
            row["status"] = "REJECT_STATIC_BACKGROUND"
        elif row["kmer"] in entropy_pass:
            row["status"] = "PASS"
            row["rejection_reason"] = "none"
            passing_rows.append(row)
        else:
            row["status"] = "REJECT_LOW_COMPLEXITY"
            row["rejection_reason"] = "low_complexity"

    _write_manifest(Path(final_manifest), rows)
    with Path(final_baits).open("w") as handle:
        for index, row in enumerate(passing_rows, start=1):
            handle.write(
                f">cc_rrna_kmer_{index:06d}|{row['rrna_classes']}\n"
                f"{row['kmer']}\n"
            )
    return len(passing_rows)


def _raw_parser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("raw")
    parser.add_argument("--target-fasta", type=Path, required=True)
    parser.add_argument("--loci-tsv", type=Path, required=True)
    parser.add_argument("--silva-print", type=Path, required=True)
    parser.add_argument("--rfam-print", type=Path, required=True)
    parser.add_argument("--other-cyclospora-print", type=Path, required=True)
    parser.add_argument("--exact-difference-print", type=Path, required=True)
    parser.add_argument("--raw-manifest", type=Path, required=True)
    parser.add_argument("--raw-baits", type=Path, required=True)
    parser.add_argument("-k", type=int, default=31)


def _finalize_parser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("finalize")
    parser.add_argument("--raw-manifest", type=Path, required=True)
    parser.add_argument("--entropy-pass-fasta", type=Path, required=True)
    parser.add_argument("--final-manifest", type=Path, required=True)
    parser.add_argument("--final-baits", type=Path, required=True)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Annotate exact target k-mers and package passing baits."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    _raw_parser(subparsers)
    _finalize_parser(subparsers)
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    args = _parser().parse_args(sys.argv[1:] if argv is None else argv)
    if args.command == "raw":
        write_raw_manifest(
            args.target_fasta,
            args.loci_tsv,
            args.silva_print,
            args.rfam_print,
            args.other_cyclospora_print,
            args.exact_difference_print,
            args.raw_manifest,
            args.raw_baits,
            args.k,
        )
    else:
        finalize_manifest(
            args.raw_manifest,
            args.entropy_pass_fasta,
            args.final_manifest,
            args.final_baits,
        )


if __name__ == "__main__":
    main()
