from __future__ import annotations

import argparse
import csv
import sys
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from rrna_bait.core import normalize_dna, reverse_complement


@dataclass(frozen=True)
class TargetRule:
    accession: str
    rrna_class: str
    query_start: int
    query_end: int
    min_identity: float
    min_query_coverage: float


@dataclass(frozen=True)
class BlastHit:
    query: str
    query_length: int
    query_start: int
    query_end: int
    subject: str
    subject_start: int
    subject_end: int
    alignment_length: int
    identity: float

    @classmethod
    def from_row(cls, row: str) -> BlastHit:
        fields = row.rstrip().split("\t")
        return cls(
            fields[0],
            int(fields[1]),
            int(fields[2]),
            int(fields[3]),
            fields[4],
            int(fields[5]),
            int(fields[6]),
            int(fields[7]),
            float(fields[8]),
        )


@dataclass(frozen=True)
class _Locus:
    rrna_class: str
    query: str
    subject: str
    start: int
    end: int
    strand: str
    sequence: str


def accept_hit(hit: BlastHit, rule: TargetRule) -> bool:
    query_coverage = (
        abs(hit.query_end - hit.query_start) + 1
    ) / hit.query_length * 100
    return (
        hit.identity >= rule.min_identity
        and query_coverage >= rule.min_query_coverage
    )


def load_rules(path: Path) -> list[TargetRule]:
    with path.open(newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return [
            TargetRule(
                accession=row["accession"],
                rrna_class=row["rrna_class"],
                query_start=int(row["query_start"]),
                query_end=int(row["query_end"]),
                min_identity=float(row["min_identity"]),
                min_query_coverage=float(row["min_query_coverage"]),
            )
            for row in rows
        ]


def _records_by_accession(records: Iterable[SeqRecord]) -> dict[str, SeqRecord]:
    return {record.id: record for record in records}


def _rules_by_accession(
    rules: Mapping[str, TargetRule] | Iterable[TargetRule],
) -> dict[str, TargetRule]:
    if isinstance(rules, Mapping):
        return dict(rules)
    return {rule.accession: rule for rule in rules}


def _validate_query_records(
    query_records: Mapping[str, SeqRecord], rules: Mapping[str, TargetRule]
) -> None:
    for accession, query_record in query_records.items():
        rule = rules.get(accession)
        if rule is None:
            continue
        if rule.query_end == 0:
            continue
        expected_length = rule.query_end - rule.query_start + 1
        if len(query_record) != expected_length:
            raise ValueError(
                f"Supplied query length mismatch for {accession}: "
                f"expected {expected_length}, got {len(query_record)}"
            )


def _validate_hit_query_lengths(
    hits: Iterable[BlastHit],
    query_records: Mapping[str, SeqRecord],
    rules: Mapping[str, TargetRule],
) -> list[BlastHit]:
    configured_hits: list[BlastHit] = []
    for hit in hits:
        query_record = query_records.get(hit.query)
        if query_record is None or hit.query not in rules:
            continue
        if hit.query_length != len(query_record):
            raise ValueError(
                f"BLAST query length mismatch for {hit.query}: "
                f"expected {len(query_record)}, got {hit.query_length}"
            )
        configured_hits.append(hit)
    return configured_hits


def prepare_queries(
    downloaded_fasta: Path,
    rules: Iterable[TargetRule],
    output_fasta: Path,
) -> None:
    source_records = _records_by_accession(SeqIO.parse(downloaded_fasta, "fasta"))
    prepared: list[SeqRecord] = []
    for rule in rules:
        try:
            source = source_records[rule.accession]
        except KeyError as error:
            raise ValueError(f"Missing query accession: {rule.accession}") from error
        end = len(source) if rule.query_end == 0 else rule.query_end
        if rule.query_start < 1 or end < rule.query_start or end > len(source):
            raise ValueError(
                f"Query interval out of range for {rule.accession}: "
                f"{rule.query_start}-{end} (length {len(source)})"
            )
        sequence = normalize_dna(str(source.seq)[rule.query_start - 1 : end])
        prepared.append(SeqRecord(Seq(sequence), id=rule.accession, description=""))
    SeqIO.write(prepared, output_fasta, "fasta")


def _accepted_loci(
    subjects: Mapping[str, SeqRecord] | Iterable[SeqRecord],
    hits: Iterable[BlastHit],
    rules: Mapping[str, TargetRule] | Iterable[TargetRule],
) -> list[_Locus]:
    subject_records = (
        dict(subjects)
        if isinstance(subjects, Mapping)
        else _records_by_accession(subjects)
    )
    rule_by_accession = _rules_by_accession(rules)
    loci: list[_Locus] = []
    for hit in hits:
        rule = rule_by_accession.get(hit.query)
        subject = subject_records.get(hit.subject)
        if rule is None or subject is None or not accept_hit(hit, rule):
            continue
        start, end = sorted((hit.subject_start, hit.subject_end))
        if start < 1 or end > len(subject):
            continue
        sequence = normalize_dna(str(subject.seq)[start - 1 : end])
        strand = "+" if hit.subject_start <= hit.subject_end else "-"
        if strand == "-":
            sequence = reverse_complement(sequence)
        if not set(sequence) <= {"A", "C", "G", "T"}:
            continue
        loci.append(
            _Locus(rule.rrna_class, hit.query, hit.subject, start, end, strand, sequence)
        )
    return sorted(
        loci,
        key=lambda locus: (
            locus.rrna_class,
            locus.sequence,
            locus.query,
            locus.subject,
            locus.start,
            locus.end,
            locus.strand,
        ),
    )


def extract_loci(
    subjects: Mapping[str, SeqRecord] | Iterable[SeqRecord],
    hits: Iterable[BlastHit],
    rules: Mapping[str, TargetRule] | Iterable[TargetRule],
) -> list[SeqRecord]:
    records: list[SeqRecord] = []
    seen_sequences: set[tuple[str, str]] = set()
    for locus in _accepted_loci(subjects, hits, rules):
        sequence_key = (locus.rrna_class, locus.sequence)
        if sequence_key in seen_sequences:
            continue
        seen_sequences.add(sequence_key)
        record_id = (
            f"{locus.rrna_class}|{locus.query}|{locus.subject}|"
            f"{locus.start}-{locus.end}|{locus.strand}"
        )
        records.append(SeqRecord(Seq(locus.sequence), id=record_id, description=""))
    return records


def _read_hits(path: Path) -> list[BlastHit]:
    with path.open() as handle:
        return [
            BlastHit.from_row(row)
            for row in handle
            if row.strip() and not row.startswith("#")
        ]


def curate_targets(
    *,
    queries_path: Path,
    subject_path: Path,
    rules_path: Path,
    hits_path: Path,
    fasta_out: Path,
    loci_out: Path,
    expected_classes: set[str] | None,
) -> None:
    rules = load_rules(rules_path)
    rule_by_accession = _rules_by_accession(rules)
    query_records = _records_by_accession(SeqIO.parse(queries_path, "fasta"))
    subject_records = _records_by_accession(SeqIO.parse(subject_path, "fasta"))
    hits = _read_hits(hits_path)
    _validate_query_records(query_records, rule_by_accession)
    configured_hits = _validate_hit_query_lengths(hits, query_records, rule_by_accession)
    loci = _accepted_loci(subject_records, configured_hits, rule_by_accession)
    records = extract_loci(subject_records, configured_hits, rule_by_accession)

    observed_classes = {record.id.split("|", 1)[0] for record in records}
    if expected_classes is not None:
        missing = sorted(expected_classes - observed_classes)
        if missing:
            raise ValueError(f"Missing expected rRNA classes: {', '.join(missing)}")

    SeqIO.write(records, fasta_out, "fasta")
    retained_sequences: set[tuple[str, str]] = {
        (record.id.split("|", 1)[0], str(record.seq)) for record in records
    }
    record_ids_by_sequence = {
        (record.id.split("|", 1)[0], str(record.seq)): record.id for record in records
    }
    with loci_out.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "target_record_id",
                "rrna_class",
                "query_accession",
                "subject_accession",
                "start",
                "end",
                "strand",
            ]
        )
        for locus in loci:
            sequence_key = (locus.rrna_class, locus.sequence)
            if sequence_key in retained_sequences:
                writer.writerow(
                    [
                        record_ids_by_sequence[sequence_key],
                        locus.rrna_class,
                        locus.query,
                        locus.subject,
                        locus.start,
                        locus.end,
                        locus.strand,
                    ]
                )


def _curate_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Curate RefSeq-backed mature rRNA targets.")
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--subject", type=Path, required=True)
    parser.add_argument("--rules", type=Path, required=True)
    parser.add_argument("--hits", type=Path, required=True)
    parser.add_argument("--fasta-out", type=Path, required=True)
    parser.add_argument("--loci-out", type=Path, required=True)
    return parser


def _prepare_queries_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare mature rRNA query FASTA records.")
    parser.add_argument("--downloaded-fasta", "--queries", dest="downloaded_fasta", type=Path, required=True)
    parser.add_argument("--rules", type=Path, required=True)
    parser.add_argument("--output-fasta", "--fasta-out", dest="output_fasta", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    arguments = list(sys.argv[1:] if argv is None else argv)
    if arguments and arguments[0] == "prepare-queries":
        args = _prepare_queries_parser().parse_args(arguments[1:])
        prepare_queries(args.downloaded_fasta, load_rules(args.rules), args.output_fasta)
        return
    args = _curate_parser().parse_args(arguments)
    curate_targets(
        queries_path=args.queries,
        subject_path=args.subject,
        rules_path=args.rules,
        hits_path=args.hits,
        fasta_out=args.fasta_out,
        loci_out=args.loci_out,
        expected_classes={"18S", "5.8S", "28S", "5S"},
    )


if __name__ == "__main__":
    main()
