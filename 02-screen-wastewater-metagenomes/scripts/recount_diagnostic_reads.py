#!/usr/bin/env python3
"""Summarize per-sample Deacon screening output into shareable tables.

Two ways to call it:

`--run-dir RUN`   summarize a screen prepared by `prepare_deacon_screen.py`.
                  Every manifest row must have a gzip-valid FASTA, a Deacon JSON
                  summary, and a worker status table, and the thresholds Deacon
                  reported must match `run_metadata.json`. This is the strict
                  path used to sign off a production screen.

`--reads DIR`     summarize any directory of published `*.fasta.gz` candidate
                  reads, with no run metadata required.

Both paths recount each retained read's diagnostic 31-mers against the
core-nt-validated bait set, because Deacon's paired-mode gate is more permissive
than the rule the threshold was calibrated under. Paired mode pools the *distinct*
k-mer hits across both mates and compares that union to `-a`, so a pair can be
retained with neither mate reaching the threshold alone (verified directly:
mates with 13 and 10 disjoint hits are retained at `-a 20`). The step-4
calibration was performed on individual reads, so every retained read is rescored
here and **only reads carrying >= threshold validated k-mers of their own are
published or counted**. A mate emitted because its partner matched, or a pair kept
on the pooled union of both mates, holds no Cyclospora sequence and is dropped: it
would otherwise overstate both the FASTA and every statistic derived from it.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_BAITS = (
    PROJECT_ROOT / "kmers/cyclospora_cayetanensis_rrna_core_nt_validated_baits.fasta"
)
DEACON_VERSION = "deacon 0.15.0"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
KMER_LENGTH = 31

# Every reported statistic describes reads that carry diagnostic k-mers of their
# own. A mate retained only because its partner matched contributes nothing here:
# it holds no Cyclospora sequence, so counting it would overstate the evidence.
SAMPLE_FIELDS = (
    "cohort",
    "sample",
    "status",
    "input_reads",
    "input_pairs",
    "diagnostic_reads",
    "diagnostic_pairs",
    "unique_diagnostic_sequences",
    "max_diagnostic_kmers",
    "diagnostic_read_proportion",
    "abs_threshold",
    "rel_threshold",
    "deacon_seconds",
    "execute_host",
)
READ_FIELDS = (
    "cohort",
    "sample",
    "read_id",
    "read_length",
    "diagnostic_kmers",
)


def canonical(sequence: str) -> str:
    """Return the lexicographically smaller of a sequence and its complement."""
    return min(sequence, sequence.translate(COMPLEMENT)[::-1])


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """Read a plain or gzipped FASTA into (header, sequence) records."""
    opener = gzip.open if path.suffix == ".gz" else open
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence: list[str] = []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence)))
                header = line[1:]
                sequence = []
            elif line:
                sequence.append(line)
    if header is not None:
        records.append((header, "".join(sequence)))
    return records


def load_baits(path: Path) -> set[str]:
    """Load the validated baits as a canonical 31-mer set."""
    baits = set()
    for _, sequence in read_fasta(path):
        sequence = sequence.upper()
        if len(sequence) != KMER_LENGTH or set(sequence) - set("ACGT"):
            raise ValueError(f"bait is not an unambiguous 31-mer: {path}")
        baits.add(canonical(sequence))
    if not baits:
        raise ValueError(f"validated bait FASTA is empty: {path}")
    return baits


def count_diagnostic_kmers(sequence: str, baits: set[str]) -> int:
    """Count positions in a read whose canonical 31-mer is a validated bait."""
    sequence = sequence.upper()
    return sum(
        canonical(sequence[offset : offset + KMER_LENGTH]) in baits
        for offset in range(len(sequence) - KMER_LENGTH + 1)
    )


def read_status(path: Path, abs_threshold: int, rel_threshold: float) -> dict[str, str]:
    """Read and validate one two-column worker status table."""
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["key", "value"]:
            raise ValueError(f"invalid status header: {path}")
        rows = list(reader)
    values = {row["key"]: row["value"] for row in rows}
    if len(values) != len(rows):
        raise ValueError(f"duplicate status key: {path}")
    if values.get("status") != "success":
        raise ValueError(f"worker status is not success: {path}")
    if values.get("deacon_version") != DEACON_VERSION:
        raise ValueError(f"unexpected Deacon version: {path}")
    if int(values["abs_threshold"]) != abs_threshold:
        raise ValueError(f"worker used a different absolute threshold: {path}")
    if float(values["rel_threshold"]) != rel_threshold:
        raise ValueError(f"worker used a different relative threshold: {path}")
    for mate in ("r1", "r2"):
        if values.get(f"{mate}_source_bytes") != values.get(f"{mate}_local_bytes"):
            raise ValueError(f"{mate.upper()} copy size mismatch: {path}")
        if values.get(f"{mate}_argument", "").startswith("/staging/"):
            raise ValueError(f"Deacon received a staging path: {path}")
    return values


def strict_targets(run_dir: Path) -> tuple[list[dict[str, str]], dict, Path]:
    """Return the manifest rows, run metadata, and results directory."""
    metadata = json.loads((run_dir / "run_metadata.json").read_text())
    with (run_dir / "samples.tsv").open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    if not manifest:
        raise ValueError("sample manifest has no rows")
    samples = [row["sample"] for row in manifest]
    if len(samples) != len(set(samples)):
        raise ValueError("sample manifest contains duplicate identifiers")
    if len(samples) != int(metadata["sample_count"]):
        raise ValueError("manifest row count disagrees with run metadata")
    return manifest, metadata, run_dir / "results"


def incomplete_reason(fasta: Path, status: Path, summary: Path) -> str:
    """Return why a sample cannot be summarized, or an empty string if it can."""
    for required in (fasta, status, summary):
        if not required.is_file() or required.stat().st_size == 0:
            return "missing"
    with status.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("key") == "status" and row.get("value") != "success":
                return "failed"
    return ""


def incomplete_row(cohort: str, sample: str, reason: str) -> dict[str, object]:
    """Return a sample row that records an unscreened sample without counts."""
    row: dict[str, object] = {field: "" for field in SAMPLE_FIELDS}
    row.update({"cohort": cohort, "sample": sample, "status": reason})
    return row


def published_targets(reads_dir: Path) -> list[dict[str, str]]:
    """Return one manifest-shaped row per published FASTA."""
    manifest: list[dict[str, str]] = []
    for path in sorted(reads_dir.glob("*.fasta.gz")):
        manifest.append({"cohort": "", "sample": path.name.split(".", 1)[0]})
    if not manifest:
        raise ValueError(f"no *.fasta.gz files were found in {reads_dir}")
    samples = [row["sample"] for row in manifest]
    if len(samples) != len(set(samples)):
        raise ValueError(f"two FASTAs map to the same sample name in {reads_dir}")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--run-dir", type=Path, help="screen directory from prepare_deacon_screen.py"
    )
    source.add_argument(
        "--reads", type=Path, help="directory of published <sample>.*.fasta.gz files"
    )
    parser.add_argument("--baits", type=Path, default=DEFAULT_BAITS)
    parser.add_argument(
        "--abs-threshold",
        type=int,
        default=None,
        help="required with --reads; taken from run metadata with --run-dir",
    )
    parser.add_argument("--sample-output", type=Path, required=True)
    parser.add_argument("--read-output", type=Path, default=None)
    parser.add_argument(
        "--publish-reads",
        type=Path,
        default=None,
        help="write <sample>.diagnostic_reads.fasta.gz holding only the reads that "
        "carry at least the absolute threshold of validated k-mers; mates retained "
        "emitted only because their partner matched are excluded",
    )
    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help="report samples whose worker failed or whose outputs are absent as "
        "status=failed/missing instead of refusing to summarize the screen",
    )
    args = parser.parse_args()

    baits = load_baits(args.baits)

    if args.run_dir:
        run_dir = args.run_dir.resolve()
        manifest, metadata, results = strict_targets(run_dir)
        abs_threshold = int(metadata["abs_threshold"])
        rel_threshold = float(metadata["rel_threshold"])
        if args.abs_threshold is not None and args.abs_threshold != abs_threshold:
            raise ValueError("--abs-threshold disagrees with run metadata")
        suffix = ".diagnostic_pairs.fasta.gz"
    else:
        results = args.reads.resolve()
        manifest = published_targets(results)
        if args.abs_threshold is None:
            raise ValueError("--abs-threshold is required with --reads")
        abs_threshold = args.abs_threshold
        rel_threshold = 0.0
        suffix = None

    sample_rows: list[dict[str, object]] = []
    read_rows: list[dict[str, object]] = []
    incomplete: list[str] = []
    for row in manifest:
        sample = row["sample"]
        cohort = row.get("cohort", "")
        if suffix is not None:
            fasta = results / f"{sample}{suffix}"
        else:
            matches = sorted(results.glob(f"{sample}.*.fasta.gz"))
            if len(matches) != 1:
                raise ValueError(f"expected one FASTA for {sample}, found {matches}")
            fasta = matches[0]

        input_reads: int | str = ""
        input_pairs: int | str = ""
        deacon_seconds: object = ""
        host = ""
        if args.run_dir:
            status_path = results / f"{sample}.status.tsv"
            summary_path = results / f"{sample}.deacon.json"
            reason = incomplete_reason(fasta, status_path, summary_path)
            if reason:
                if not args.allow_incomplete:
                    raise ValueError(
                        f"{sample} is {reason}; pass --allow-incomplete to report "
                        "the rest of the screen without it"
                    )
                incomplete.append(f"{sample}={reason}")
                sample_rows.append(incomplete_row(cohort, sample, reason))
                continue
            status = read_status(status_path, abs_threshold, rel_threshold)
            summary = json.loads(summary_path.read_text())
            if int(summary["abs_threshold"]) != abs_threshold:
                raise ValueError(f"unexpected absolute threshold: {summary_path}")
            if float(summary["rel_threshold"]) != rel_threshold:
                raise ValueError(f"unexpected relative threshold: {summary_path}")
            input_reads = int(summary["seqs_in"])
            if input_reads % 2:
                raise ValueError(f"paired input read count is odd: {summary_path}")
            input_pairs = input_reads // 2
            deacon_seconds = summary["time"]
            host = status["hostname"]
        elif not fasta.is_file() or fasta.stat().st_size == 0:
            raise ValueError(f"result is absent or empty: {fasta}")

        records = read_fasta(fasta)
        if args.run_dir and len(records) != int(summary["seqs_out"]):
            raise ValueError(
                f"FASTA holds {len(records)} reads but Deacon reported "
                f"{summary['seqs_out']}: {fasta}"
            )

        # Deacon's paired mode emits the mate of any hit, and can even retain a
        # pair on the pooled union of both mates. Every retained read is therefore
        # rescored here and only those clearing the threshold on their own are
        # kept, published, or counted.
        diagnostic = 0
        maximum = 0
        sequences: Counter[str] = Counter()
        pairs: set[str] = set()
        published: list[tuple[str, str]] = []
        for header, sequence in records:
            hits = count_diagnostic_kmers(sequence, baits)
            if hits < abs_threshold:
                continue
            read_id = header.split()[0]
            diagnostic += 1
            maximum = max(maximum, hits)
            sequences[canonical(sequence.upper())] += 1
            pairs.add(read_id.removesuffix("/1").removesuffix("/2"))
            published.append(
                (
                    f"{header} diagnostic_kmers={hits} classification=diagnostic",
                    sequence,
                )
            )
            read_rows.append(
                {
                    "cohort": cohort,
                    "sample": sample,
                    "read_id": read_id,
                    "read_length": len(sequence),
                    "diagnostic_kmers": hits,
                }
            )

        if args.publish_reads:
            write_fasta_gz(
                args.publish_reads / f"{sample}.diagnostic_reads.fasta.gz", published
            )

        sample_rows.append(
            {
                "cohort": cohort,
                "sample": sample,
                "status": "success",
                "input_reads": input_reads,
                "input_pairs": input_pairs,
                "diagnostic_reads": diagnostic,
                "diagnostic_pairs": len(pairs),
                "unique_diagnostic_sequences": len(sequences),
                "max_diagnostic_kmers": maximum,
                "diagnostic_read_proportion": (
                    f"{diagnostic / int(input_reads):.6g}"
                    if str(input_reads).isdigit() and int(input_reads)
                    else ""
                ),
                "abs_threshold": abs_threshold,
                "rel_threshold": rel_threshold,
                "deacon_seconds": deacon_seconds,
                "execute_host": host,
            }
        )

    write_tsv(args.sample_output, SAMPLE_FIELDS, sample_rows)
    if args.read_output:
        write_tsv(args.read_output, READ_FIELDS, read_rows)

    screened = [row for row in sample_rows if row["status"] == "success"]

    def total(field: str) -> int:
        return sum(int(row[field]) for row in screened)

    print(f"samples={len(sample_rows)} screened={len(screened)}")
    print(f"positive_samples={sum(1 for row in screened if row['diagnostic_reads'])}")
    if args.run_dir:
        print(f"input_reads={total('input_reads')}")
    print(f"diagnostic_reads={total('diagnostic_reads')}")
    print(f"diagnostic_pairs={total('diagnostic_pairs')}")
    if incomplete:
        print(f"incomplete={len(incomplete)}: {', '.join(incomplete)}")
    print(f"Wrote {args.sample_output}")
    if args.read_output:
        print(f"Wrote {args.read_output}")
    return 0


def write_fasta_gz(path: Path, records: list[tuple[str, str]]) -> None:
    """Write records to a deterministic gzipped FASTA."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            for header, sequence in records:
                zipped.write(f">{header}\n{sequence}\n".encode())
    temporary.replace(path)


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    """Write rows atomically as a tab-separated table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


if __name__ == "__main__":
    raise SystemExit(main())
