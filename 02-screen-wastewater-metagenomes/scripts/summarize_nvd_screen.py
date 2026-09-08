#!/usr/bin/env python3
"""Recount the frozen public NVD screen and write publication inputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import shutil
import tarfile
import tempfile
from pathlib import Path

STAGE_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_ARCHIVE = STAGE_ROOT / "results/screen-source/prjna1247874_nvd_source.tar.gz"
DEFAULT_MANIFEST = STAGE_ROOT / "results/screen-source/prjna1247874_nvd_source.manifest.tsv"
DEFAULT_COHORT = STAGE_ROOT / "config/wastewater_cohort.tsv"
DEFAULT_BAITS = (
    STAGE_ROOT.parent
    / "01-identify-cyclospora-specific-kmers/baits"
    / "cyclospora_cayetanensis_rrna_core_nt_validated_baits.fasta"
)
THRESHOLD = 24
K = 31
COMPLEMENT = str.maketrans("ACGT", "TGCA")
SUMMARY_FIELDS = (
    "public_id",
    "casper_code",
    "sra_accession",
    "collection_date",
    "input_reads",
    "diagnostic_reads",
    "distinct_diagnostic_reads",
    "max_diagnostic_kmers",
    "distinct_diagnostic_reads_per_billion",
    "positive",
)
EXPECTED = {
    "runs": 2_333,
    "candidate_reads": 1_476,
    "diagnostic_reads": 1_244,
    "positive_runs": 81,
}
EXPECTED_ARCHIVE_SHA256 = "5313940b512b9bd8010e96cd6c8f6a2448cc18f57589c55d78f9e22cb95c6cdd"
EXPECTED_MANIFEST_SHA256 = "f72328095f2ac921bec249253b470fad1db97f1341cdc87cd025fe4b59de733f"
EXPECTED_COHORT_SHA256 = "1df452f80b3f8837ee45f0fcee6583b8e0f875ec8d5efbed5964b8ee969a4a61"


def digest(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def canonical(sequence: str) -> str:
    reverse = sequence.translate(COMPLEMENT)[::-1]
    return min(sequence, reverse)


def fasta_sequences(content: bytes) -> list[str]:
    sequences = []
    parts: list[str] = []
    for raw_line in content.decode().splitlines():
        line = raw_line.strip()
        if line.startswith(">"):
            if parts:
                sequences.append("".join(parts).upper())
                parts = []
        elif line:
            parts.append(line)
    if parts:
        sequences.append("".join(parts).upper())
    if not sequences:
        raise ValueError("FASTA contains no sequences")
    return sequences


def read_fastq(content: bytes) -> list[tuple[str, str]]:
    text = gzip.decompress(content).decode()
    lines = text.splitlines()
    if len(lines) % 4:
        raise ValueError("candidate FASTQ does not contain complete records")
    records = []
    for offset in range(0, len(lines), 4):
        header, sequence, plus, quality = lines[offset:offset + 4]
        if not header.startswith("@") or not plus.startswith("+") or len(sequence) != len(quality):
            raise ValueError(f"malformed FASTQ record: {header}")
        records.append((header[1:].split()[0], sequence.upper()))
    return records


def matching_baits(sequence: str, baits: set[str]) -> set[str]:
    matches = set()
    for start in range(len(sequence) - K + 1):
        kmer = sequence[start:start + K]
        if "N" not in kmer:
            kmer = canonical(kmer)
            if kmer in baits:
                matches.add(kmer)
    return matches


def load_source(archive_path: Path, manifest_path: Path) -> tuple[dict[str, bytes], list[dict[str, str]]]:
    with manifest_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected_fields = {
            "archive_path", "kind", "srr", "bytes", "sha256",
        }
        if set(reader.fieldnames or ()) != expected_fields:
            raise ValueError("source manifest has unexpected columns")
        manifest = list(reader)
    if not manifest or len({row["archive_path"] for row in manifest}) != len(manifest):
        raise ValueError("source manifest is empty or contains duplicate paths")
    expected_paths = {row["archive_path"] for row in manifest}
    with tarfile.open(archive_path, "r:gz") as archive:
        members = archive.getmembers()
        member_names = [member.name for member in members]
        if any(not member.isfile() for member in members):
            raise ValueError("source archive may contain only regular files")
        if len(set(member_names)) != len(member_names) or set(member_names) != expected_paths:
            raise ValueError("archive members do not match the source manifest")
        source = {}
        for member in members:
            if member.mtime != 0 or member.uid != 0 or member.gid != 0:
                raise ValueError(f"archive metadata is not normalized: {member.name}")
            extracted = archive.extractfile(member)
            if extracted is None:
                raise ValueError(f"cannot read archive member {member.name}")
            source[member.name] = extracted.read()
    for row in manifest:
        content = source[row["archive_path"]]
        if len(content) != int(row["bytes"]) or digest(content) != row["sha256"]:
            raise ValueError(f"archive member fails its manifest: {row['archive_path']}")
    return source, manifest


def fasta_bytes(records: list[tuple[str, str, int]]) -> bytes:
    return "".join(
        f">{read_id} diagnostic_kmers={count}\n{sequence}\n"
        for read_id, sequence, count in records
    ).encode()


def archive_bytes(entries: dict[str, bytes]) -> bytes:
    output = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=output, mtime=0) as compressed:
        with tarfile.open(fileobj=compressed, mode="w|", format=tarfile.PAX_FORMAT) as archive:
            for name in sorted(entries):
                content = entries[name]
                info = tarfile.TarInfo(name)
                info.size = len(content)
                info.mtime = 0
                info.mode = 0o644
                info.uid = info.gid = 0
                info.uname = info.gname = ""
                archive.addfile(info, io.BytesIO(content))
    return output.getvalue()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=STAGE_ROOT / "results")
    args = parser.parse_args()

    archive_hash = digest(DEFAULT_ARCHIVE.read_bytes())
    manifest_hash = digest(DEFAULT_MANIFEST.read_bytes())
    if archive_hash != EXPECTED_ARCHIVE_SHA256 or manifest_hash != EXPECTED_MANIFEST_SHA256:
        raise SystemExit("source archive or manifest does not match the accepted screen")
    source, manifest = load_source(DEFAULT_ARCHIVE, DEFAULT_MANIFEST)
    report_list = [row for row in manifest if row["kind"] == "deacon_report"]
    read_list = [row for row in manifest if row["kind"] == "candidate_reads"]
    if len(report_list) != EXPECTED["runs"] or len(read_list) != EXPECTED["runs"]:
        raise SystemExit("source bundle does not contain one report/FASTQ pair per public run")
    if len({row["srr"] for row in report_list}) != len(report_list):
        raise SystemExit("source manifest contains duplicate report accessions")
    if len({row["srr"] for row in read_list}) != len(read_list):
        raise SystemExit("source manifest contains duplicate FASTQ accessions")
    report_rows = {row["srr"]: row for row in report_list}
    read_rows = {row["srr"]: row for row in read_list}
    if set(report_rows) != set(read_rows):
        raise SystemExit("source bundle does not contain one report/FASTQ pair per public run")

    with DEFAULT_COHORT.open(newline="") as handle:
        cohort_rows = list(csv.DictReader(handle, delimiter="\t"))
    if digest(DEFAULT_COHORT.read_bytes()) != EXPECTED_COHORT_SHA256:
        raise SystemExit("cohort metadata does not match the accepted screen")
    cohort = {row["srr"]: row for row in cohort_rows}
    if len(cohort) != len(cohort_rows) or set(report_rows) != set(cohort):
        raise SystemExit("cohort metadata does not cover each screened accession exactly once")
    for row in cohort_rows:
        public_id = f"{row['casper_code']}_{row['collection_date'].replace('-', '')}"
        if row["sample"] != f"{public_id}__{row['srr']}":
            raise SystemExit(f"invalid canonical cohort identity for {row['srr']}")

    executed_bait_content = source["provenance/taxonomically_screened_baits.fasta"]
    executed_baits = {canonical(sequence) for sequence in fasta_sequences(executed_bait_content)}
    committed_baits = {canonical(sequence) for sequence in fasta_sequences(DEFAULT_BAITS.read_bytes())}
    if executed_baits != committed_baits or len(executed_baits) != 1_464:
        raise SystemExit("executed and committed bait sequence sets differ")

    summary = []
    published: dict[str, bytes] = {}
    candidate_reads = 0
    positive_runs = set()
    for srr in sorted(report_rows):
        report_row, read_row = report_rows[srr], read_rows[srr]
        report_sample = Path(report_row["archive_path"]).name.removesuffix(".deacon_filter.json")
        read_sample = Path(read_row["archive_path"]).name.removesuffix(".target_enriched.fastq.gz")
        if report_sample != read_sample:
            raise SystemExit(f"source sample identifiers differ for {srr}")
        metadata = cohort[srr]
        report = json.loads(source[report_row["archive_path"]])
        if (
            report.get("version") != "deacon 0.16.0"
            or report.get("k") != K
            or report.get("w") != 1
            or report.get("abs_threshold") != THRESHOLD
            or report.get("rel_threshold") != 0.0
        ):
            raise SystemExit(f"unexpected Deacon settings for {srr}")
        records = read_fastq(source[read_row["archive_path"]])
        if len(records) != report.get("seqs_out"):
            raise SystemExit(f"report/FASTQ read-count mismatch for {srr}")
        candidate_reads += len(records)
        accepted = []
        for read_id, sequence in records:
            hits = matching_baits(sequence, executed_baits)
            if len(hits) >= THRESHOLD:
                accepted.append((read_id, sequence, len(hits)))
        if accepted:
            positive_runs.add(srr)

        public_id = f"{metadata['casper_code']}_{metadata['collection_date'].replace('-', '')}"
        output_name = f"reads/{public_id}__{srr}.diagnostic_reads.fasta"
        if accepted:
            published[output_name] = fasta_bytes(accepted)
        distinct = len({canonical(sequence) for _, sequence, _ in accepted})
        maximum = max((count for _, _, count in accepted), default=0)
        input_reads = int(report["seqs_in"])
        summary.append(
            {
                "public_id": public_id,
                "casper_code": metadata["casper_code"],
                "sra_accession": srr,
                "collection_date": metadata["collection_date"],
                "input_reads": input_reads,
                "diagnostic_reads": len(accepted),
                "distinct_diagnostic_reads": distinct,
                "max_diagnostic_kmers": maximum,
                "distinct_diagnostic_reads_per_billion": f"{distinct / input_reads * 1e9:.3f}",
                "positive": "yes" if accepted else "no",
            }
        )

    observed = {
        "runs": len(summary),
        "candidate_reads": candidate_reads,
        "diagnostic_reads": sum(int(row["diagnostic_reads"]) for row in summary),
        "positive_runs": len(positive_runs),
    }
    if observed != EXPECTED:
        raise SystemExit(f"screen recount differs from the accepted result: {observed}")

    summary.sort(key=lambda row: (row["casper_code"], row["collection_date"], row["sra_accession"]))
    summary_output = io.StringIO(newline="")
    writer = csv.DictWriter(summary_output, fieldnames=SUMMARY_FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary)
    summary_content = summary_output.getvalue().encode()
    read_archive_content = archive_bytes(published)

    receipt = {
        "inputs": {
            "source_archive_sha256": archive_hash,
            "source_manifest_sha256": manifest_hash,
            "cohort_sha256": digest(DEFAULT_COHORT.read_bytes()),
            "committed_baits_sha256": digest(DEFAULT_BAITS.read_bytes()),
            "executed_baits_sha256": digest(executed_bait_content),
        },
        "rule": {
            "k": K,
            "threshold": THRESHOLD,
        },
        "results": {
            **observed,
            "negative_runs": len(summary) - len(positive_runs),
        },
        "outputs": {
            "sra_sample_summary.tsv": digest(summary_content),
            "diagnostic_reads.tar.gz": digest(read_archive_content),
            "diagnostic_read_members": {
                name: digest(content) for name, content in sorted(published.items())
            },
        },
    }
    receipt_content = (json.dumps(receipt, indent=2, sort_keys=True) + "\n").encode()

    args.outdir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".nvd-screen.", dir=args.outdir.parent))
    previous = temporary / "previous"
    installed: list[Path] = []
    moved: list[Path] = []
    preserve_backup = False
    try:
        generated = temporary / "generated"
        generated.mkdir()
        (generated / "sra_sample_summary.tsv").write_bytes(summary_content)
        (generated / "nvd_screen_receipt.json").write_bytes(receipt_content)
        (generated / "diagnostic_reads.tar.gz").write_bytes(read_archive_content)
        args.outdir.mkdir(parents=True, exist_ok=True)
        previous.mkdir()
        targets = (
            "diagnostic_reads.tar.gz",
            "sra_sample_summary.tsv",
            "nvd_screen_receipt.json",
        )
        try:
            for name in targets:
                target = args.outdir / name
                if target.exists():
                    target.replace(previous / name)
                    moved.append(previous / name)
            for name in targets:
                target = args.outdir / name
                (generated / name).replace(target)
                installed.append(target)
        except BaseException:
            try:
                for target in installed:
                    target.unlink()
                for old in moved:
                    old.replace(args.outdir / old.name)
            except BaseException as error:
                preserve_backup = True
                error.add_note(f"Rollback failed; recovery files preserved at {temporary}")
                raise
            raise
    finally:
        if not preserve_backup:
            shutil.rmtree(temporary, ignore_errors=True)

    print(f"public SRA runs: {len(summary):,}")
    print(f"candidate reads: {candidate_reads:,}")
    print(f"diagnostic reads reaching 24: {observed['diagnostic_reads']:,} across {len(positive_runs)} runs")
    print(f"wrote {args.outdir / 'sra_sample_summary.tsv'} and an archive with {len(published)} FASTAs")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, tarfile.TarError) as error:
        raise SystemExit("\n".join((str(error), *getattr(error, "__notes__", ())))) from error
