#!/usr/bin/env python3
"""Freeze the primary public-SRA outputs from the returned NVD screen."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import tarfile
import tempfile
from pathlib import Path

STAGE_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTDIR = STAGE_ROOT / "results/screen-source"
ARCHIVE_NAME = "prjna1247874_nvd_source.tar.gz"
MANIFEST_NAME = "prjna1247874_nvd_source.manifest.tsv"
EXPECTED_PUBLIC_RUNS = 2_333
EXPECTED_LOCAL_RUNS = 251
EXPECTED_CANDIDATE_READS = 1_476
EXPECTED_SELECTED_RUNS = 81
REPORT_SUFFIX = ".deacon_filter.json"
READ_SUFFIX = ".target_enriched.fastq.gz"
MANIFEST_FIELDS = (
    "archive_path",
    "kind",
    "srr",
    "bytes",
    "sha256",
)
READ_ID = re.compile(r"@(SRR[0-9]+)\.[^/\s]+/([12])(?:\s|$)")
EXPECTED_EXCLUDED = {
    "01_target_enrichment/reads/Columbia_MO_20240220__SRR35904087.target_enriched.deacon-0.15.0-a20.fastq.gz":
        "ac73670af3abed54ac6fb4695131f4099be9fbe39d6076c5d0264a6bbdae9d83",
    "01_target_enrichment/reads/Columbia_MO_20240220__SRR35904087.target_enriched.superseded-a20-20260828.fastq.gz":
        "ac73670af3abed54ac6fb4695131f4099be9fbe39d6076c5d0264a6bbdae9d83",
    "01_target_enrichment/summaries/Columbia_MO_20240220__SRR35904087.deacon_filter.deacon-0.15.0-a20.json":
        "263c0c7c96ec7d53ffcdf53c896e6c442340ff9363090cd5ffd93305bda6ce4f",
    "01_target_enrichment/summaries/Columbia_MO_20240220__SRR35904087.deacon_filter.superseded-a20-20260828.json":
        "263c0c7c96ec7d53ffcdf53c896e6c442340ff9363090cd5ffd93305bda6ce4f",
    "01_target_enrichment/summaries/target_enrichment_summary.tsv":
        "f0b8b17d6558065b4b7607c30a618eb04298c630ff989f0fa52691e6f4ffd67c",
}


def digest(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def read_fastq(path: Path, srr: str) -> int:
    records = []
    with gzip.open(path, "rt", newline="") as handle:
        while header := handle.readline():
            sequence = handle.readline().rstrip("\r\n")
            plus = handle.readline()
            quality = handle.readline().rstrip("\r\n")
            header = header.rstrip("\r\n")
            if not plus or not quality or not plus.startswith("+"):
                raise ValueError(f"malformed FASTQ record in {path}")
            if len(sequence) != len(quality):
                raise ValueError(f"sequence and quality lengths differ in {path}")
            match = READ_ID.match(header)
            if match is None or match.group(1) != srr:
                raise ValueError(f"read header does not match {srr}: {header}")
            records.append((header.split()[0], match.group(2)))
    if len(records) % 2:
        raise ValueError(f"odd number of retained reads in {path}")
    for first, second in zip(records[::2], records[1::2], strict=True):
        if first[1] != "1" or second[1] != "2" or second[0] != f"{first[0][:-2]}/2":
            raise ValueError(f"retained reads are not intact consecutive pairs in {path}")
    return len(records)


def archive_bytes(entries: dict[str, bytes], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
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
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def write_manifest(path: Path, rows: list[dict[str, object]]) -> None:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(output.getvalue())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--replace", action="store_true")
    args = parser.parse_args()
    if args.outdir.exists() and not args.replace:
        raise SystemExit("screen-source output exists; pass --replace to rebuild it")

    samplesheet = args.run_dir / "prjna1247874_plus_local_fastqs_samplesheet.csv"
    with samplesheet.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    public = [row for row in rows if row.get("srr")]
    local = [row for row in rows if not row.get("srr")]
    if len(public) != EXPECTED_PUBLIC_RUNS or len(local) != EXPECTED_LOCAL_RUNS:
        raise SystemExit(f"expected {EXPECTED_PUBLIC_RUNS} public and {EXPECTED_LOCAL_RUNS} local rows")
    if len({row["sample_id"] for row in rows}) != len(rows):
        raise SystemExit("samplesheet contains duplicate sample identifiers")
    if len({row["srr"] for row in public}) != len(public):
        raise SystemExit("samplesheet contains duplicate SRA accessions")
    if {row["platform"] for row in public} != {"illumina"}:
        raise SystemExit("public screen must contain only Illumina runs")

    entries: dict[str, bytes] = {}
    details: dict[str, tuple[str, str]] = {}
    report_dir = args.run_dir / "01_target_enrichment/summaries"
    read_dir = args.run_dir / "01_target_enrichment/reads"
    candidate_reads = selected_runs = 0
    for row in sorted(public, key=lambda item: item["srr"]):
        sample, srr = row["sample_id"], row["srr"]
        report_path = report_dir / f"{sample}{REPORT_SUFFIX}"
        reads_path = read_dir / f"{sample}{READ_SUFFIX}"
        report = json.loads(report_path.read_text())
        expected = {
            "version": "deacon 0.16.0",
            "index": "taxonomically_screened.idx",
            "input": "-",
            "input2": None,
            "output": reads_path.name,
            "output2": None,
            "k": 31,
            "w": 1,
            "abs_threshold": 24,
            "rel_threshold": 0.0,
            "check_pairs": False,
        }
        for field, value in expected.items():
            if report.get(field) != value:
                raise SystemExit(f"unexpected {field} in {report_path}: {report.get(field)!r}")
        records = read_fastq(reads_path, srr)
        if report.get("seqs_out") != records:
            raise SystemExit(f"JSON/FASTQ read-count mismatch for {sample}")
        candidate_reads += records
        selected_runs += records > 0
        for kind, source, name in (
            ("deacon_report", report_path, f"reports/{report_path.name}"),
            ("candidate_reads", reads_path, f"reads/{reads_path.name}"),
        ):
            entries[name] = source.read_bytes()
            details[name] = (kind, srr)

    if candidate_reads != EXPECTED_CANDIDATE_READS or selected_runs != EXPECTED_SELECTED_RUNS:
        raise SystemExit(
            f"expected {EXPECTED_CANDIDATE_READS} reads across {EXPECTED_SELECTED_RUNS} runs; "
            f"found {candidate_reads} across {selected_runs}"
        )

    provenance = {
        "provenance/run_command.sh": args.run_dir / "run_command.sh",
        "provenance/nvd_version.txt": args.run_dir / "nvd_version.txt",
        "provenance/taxonomically_screened_baits.fasta": args.run_dir / "cyclospora_cayetanensis_rrna/04_bait_sets/taxonomically_screened_baits.fasta",
        "provenance/taxonomically_screened.idx": args.run_dir / "cyclospora_cayetanensis_rrna/05_deacon_index/taxonomically_screened.idx",
    }
    for name, source in provenance.items():
        entries[name] = source.read_bytes()
        details[name] = ("provenance", "")

    expected_samples = {row["sample_id"] for row in rows}
    expected_outputs = {
        *(report_dir / f"{sample}{REPORT_SUFFIX}" for sample in expected_samples),
        *(read_dir / f"{sample}{READ_SUFFIX}" for sample in expected_samples),
    }
    # Postmerge recounts are downstream audit work, not returned NVD output.
    postmerge = report_dir / "postmerge"
    outputs = set(report_dir.iterdir()) | set(read_dir.iterdir())
    extras = sorted(outputs - expected_outputs - ({postmerge} if postmerge.is_dir() else set()))
    if any(not path.is_file() for path in extras):
        raise SystemExit(f"unexpected non-primary outputs: {[str(path) for path in extras]}")
    excluded = {
        str(path.relative_to(args.run_dir)): digest(path.read_bytes())
        for path in extras
    }
    if excluded != EXPECTED_EXCLUDED:
        raise SystemExit(f"unexpected non-primary outputs: {sorted(excluded)}")

    manifest = []
    for name in sorted(entries):
        kind, srr = details[name]
        manifest.append(
            {
                "archive_path": name,
                "kind": kind,
                "srr": srr,
                "bytes": len(entries[name]),
                "sha256": digest(entries[name]),
            }
        )

    args.outdir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{args.outdir.name}.", dir=args.outdir.parent))
    backup = None
    try:
        archive = temporary / ARCHIVE_NAME
        manifest_path = temporary / MANIFEST_NAME
        archive_bytes(entries, archive)
        write_manifest(manifest_path, manifest)
        if args.outdir.exists():
            backup = Path(tempfile.mkdtemp(prefix=f".{args.outdir.name}.old.", dir=args.outdir.parent))
            backup.rmdir()
            args.outdir.replace(backup)
        try:
            temporary.replace(args.outdir)
        except BaseException:
            if backup is not None:
                if args.outdir.exists():
                    shutil.rmtree(args.outdir)
                backup.replace(args.outdir)
            raise
        if backup is not None:
            shutil.rmtree(backup)
    finally:
        shutil.rmtree(temporary, ignore_errors=True)
    print(f"public SRA runs: {len(public):,}; local rows excluded: {len(local):,}")
    print(f"candidate reads: {candidate_reads:,} across {selected_runs} runs")
    print(f"wrote {args.outdir / ARCHIVE_NAME} ({(args.outdir / ARCHIVE_NAME).stat().st_size:,} bytes)")
    print(f"wrote {args.outdir / MANIFEST_NAME}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        raise SystemExit(str(error)) from error
