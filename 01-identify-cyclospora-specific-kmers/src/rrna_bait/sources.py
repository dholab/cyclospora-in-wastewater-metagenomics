from __future__ import annotations

import argparse
import csv
import gzip
import io
import shutil
import sys
from collections.abc import Callable, Sequence
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from rrna_bait.core import normalize_dna, open_text, sha256_file

_PROVENANCE_FIELDS = [
    "name",
    "release",
    "url",
    "local_path",
    "retrieved_at_utc",
    "bytes",
    "sha256",
]


def accession_root(header: str) -> str:
    token = header.split()[0].split("/")[0]
    if token.count(".") >= 2 and token.rsplit(".", 1)[1].isdigit():
        token = token.split(".", 1)[0]
    else:
        token = token.rsplit(".", 1)[0]
    return token


def verify_expected_identity(
    path: Path, *, expected_bytes: int, expected_sha256: str
) -> None:
    """Require a file to match immutable identity recorded in configuration."""
    path = Path(path)
    observed_bytes = path.stat().st_size
    observed_sha256 = sha256_file(path)
    if (
        observed_bytes != expected_bytes
        or observed_sha256 != expected_sha256
    ):
        raise ValueError(
            f"expected identity mismatch for {path}: expected "
            f"{expected_bytes} bytes / {expected_sha256}, observed "
            f"{observed_bytes} bytes / {observed_sha256}"
        )


def atomic_download(
    url: str,
    destination: Path,
    *,
    expected_bytes: int | None = None,
    expected_sha256: str | None = None,
) -> None:
    """Download a non-empty file without replacing a valid destination early."""
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    try:
        with urlopen(url) as response, temporary.open("wb") as handle:
            shutil.copyfileobj(response, handle)
        if temporary.stat().st_size == 0:
            raise ValueError(f"Downloaded response from {url} is empty")
        if expected_bytes is not None or expected_sha256 is not None:
            if expected_bytes is None or expected_sha256 is None:
                raise ValueError(
                    "expected_bytes and expected_sha256 must be supplied together"
                )
            verify_expected_identity(
                temporary,
                expected_bytes=expected_bytes,
                expected_sha256=expected_sha256,
            )
        temporary.replace(destination)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def reconstruct_gzip_concat(
    sources: Sequence[Path],
    destination: Path,
    *,
    expected_bytes: int,
    expected_sha256: str,
) -> None:
    """Atomically decompress and concatenate pinned gzip inputs in order."""
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".tmp")
    temporary.unlink(missing_ok=True)
    try:
        with temporary.open("wb") as output_handle:
            for source in sources:
                with gzip.open(Path(source), "rb") as input_handle:
                    shutil.copyfileobj(input_handle, output_handle)
        verify_expected_identity(
            temporary,
            expected_bytes=expected_bytes,
            expected_sha256=expected_sha256,
        )
        temporary.replace(destination)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def _write_filtered_fasta(
    source: Path,
    destination: Path,
    should_remove: Callable[[SeqRecord], bool],
) -> tuple[int, int]:
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    removed = 0
    with open_text(Path(source)) as input_handle, destination.open("wb") as raw_output:
        with gzip.GzipFile(
            filename="", fileobj=raw_output, mode="wb", mtime=0
        ) as gzip_output:
            with io.TextIOWrapper(gzip_output, encoding="utf-8", newline="\n") as output_handle:
                for record in SeqIO.parse(input_handle, "fasta"):
                    if should_remove(record):
                        removed += 1
                        continue
                    record.seq = Seq(normalize_dna(str(record.seq)))
                    SeqIO.write(record, output_handle, "fasta")
                    kept += 1
    return kept, removed


def _silva_lineage_is_target(description: str, excluded_taxon: str) -> bool:
    fields = description.split(maxsplit=1)
    lineage = fields[1] if len(fields) == 2 else ""
    return lineage.split(";")[-1].strip() == excluded_taxon


def filter_silva(
    source: Path,
    destination: Path,
    excluded_accessions: set[str],
    excluded_taxon: str = "Cyclospora cayetanensis",
) -> tuple[int, int]:
    """Write non-target SILVA records and return their kept/removed counts."""
    excluded_roots = {accession_root(accession) for accession in excluded_accessions}
    return _write_filtered_fasta(
        source,
        destination,
        lambda record: accession_root(record.description) in excluded_roots
        or _silva_lineage_is_target(record.description, excluded_taxon),
    )


def filter_rfam(
    source: Path,
    destination: Path,
    excluded_accessions: set[str],
) -> tuple[int, int]:
    """Write non-target Rfam records and return their kept/removed counts."""
    excluded_roots = {accession_root(accession) for accession in excluded_accessions}
    return _write_filtered_fasta(
        source,
        destination,
        lambda record: accession_root(record.description) in excluded_roots,
    )


def write_provenance_row(
    path: Path,
    name: str,
    release: str,
    url: str,
    local_path: Path,
    project_root: Path | None = None,
) -> None:
    """Upsert a source's local acquisition details in a name-sorted TSV."""
    path = Path(path)
    local_path = Path(local_path)
    rows_by_name: dict[str, dict[str, str]] = {}
    if path.exists():
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("name"):
                    rows_by_name[row["name"]] = {
                        field: row.get(field, "") for field in _PROVENANCE_FIELDS
                    }

    stored_path = (
        str(local_path.resolve().relative_to(Path(project_root).resolve()))
        if project_root is not None
        else str(local_path)
    )
    rows_by_name = {
        source_name: row
        for source_name, row in rows_by_name.items()
        if source_name == name or row.get("local_path") != stored_path
    }
    rows_by_name[name] = {
        "name": name,
        "release": release,
        "url": url,
        "local_path": stored_path,
        "retrieved_at_utc": datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z"),
        "bytes": str(local_path.stat().st_size),
        "sha256": sha256_file(local_path),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=_PROVENANCE_FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows_by_name[source_name] for source_name in sorted(rows_by_name))


def verify_provenance_row(
    path: Path,
    name: str,
    release: str,
    url: str,
    local_path: Path,
    project_root: Path | None = None,
    expected_bytes: int | None = None,
    expected_sha256: str | None = None,
) -> bool:
    """Return whether a cache exactly matches its named provenance row."""
    path = Path(path)
    local_path = Path(local_path)
    if not path.is_file() or not local_path.is_file():
        return False

    try:
        with path.open(newline="") as handle:
            matching_rows = [
                row
                for row in csv.DictReader(handle, delimiter="\t")
                if row.get("name") == name
            ]
        if len(matching_rows) != 1:
            return False
        row = matching_rows[0]
        stored_path = (
            str(local_path.resolve().relative_to(Path(project_root).resolve()))
            if project_root is not None
            else str(local_path)
        )
        if expected_bytes is not None or expected_sha256 is not None:
            if expected_bytes is None or expected_sha256 is None:
                return False
            verify_expected_identity(
                local_path,
                expected_bytes=expected_bytes,
                expected_sha256=expected_sha256,
            )
        return (
            row.get("release") == release
            and row.get("url") == url
            and row.get("local_path") == stored_path
            and bool(row.get("retrieved_at_utc"))
            and row.get("bytes") == str(local_path.stat().st_size)
            and row.get("sha256") == sha256_file(local_path)
        )
    except (OSError, ValueError):
        return False


def _filter_silva_parser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("filter-silva")
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--excluded-accession", action="append", default=[])
    parser.add_argument("--excluded-taxon", default="Cyclospora cayetanensis")


def _filter_rfam_parser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("filter-rfam")
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--excluded-accession", action="append", default=[])


def _record_provenance_parser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("record-provenance")
    parser.add_argument("path", type=Path)
    parser.add_argument("--name", required=True)
    parser.add_argument("--release", required=True)
    parser.add_argument("--url", required=True)
    parser.add_argument("--local-path", type=Path, required=True)
    parser.add_argument("--project-root", type=Path)


def _verify_provenance_parser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("verify-provenance")
    parser.add_argument("path", type=Path)
    parser.add_argument("--name", required=True)
    parser.add_argument("--release", required=True)
    parser.add_argument("--url", required=True)
    parser.add_argument("--local-path", type=Path, required=True)
    parser.add_argument("--project-root", type=Path)
    parser.add_argument("--expected-bytes", type=int)
    parser.add_argument("--expected-sha256")


def _verify_identity_parser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("verify-identity")
    parser.add_argument("path", type=Path)
    parser.add_argument("--expected-bytes", type=int, required=True)
    parser.add_argument("--expected-sha256", required=True)


def _reconstruct_parser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("reconstruct-gzip-concat")
    parser.add_argument("destination", type=Path)
    parser.add_argument("sources", nargs="+", type=Path)
    parser.add_argument("--expected-bytes", type=int, required=True)
    parser.add_argument("--expected-sha256", required=True)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare pinned non-target rRNA sources.")
    subparsers = parser.add_subparsers(dest="command", required=True)
    _filter_silva_parser(subparsers)
    _filter_rfam_parser(subparsers)
    _record_provenance_parser(subparsers)
    _verify_provenance_parser(subparsers)
    _verify_identity_parser(subparsers)
    _reconstruct_parser(subparsers)
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    args = _parser().parse_args(sys.argv[1:] if argv is None else argv)
    if args.command == "filter-silva":
        filter_silva(
            args.source,
            args.destination,
            set(args.excluded_accession),
            args.excluded_taxon,
        )
    elif args.command == "filter-rfam":
        filter_rfam(args.source, args.destination, set(args.excluded_accession))
    elif args.command == "record-provenance":
        write_provenance_row(
            args.path,
            args.name,
            args.release,
            args.url,
            args.local_path,
            args.project_root,
        )
    elif args.command == "verify-identity":
        try:
            verify_expected_identity(
                args.path,
                expected_bytes=args.expected_bytes,
                expected_sha256=args.expected_sha256,
            )
        except (OSError, ValueError) as error:
            raise SystemExit(str(error)) from error
    elif args.command == "reconstruct-gzip-concat":
        try:
            reconstruct_gzip_concat(
                args.sources,
                args.destination,
                expected_bytes=args.expected_bytes,
                expected_sha256=args.expected_sha256,
            )
        except (OSError, ValueError) as error:
            raise SystemExit(str(error)) from error
    elif not verify_provenance_row(
        args.path,
        args.name,
        args.release,
        args.url,
        args.local_path,
        args.project_root,
        args.expected_bytes,
        args.expected_sha256,
    ):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
