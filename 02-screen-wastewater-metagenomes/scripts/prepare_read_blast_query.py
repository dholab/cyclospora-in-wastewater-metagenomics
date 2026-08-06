#!/usr/bin/env python3
"""Build the core-nt BLAST query from `-a 1` candidate reads.

A validated 31-mer is necessary but not sufficient evidence about a read (see the
README, step 4), so candidate reads are classified against core-nt along their
full length before any threshold is chosen. This prepares that query, and in
doing so fixes the read population the calibration is computed over.

Deacon's paired mode emits two kinds of sequence that must not reach the sweep:

1. **Mates carrying no diagnostic 31-mer of their own.** Paired mode pools the
   distinct k-mer hits across both mates, so a read with nothing can be emitted
   because its partner matched. The `-a` rule is per-read, so such a read could
   never pass any threshold; leaving it in the denominator would understate the
   fraction of candidates that are genuinely *Cyclospora*, and BLASTing it spends
   database time on a read no threshold could keep.
2. **Duplicate molecules.** A fragment that was amplified or resequenced is one
   observation. Two fragments count as one only when *both* mates match base for
   base — a shared R1 with different R2s is two molecules and stays two.

Order matters: fragments are collapsed first, then the k-mer-free mates are
dropped. Dropping first would make two distinct molecules that happen to share
their k-mer-bearing mate look like a single duplicate.

Each read is recounted here against the bait FASTA from stage 01 rather than
trusting the screen, so this runs on a bare checkout with no Deacon, no network,
and no database. The counts it produces are the same ones `deacon filter --debug`
reports, because the index is built with `w = 1` and both count distinct
canonical bait 31-mers per read.

Standard library only.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
from pathlib import Path

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def canonical(sequence: str) -> str:
    """Return the lexicographically smaller of a sequence and its complement."""
    sequence = sequence.upper()
    return min(sequence, sequence.translate(COMPLEMENT)[::-1])


def read_fasta(path: Path):
    """Yield (header, sequence) from a plain or gzipped FASTA."""
    opener = gzip.open if path.suffix == ".gz" else open
    header: str | None = None
    sequence: list[str] = []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line[1:]
                sequence = []
            elif line:
                sequence.append(line)
    if header is not None:
        yield header, "".join(sequence)


def load_baits(path: Path) -> tuple[set[str], int]:
    """Return the canonical bait k-mers and k, from the stage 01 FASTA."""
    baits: set[str] = set()
    lengths: set[int] = set()
    for _, sequence in read_fasta(path):
        baits.add(canonical(sequence))
        lengths.add(len(sequence))
    if not baits:
        raise SystemExit(f"no baits in {path}")
    if len(lengths) != 1:
        raise SystemExit(f"{path} mixes k-mer lengths: {sorted(lengths)}")
    return baits, lengths.pop()


def count_kmers(sequence: str, baits: set[str], k: int) -> int:
    """Distinct diagnostic k-mers the read carries on its own."""
    sequence = sequence.upper()
    return len({canonical(sequence[i:i + k])
                for i in range(len(sequence) - k + 1)} & baits)


def mate_suffix(header: str) -> str:
    """Return `/1` or `/2` for a read, from its name or its Illumina comment."""
    name = header.split()[0]
    if name.endswith("/1") or name.endswith("/2"):
        return ""
    fields = header.split()
    if len(fields) > 1:
        for field in fields[1:]:
            # SRA writes the original instrument name, which ends `/1` or `/2`.
            if field.endswith("/1") or field.endswith("/2"):
                return f"/{field[-1]}"
        if fields[1][:1] in {"1", "2"}:
            return f"/{fields[1][0]}"
    return "/1"


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    stage = Path(__file__).resolve().parent.parent
    parser.add_argument("--reads", type=Path, nargs="+", required=True,
                        help="directories of <sample>.*.fasta.gz candidate reads")
    parser.add_argument(
        "--baits", type=Path,
        default=stage.parent / "01-identify-cyclospora-specific-kmers/baits"
                               "/cyclospora_cayetanensis_rrna_core_nt_validated_baits.fasta",
        help="the core-nt-validated bait FASTA from stage 01")
    parser.add_argument("--query-output", type=Path, required=True)
    parser.add_argument("--map-output", type=Path, required=True)
    parser.add_argument("--min-kmers", type=int, default=1,
                        help="a read must carry at least this many diagnostic "
                             "31-mers of its own to be kept at all")
    args = parser.parse_args()

    baits, k = load_baits(args.baits)

    # --- read every candidate, recounting its k-mers -------------------------
    fragments: dict[tuple[str, str], dict[str, dict]] = collections.defaultdict(dict)
    total_reads = 0
    for directory in args.reads:
        for path in sorted(directory.glob("*.fasta.gz")):
            sample = path.name.split(".", 1)[0]
            for header, sequence in read_fasta(path):
                suffix = mate_suffix(header)
                name = header.split()[0]
                read_id = f"{sample}|{name}{suffix}"
                mate = read_id[-1] if read_id[-2] == "/" else "1"
                stem = read_id.rsplit("/", 1)[0]
                fragments[(sample, stem)][mate] = {
                    "read_id": read_id,
                    "sample": sample,
                    "kmers": count_kmers(sequence, baits, k),
                    "canonical": canonical(sequence),
                    "length": len(sequence),
                }
                total_reads += 1
    if not total_reads:
        raise SystemExit("no candidate reads were found in the given directories")

    # --- collapse duplicate fragments, then drop the empty mates -------------
    seen: set[tuple] = set()
    kept: list[dict] = []
    duplicate_fragments = no_kmer = 0
    for (sample, stem), mates in sorted(fragments.items()):
        signature = (sample,) + tuple(mates[m]["canonical"] for m in sorted(mates))
        if signature in seen:
            duplicate_fragments += 1
            continue
        seen.add(signature)
        for mate in sorted(mates):
            record = mates[mate]
            if record["kmers"] < args.min_kmers:
                no_kmer += 1
                continue
            kept.append(record)
    if not kept:
        raise SystemExit("no read carried a diagnostic k-mer of its own")

    # --- one BLAST query per distinct sequence -------------------------------
    representative: dict[str, str] = {}
    for record in kept:
        representative.setdefault(record["canonical"], record["read_id"])

    args.query_output.parent.mkdir(parents=True, exist_ok=True)
    with args.query_output.open("w") as handle:
        for sequence, read_id in representative.items():
            handle.write(f">{read_id}\n{sequence}\n")

    args.map_output.parent.mkdir(parents=True, exist_ok=True)
    with args.map_output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("sample", "read_id", "representative_id", "deacon_hits",
                        "read_length"),
            delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for record in kept:
            writer.writerow({
                "sample": record["sample"],
                "read_id": record["read_id"],
                "representative_id": representative[record["canonical"]],
                "deacon_hits": record["kmers"],
                "read_length": record["length"],
            })

    print(f"baits={len(baits)} k={k}")
    print(f"candidate reads        {total_reads:>8,}")
    print(f"fragments              {len(fragments):>8,}")
    print(f"duplicate fragments    {duplicate_fragments:>8,}  (dropped)")
    print(f"mates without a k-mer  {no_kmer:>8,}  (dropped)")
    print(f"reads kept             {len(kept):>8,}")
    print(f"unique sequences       {len(representative):>8,}  -> BLAST query")
    print(f"Wrote {args.query_output}")
    print(f"Wrote {args.map_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
