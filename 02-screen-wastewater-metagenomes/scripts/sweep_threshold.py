#!/usr/bin/env python3
"""Re-derive the threshold sweep from committed per-read evidence.

Reads `results/calibration/read_blast_deacon.tsv`, which carries one row per
candidate read with its diagnostic 31-mer count and the independent
local-alignment BLAST classification, and rewrites the threshold table the
choice rests on.

  threshold_read_counts.tsv   reads surviving each threshold, by class, and
                              runs with at least one surviving read

The table is regenerated from scratch, so running this against a clean checkout
reproduces the committed file exactly. No database, no cluster, no network.

A read counts toward a threshold when it carries at least that many diagnostic
31-mers on its own. That is the per-read rule used throughout this work, and it
is deliberately stricter than Deacon's paired mode, which pools distinct hits
across both mates.
"""

from __future__ import annotations

import argparse
import collections
import csv
from pathlib import Path

STAGE_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_EVIDENCE = STAGE_ROOT / "results/calibration/read_blast_deacon.tsv"
CLASSES = ("target", "non_target", "top_tie", "no_hit")


def load(evidence: Path) -> list[dict]:
    with evidence.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit(f"no rows in {evidence}")
    missing = {"read_id", "blast_class", "sample", "deacon_hits"} - set(rows[0])
    if missing:
        raise SystemExit(f"{evidence} is missing columns: {sorted(missing)}")
    return rows


def write_tsv(path: Path, header: tuple[str, ...], rows: list[tuple]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    parser.add_argument("--outdir", type=Path,
                        default=STAGE_ROOT / "results/calibration")
    parser.add_argument("--max-threshold", type=int, default=0,
                        help="highest threshold to report (default: the highest "
                             "k-mer count observed on any read)")
    args = parser.parse_args()

    rows = load(args.evidence)
    hits = {r["read_id"]: int(r["deacon_hits"]) for r in rows}
    top = args.max_threshold or max(hits.values())

    read_rows = []
    for t in range(1, top + 1):
        retained = [r for r in rows if int(r["deacon_hits"]) >= t]
        counts = collections.Counter(
            r["blast_class"] for r in retained)
        retained_runs = len({r["sample"] for r in retained})
        read_rows.append((t, *(counts.get(c, 0) for c in CLASSES), retained_runs))

    write_tsv(args.outdir / "threshold_read_counts.tsv",
              ("threshold", "target_reads", "non_target_reads",
               "top_tie_reads", "no_hit_reads", "runs_with_retained_reads"),
              read_rows)

    # The chosen threshold is the lowest supported by at least one target read
    # and no non-target, tied, or no-hit reads.
    chosen = next((t for t, target, nt, tie, no_hit, _ in read_rows
                   if target > 0 and nt == 0 and tie == 0 and no_hit == 0), None)
    worst_nt = max((int(r["deacon_hits"]) for r in rows
                    if r["blast_class"] == "non_target"), default=0)
    print(f"reads: {len(rows)}  runs: {len({r['sample'] for r in rows})}")
    print(f"highest diagnostic k-mer count on a non-target read: {worst_nt}")
    print(f"lowest supported threshold: {chosen}")
    if chosen:
        row = read_rows[chosen - 1]
        print(f"at that threshold, target reads retained: {row[1]}")
    print(f"wrote {args.outdir / 'threshold_read_counts.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
