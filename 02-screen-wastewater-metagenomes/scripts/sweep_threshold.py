#!/usr/bin/env python3
"""Re-derive the threshold sweep from committed per-read evidence.

Reads `results/calibration/read_blast_deacon.tsv`, which carries one row per
candidate read with its diagnostic 31-mer count and the independent whole-read
BLAST classification, and rewrites the two sweep tables the threshold choice
rests on.

  threshold_blast_read_counts.tsv   reads surviving each threshold, by class
  absolute_threshold_curve.tsv      samples and pairs surviving each threshold

Both are regenerated from scratch, so running this against a clean checkout
reproduces the committed files exactly. No database, no cluster, no network.

A read counts toward a threshold when it carries at least that many diagnostic
31-mers on its own. A pair counts when either of its mates does. That is the
per-read rule used throughout this work, and it is deliberately stricter than
Deacon's paired mode, which pools distinct hits across both mates.
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

    # A pair is keyed by the read id with its mate suffix removed.
    pair_best: dict[str, int] = {}
    pair_sample: dict[str, str] = {}
    for r in rows:
        pair = r["read_id"].rsplit("/", 1)[0]
        pair_best[pair] = max(pair_best.get(pair, 0), int(r["deacon_hits"]))
        pair_sample[pair] = r["sample"]

    read_rows, curve_rows = [], []
    for t in range(1, top + 1):
        counts = collections.Counter(
            r["blast_class"] for r in rows if int(r["deacon_hits"]) >= t)
        read_rows.append((t, *(counts.get(c, 0) for c in CLASSES)))

        kept = [p for p, best in pair_best.items() if best >= t]
        curve_rows.append((t, len({pair_sample[p] for p in kept}), len(kept)))

    write_tsv(args.outdir / "threshold_blast_read_counts.tsv",
              ("threshold", "target_reads", "non_target_reads",
               "top_tie_reads", "no_hit_reads"), read_rows)
    write_tsv(args.outdir / "absolute_threshold_curve.tsv",
              ("abs_threshold", "retained_samples", "retained_pairs"), curve_rows)

    # The chosen threshold is the lowest at which neither a confidently
    # non-target read nor an ambiguous one survives.
    chosen = next((t for t, _, nt, tie, _ in read_rows if nt == 0 and tie == 0), None)
    worst_nt = max((int(r["deacon_hits"]) for r in rows
                    if r["blast_class"] == "non_target"), default=0)
    print(f"reads: {len(rows)}  pairs: {len(pair_best)}  "
          f"samples: {len(set(pair_sample.values()))}")
    print(f"highest diagnostic k-mer count on a non-target read: {worst_nt}")
    print(f"lowest threshold with no non-target and no tie reads: {chosen}")
    if chosen:
        row = read_rows[chosen - 1]
        print(f"at that threshold, target reads retained: {row[1]}")
    print(f"wrote 2 tables to {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
