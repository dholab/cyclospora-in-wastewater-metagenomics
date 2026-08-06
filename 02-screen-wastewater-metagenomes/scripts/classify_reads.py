#!/usr/bin/env python3
"""Classify each candidate read against core-nt and assemble the sweep evidence.

Joins two independent facts for every read kept by `prepare_read_blast_query.py`:

  * `deacon_hits`  — diagnostic 31-mers the read carries on its own
  * `blast_class`  — what the read's full length says it is, from core-nt

and writes `results/calibration/read_blast_deacon.tsv`, the one table the
threshold sweep reads.

Four classes, from the best bit score on each side:

  target      the best alignment anywhere in core-nt is *Cyclospora*
  non_target  the best alignment is to something else
  top_tie     the two tie for the top bit score — the evidence does not decide
  no_hit      the read has no core-nt alignment at all

**Target is the genus, not the species.** GenBank carries genus-level deposits
such as U40261.1 "Cyclospora sp." whose 18S is identical to *C. cayetanensis*.
Under a species-strict rule a read matching both is scored a tie and looks like a
specificity failure, when the only thing that happened is that one depositor did
not name a species. That artifact moves the apparent clean threshold from 7 to 26
without a single genuinely non-Cyclospora read being involved. The genus is the
biologically meaningful unit for "is this Cyclospora", so it is the default;
`--target-taxids 88456` reproduces the species-strict variant for comparison.

Ties are held apart rather than assigned, because a tie is exactly the case where
the evidence does not decide.

Standard library only. Needs the BLAST output, not the database.
"""

from __future__ import annotations

import argparse
import collections
import csv
from pathlib import Path

STAGE = Path(__file__).resolve().parent.parent
# outfmt 6: qseqid qlen saccver staxids pident length mismatch gapopen qstart
#           qend sstart send evalue bitscore qcovhsp stitle
Q, STAXIDS, BITSCORE = 0, 3, 13
TIE_EPS = 0.1


def load_target_taxids(path: Path, override: str | None) -> set[str]:
    if override:
        return {t.strip() for t in override.replace(",", " ").split() if t.strip()}
    taxids = set()
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            taxids.add(line)
    if not taxids:
        raise SystemExit(f"no taxids in {path}")
    return taxids


def classify(blast: Path, targets: set[str]) -> dict[str, str]:
    """Best target vs best non-target bit score, per BLAST query."""
    best_target: dict[str, float] = collections.defaultdict(lambda: float("-inf"))
    best_other: dict[str, float] = collections.defaultdict(lambda: float("-inf"))
    seen: set[str] = set()
    with blast.open(newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) <= BITSCORE:
                continue
            query = row[Q]
            seen.add(query)
            try:
                bits = float(row[BITSCORE])
            except ValueError:
                continue
            taxa = {t for t in row[STAXIDS].replace(",", ";").split(";")
                    if t and t != "N/A"}
            if taxa & targets:
                best_target[query] = max(best_target[query], bits)
            # An unknown taxid counts as non-target: it is not evidence for us.
            if taxa - targets or not taxa:
                best_other[query] = max(best_other[query], bits)
    out = {}
    for query in seen:
        target, other = best_target[query], best_other[query]
        if target == float("-inf") and other == float("-inf"):
            out[query] = "no_hit"
        elif target == float("-inf"):
            out[query] = "non_target"
        elif other == float("-inf"):
            out[query] = "target"
        elif abs(target - other) < TIE_EPS:
            out[query] = "top_tie"
        else:
            out[query] = "target" if target > other else "non_target"
    return out


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    calibration = STAGE / "results/calibration"
    parser.add_argument("--blast", type=Path, required=True,
                        help="core-nt BLAST output, outfmt 6 as in the README")
    parser.add_argument("--map", type=Path,
                        default=calibration / "read_blast_map.tsv",
                        help="read -> BLAST representative, from "
                             "prepare_read_blast_query.py")
    parser.add_argument("--taxids", type=Path,
                        default=calibration / "cyclospora_genus_taxids.txt")
    parser.add_argument("--target-taxids", default=None,
                        help="override the target set, e.g. '88456' for the "
                             "species-strict variant")
    parser.add_argument("--output", type=Path,
                        default=calibration / "read_blast_deacon.tsv")
    args = parser.parse_args()

    targets = load_target_taxids(args.taxids, args.target_taxids)
    class_of = classify(args.blast, targets)

    rows = list(csv.DictReader(args.map.open(newline=""), delimiter="\t"))
    if not rows:
        raise SystemExit(f"no rows in {args.map}")

    out_rows = []
    for row in rows:
        # A representative with no row in the BLAST output had no hit at all.
        out_rows.append((row["read_id"], row["sample"], row["deacon_hits"],
                         class_of.get(row["representative_id"], "no_hit")))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("read_id", "sample", "deacon_hits", "blast_class"))
        writer.writerows(out_rows)

    counts = collections.Counter(row[3] for row in out_rows)
    print(f"target taxids={len(targets)} reads={len(out_rows):,}")
    print(f"class distribution: {dict(counts)}")
    # The number the threshold rests on: nothing that is not Cyclospora may
    # reach it.
    ceiling = max((int(k) for _, _, k, c in out_rows
                   if c in ("non_target", "top_tie")), default=0)
    print(f"highest k-mer count on a non-target or tied read: {ceiling}")
    print(f"lowest fully specific threshold: {ceiling + 1}")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
