#!/usr/bin/env python3
"""Re-count every published diagnostic read against the bait FASTA.

The reads in `results/reads/` are the evidence behind every positive cell in the
heatmap. Each file is one SRA run, named `<CODE>_<YYYYMMDD>__<accession>`. This
script recomputes, from scratch, the number of distinct diagnostic 31-mers each
read carries, and checks three things:

  1. every read reaches the threshold of 20 on its own, so nothing survives on a
     mate's hits (Deacon pools k-mer hits across mates in paired mode, which is
     why the counts are recomputed here rather than trusted);
  2. the count recorded in each read's header matches the recount;
  3. the per-run read totals match the `diagnostic_reads` column of
     `results/sra_sample_summary.tsv`.

Nothing is taken on faith from the screening run. The only inputs are the bait
FASTA from stage 01, the published reads, and the summary table.

Standard library only. No network, no environment.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import re
import sys
from pathlib import Path

STAGE_ROOT = Path(__file__).resolve().parent.parent
BAITS = (STAGE_ROOT.parent / "01-identify-cyclospora-specific-kmers/baits"
         / "cyclospora_cayetanensis_rrna_core_nt_validated_baits.fasta")
K = 31
COMPLEMENT = str.maketrans("ACGT", "TGCA")
SUFFIX = ".diagnostic_reads.fasta.gz"


def canonical(kmer: str) -> str:
    rev = kmer.translate(COMPLEMENT)[::-1]
    return kmer if kmer < rev else rev


def read_fasta(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    header, parts = None, []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header, parts = line[1:], []
            elif line:
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)


def count_diagnostic(sequence: str, baits: set[str]) -> int:
    """Distinct baits found in the read, counting both strands."""
    seen = set()
    upper = sequence.upper()
    for i in range(len(upper) - K + 1):
        kmer = upper[i:i + K]
        if "N" in kmer:
            continue
        c = canonical(kmer)
        if c in baits:
            seen.add(c)
    return len(seen)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--baits", type=Path, default=BAITS)
    p.add_argument("--reads", type=Path, default=STAGE_ROOT / "results/reads")
    p.add_argument("--summary", type=Path,
                   default=STAGE_ROOT / "results/sra_sample_summary.tsv")
    p.add_argument("--threshold", type=int, default=20)
    args = p.parse_args()

    baits = {canonical(s.upper()) for _, s in read_fasta(args.baits)}
    print(f"baits loaded: {len(baits):,} distinct canonical 31-mers")

    per_sample = collections.Counter()
    counts, below, mismatched = [], [], []
    for path in sorted(args.reads.glob("*" + SUFFIX)):
        sample = path.name[: -len(SUFFIX)]
        for header, sequence in read_fasta(path):
            n = count_diagnostic(sequence, baits)
            counts.append(n)
            per_sample[sample] += 1
            if n < args.threshold:
                below.append((sample, header.split()[0], n))
            claimed = re.search(r"diagnostic_kmers=(\d+)", header)
            if claimed and int(claimed.group(1)) != n:
                mismatched.append((header.split()[0], int(claimed.group(1)), n))

    if not counts:
        print("no reads found", file=sys.stderr)
        return 1
    counts.sort()
    print(f"reads recounted: {len(counts):,}")
    print(f"  diagnostic 31-mers per read: min {counts[0]}, "
          f"median {counts[len(counts) // 2]}, max {counts[-1]}")

    ok = True
    print(f"\n1. every read reaches {args.threshold} on its own")
    if below:
        ok = False
        print(f"   FAIL: {len(below)} read(s) below threshold")
        for s, r, n in below[:10]:
            print(f"     {s} {r} has {n}")
    else:
        print(f"   PASS: 0 of {len(counts):,} reads fall below")

    print("\n2. header counts match the recount")
    if mismatched:
        ok = False
        print(f"   FAIL: {len(mismatched)} header(s) disagree")
        for r, c, n in mismatched[:10]:
            print(f"     {r}: header says {c}, recount says {n}")
    else:
        print(f"   PASS: all {len(counts):,} headers agree")

    print("\n3. per-run totals match the summary table")
    # The reads are one file per run, named <public_id>__<accession>; the summary
    # carries those two fields separately, so the key is rebuilt to join them.
    with args.summary.open(newline="") as handle:
        expected = {f"{r['public_id']}__{r['sra_accession']}":
                    int(r["diagnostic_reads"] or 0)
                    for r in csv.DictReader(handle, delimiter="\t")}
    bad = []
    for sample, n in sorted(per_sample.items()):
        if expected.get(sample) != n:
            bad.append((sample, expected.get(sample), n))
    # every run the summary calls positive must have published reads
    missing = [s for s, n in expected.items() if n > 0 and s not in per_sample]
    if bad or missing:
        ok = False
        print(f"   FAIL: {len(bad)} disagreement(s), {len(missing)} run(s) "
              f"positive in the summary with no published reads")
        for s, e, g in bad[:10]:
            print(f"     {s}: summary says {e}, reads directory has {g}")
        for s in missing[:10]:
            print(f"     {s}: summary says {expected[s]}, no reads published")
    else:
        print(f"   PASS: {len(per_sample)} runs with reads, all totals agree")

    print("\n" + ("ALL CHECKS PASSED" if ok else "CHECKS FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
