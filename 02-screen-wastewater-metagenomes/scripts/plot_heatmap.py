#!/usr/bin/env python3
"""Render the sewershed-by-fortnight *Cyclospora* heatmap from the SRA screen.

Writes four files from `results/sra_sample_summary.tsv` and `results/casper_sites.tsv`:

  results/figures/cyclospora_heatmap.svg      static figure, embedded in the README
  results/figures/cyclospora_heatmap.vl.json  Vega-Lite spec, the portable artifact
  results/figures/cyclospora_heatmap.html     vega-embed wrapper around that spec
  results/site_fortnight_matrix_per_billion.tsv   the matrix behind all three

Every column is a sewershed named by its `public_code_CASPER_SRA` — the identity
under which it is deposited in SRA — so nothing here depends on an internal name.

Encoding decisions, in the order the choices were made.

* **The value is a depth-normalised rate**, distinct diagnostic reads per billion
  reads sequenced, so samples of very different depth are comparable. A cell
  covering more than one sample pools summed reads over summed depth. It is never
  a mean of rates, which would let a shallow sample dominate a deep one.
* **Reads are already deduplicated.** The value is Deacon's count of distinct
  diagnostic read sequences per sample, so PCR and optical copies are collapsed
  before anything is pooled. Raw retained counts are in the summary alongside.
* **Rows are fortnights.** Most sewersheds are sampled weekly or fortnightly, so
  a two-week bin is the coarsest that still shows a rise beginning and the finest
  that stays legible across the multi-year public record. Bins are counted from a
  fixed Monday so a cell always holds the same two calendar weeks.
* **Columns run west to east**, from the coordinates in `casper_sites.tsv`, so a
  regional pattern reads as one band rather than being scattered by rank.
* **Magnitude gets a sequential single-hue ramp**, faint to dark blue, binned
  rather than mapped linearly, because the rate distribution is strongly right
  skewed. A linear ramp would collapse most positives into the lightest step.
* **Four cell states stay visually distinct** so that "no sample" can never be
  misread as "looked and found nothing". Grey is no sample, white is screened
  with nothing found, and blue is abundance.

Standard library only. No network, no plotting dependency.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import html
import json
from collections import defaultdict
from pathlib import Path

STAGE_ROOT = Path(__file__).resolve().parent.parent

# A Monday, so fortnight bins are anchored to the calendar rather than to the
# window: a cell always holds the same two calendar weeks.
FORTNIGHT_EPOCH = dt.date(1970, 1, 5)

# Faint to dark blue. White is a screened zero, grey is an unscreened cell.
BINS = [(2.0, "0 to 2", "#cde2fb"), (5.0, "2 to 5", "#9ec5f4"),
        (10.0, "5 to 10", "#6da7ec"), (20.0, "10 to 20", "#3987e5"),
        (40.0, "20 to 40", "#1c5cab"), (float("inf"), "40 or more", "#0d366b")]
ZERO, NODATA = "#ffffff", "#d5d4cd"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e0"


def colour(rate: float) -> str:
    if rate <= 0:
        return ZERO
    return next(c for upper, _, c in BINS if rate < upper)


def fortnight_of(iso: str) -> str:
    d = dt.date.fromisoformat(iso)
    monday = d - dt.timedelta(days=d.weekday())
    blocks = (monday - FORTNIGHT_EPOCH).days // 14
    return (FORTNIGHT_EPOCH + dt.timedelta(days=14 * blocks)).isoformat()


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--summary", type=Path,
                   default=STAGE_ROOT / "results/sra_sample_summary.tsv")
    p.add_argument("--sites", type=Path,
                   default=STAGE_ROOT / "results/casper_sites.tsv")
    p.add_argument("--outdir", type=Path, default=STAGE_ROOT / "results")
    p.add_argument("--min-timepoints", type=int, default=10,
                   help="sewersheds sampled on fewer distinct dates are excluded")
    p.add_argument("--raw", action="store_true",
                   help="plot raw retained reads instead of the distinct count")
    args = p.parse_args()

    place: dict[str, tuple[float, float]] = {}
    role: dict[str, str] = {}
    with args.sites.open(newline="") as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            role[r["casper_code"]] = r["role"]
            if r["longitude"] and r["latitude"]:
                place[r["casper_code"]] = (float(r["longitude"]),
                                           float(r["latitude"]))

    value_field = "diagnostic_reads" if args.raw else "distinct_diagnostic_reads"
    with args.summary.open(newline="") as handle:
        rows = [r for r in csv.DictReader(handle, delimiter="\t")
                if r["collection_date"]
                and role.get(r["casper_code"]) == "surveillance"]

    dates: dict[str, set] = defaultdict(set)
    for r in rows:
        dates[r["casper_code"]].add(r["collection_date"])
    keep = {c for c, d in dates.items() if len(d) >= args.min_timepoints}
    kept = [r for r in rows if r["casper_code"] in keep]
    dropped = len(dates) - len(keep)

    # West to east; a code with no coordinates sorts last rather than at zero.
    sites = sorted(keep, key=lambda c: (place.get(c, (float("inf"), 0.0))[0],
                                        -place.get(c, (0.0, 0.0))[1], c))

    # Pool reads over depth per cell; never average the per-sample rates.
    hits: dict[tuple, int] = defaultdict(int)
    depth: dict[tuple, int] = defaultdict(int)
    samples: dict[tuple, list] = defaultdict(list)
    for r in kept:
        key = (r["casper_code"], fortnight_of(r["collection_date"]))
        hits[key] += int(r[value_field] or 0)
        depth[key] += int(r["input_reads"] or 0)
        samples[key].append(r["public_id"])
    fortnights = sorted({f for _, f in hits})
    rate = {k: (hits[k] / depth[k] * 1e9 if depth[k] else 0.0) for k in hits}

    args.outdir.mkdir(parents=True, exist_ok=True)
    matrix = args.outdir / "site_fortnight_matrix_per_billion.tsv"
    with matrix.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["fortnight_beginning", *sites])
        for f in fortnights:
            w.writerow([f] + [f"{rate[(s, f)]:.3f}" if (s, f) in rate else ""
                              for s in sites])

    # --- geometry -----------------------------------------------------------
    cell, gap = 18, 2
    pitch = cell + gap
    left, header = 108, 168
    grid_w, grid_h = len(sites) * pitch, len(fortnights) * pitch
    overhang = max(len(s) for s in sites) * 5.6 * 0.7071
    width = left + grid_w + overhang + 16
    height = header + grid_h + 108

    def cells_svg() -> str:
        out = []
        for yi, f in enumerate(fortnights):
            y = header + yi * pitch
            for xi, s in enumerate(sites):
                x = left + xi * pitch
                key = (s, f)
                fill = colour(rate[key]) if key in rate else NODATA
                out.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" '
                           f'rx="3" fill="{fill}"/>')
        return "\n".join(out)

    site_labels = "\n".join(
        f'<text x="{left + xi * pitch + cell / 2}" y="{header - 8}" fill="{MUTED}" '
        f'font-size="10" text-anchor="start" font-family="ui-monospace,monospace" '
        f'transform="rotate(-45 {left + xi * pitch + cell / 2} {header - 8})">'
        f'{html.escape(s)}</text>' for xi, s in enumerate(sites))

    row_labels = "\n".join(
        f'<text x="{left - 10}" y="{header + yi * pitch + cell - 5}" fill="{MUTED}" '
        f'font-size="10" text-anchor="end" font-family="ui-monospace,monospace">{f}</text>'
        for yi, f in enumerate(fortnights))

    key_items = [("No sample", NODATA), ("Screened, none found", ZERO)]
    key_items += [(lab, col) for _, lab, col in BINS]
    lx, ly = left, header + grid_h + 34
    legend = []
    for lab, col in key_items:
        legend.append(f'<rect x="{lx}" y="{ly}" width="12" height="12" rx="3" '
                      f'fill="{col}" stroke="{GRID}"/>')
        legend.append(f'<text x="{lx + 17}" y="{ly + 10}" fill="{MUTED}" font-size="10" '
                      f'font-family="ui-sans-serif,system-ui,sans-serif">{lab}</text>')
        lx += 24 + len(lab) * 5.6
    caption = [
        "Distinct diagnostic reads per billion reads sequenced. Columns are "
        "sewersheds named by their SRA code, ordered west to east.",
        f"{len(sites)} sewersheds with at least {args.min_timepoints} sampled "
        "timepoints in SRA. Rows are fortnights in which any of them was sampled.",
    ]
    for i, line in enumerate(caption):
        legend.append(f'<text x="{left}" y="{ly + 34 + i * 15}" fill="{MUTED}" '
                      f'font-size="10" '
                      f'font-family="ui-sans-serif,system-ui,sans-serif">{line}</text>')

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{height}" \
viewBox="0 0 {width:.0f} {height}" role="img" \
aria-label="Cyclospora distinct diagnostic reads per billion, by sewershed and fortnight">
<rect width="{width:.0f}" height="{height}" fill="#fcfcfb"/>
<text x="{left}" y="26" fill="{INK}" font-size="15" font-weight="600" \
font-family="ui-sans-serif,system-ui,sans-serif">\
Cyclospora cayetanensis in public wastewater sequencing, by sewershed and fortnight</text>
{site_labels}
{row_labels}
{cells_svg()}
{chr(10).join(legend)}
</svg>
"""
    (args.outdir / "figures").mkdir(parents=True, exist_ok=True)
    (args.outdir / "figures/cyclospora_heatmap.svg").write_text(svg)

    # --- interactive ----------------------------------------------------
    records = [{
        "site": s, "fortnight": f,
        "rate": round(rate[(s, f)], 3),
        "reads": hits[(s, f)], "depth": depth[(s, f)],
        "samples": len(samples[(s, f)]),
        "sample_ids": ", ".join(sorted(samples[(s, f)])),
    } for (s, f) in sorted(hits)]

    spec = {
        "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
        "title": {
            "text": "Cyclospora cayetanensis in public wastewater sequencing, "
                    "by sewershed and fortnight",
            "subtitle": [
                "Distinct diagnostic reads per billion reads sequenced. White is "
                "screened with nothing found; grey is no sample that fortnight.",
                "Columns are SRA codes, west to east; sewersheds sampled at "
                f"{args.min_timepoints} or more timepoints in SRA.",
            ],
            "anchor": "start", "fontSize": 15, "subtitleColor": MUTED,
        },
        "data": {"values": records},
        "mark": {"type": "rect", "stroke": "#ffffff", "strokeWidth": 1},
        "width": {"step": 18}, "height": {"step": 18},
        "encoding": {
            "x": {"field": "site", "type": "nominal", "sort": sites,
                  "axis": {"orient": "top", "labelAngle": -45, "title": None,
                           "labelColor": INK,
                           "labelFont": "ui-monospace,monospace",
                           "domain": False, "ticks": False}},
            "y": {"field": "fortnight", "type": "ordinal", "sort": fortnights,
                  "axis": {"title": "fortnight beginning", "labelColor": MUTED,
                           "domain": False, "ticks": False}},
            "color": {
                "field": "rate", "type": "quantitative",
                "scale": {"type": "threshold",
                          "domain": [0.0001] + [u for u, _, _ in BINS[:-1]],
                          "range": [ZERO] + [c for _, _, c in BINS]},
                "legend": {"title": "distinct reads per billion",
                           "titleColor": MUTED},
            },
            "tooltip": [
                {"field": "site", "title": "sewershed (SRA code)"},
                {"field": "fortnight", "title": "fortnight beginning"},
                {"field": "rate", "title": "distinct reads per billion",
                 "format": ".2f"},
                {"field": "reads", "title": "distinct diagnostic reads",
                 "format": ","},
                {"field": "depth", "title": "reads screened", "format": ","},
                {"field": "samples", "title": "runs pooled"},
                {"field": "sample_ids", "title": "run IDs"},
            ],
        },
        "config": {
            "view": {"fill": NODATA, "stroke": None},
            "background": "#fcfcfb",
            "font": "ui-sans-serif,system-ui,-apple-system,sans-serif",
        },
    }

    figures = args.outdir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    (figures / "cyclospora_heatmap.vl.json").write_text(
        json.dumps(spec, indent=1) + "\n")
    (figures / "cyclospora_heatmap.html").write_text(
        HTML.replace("__SPEC__", json.dumps(spec)))

    pos = sum(1 for r in kept if r["positive"] == "yes")
    print(f"sewersheds kept: {len(sites)} (excluded {dropped} with fewer than "
          f"{args.min_timepoints} timepoints)")
    print(f"SRA runs: {len(kept)}   positive: {pos}   fortnights: {len(fortnights)}")
    print(f"distinct diagnostic reads: {sum(hits.values()):,}   "
          f"reads screened: {sum(depth.values()):,}")
    print(f"wrote figures and matrix under {args.outdir}")
    return 0


# A thin vega-embed wrapper. The spec is inlined so the file opens from disk as
# well as over HTTP; only the Vega runtime comes from the CDN.
HTML = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Cyclospora in public wastewater sequencing, by sewershed and fortnight</title>
<script src="https://cdn.jsdelivr.net/npm/vega@5"></script>
<script src="https://cdn.jsdelivr.net/npm/vega-lite@5"></script>
<script src="https://cdn.jsdelivr.net/npm/vega-embed@6"></script>
<style>
:root { color-scheme: light; }
body { margin:0; padding:28px; background:#fcfcfb; color:#0b0b0b;
  font:14px/1.5 ui-sans-serif,system-ui,-apple-system,sans-serif; }
p.sub { margin:0 0 20px; color:#52514e; max-width:74ch; }
#chart { overflow:auto; }
</style></head><body>
<p class="sub">Hover any cell for the underlying read counts, sequencing depth,
and contributing runs. The static version of this figure, and the matrix behind
both, are in the same directory. Every column is a sewershed named by its SRA
code.</p>
<div id="chart"></div>
<script>vegaEmbed('#chart', __SPEC__, {actions: {editor: false}});</script>
</body></html>
"""


if __name__ == "__main__":
    raise SystemExit(main())
