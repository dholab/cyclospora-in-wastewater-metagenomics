#!/usr/bin/env python3
"""Build the GitHub Pages reading version of the manuscript from `README.md`.

The README is the single source of the manuscript text. This script renders it to
`docs/index.html` and makes the changes the web version needs:

1. **The figure becomes interactive.** The static SVG in the README is replaced by
   the Vega-Lite specification behind it, rendered in place with vega-embed, so
   read counts, sequencing depth, and contributing runs appear on hover. The SVG
   stays as the fallback if the chart cannot load.
2. **Authors move under the title**, where a reader expects them, rather than
   sitting in a section near the end.
3. **Display items are set off** from the text by rules, and every "Table 1",
   "Table 2", or "Figure 1" in the prose becomes a link that jumps to them.
4. **Repository links resolve.** Pages serves only `docs/`, so links to files in
   the repository are rewritten to `blob`/`tree` URLs on GitHub. Anchors and
   external links are left alone.

Run `pixi run build` in this directory. Regenerate and commit `index.html`
whenever the README changes.
"""

from __future__ import annotations

import gzip
import html
import json
import re
import shutil
from pathlib import Path

import markdown

DOCS = Path(__file__).resolve().parent
ROOT = DOCS.parent
REPO = "https://github.com/dholab/cyclospora-in-wastewater-metagenomics"
BRANCH = "main"

FIGURES = ROOT / "02-screen-wastewater-metagenomes/results/figures"
SPEC = FIGURES / "cyclospora_heatmap.vl.json"
SVG = FIGURES / "cyclospora_heatmap.svg"
READS = ROOT / "02-screen-wastewater-metagenomes/results/reads"
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# The README line that places the static figure, replaced by the live chart.
FIGURE_IMAGE = re.compile(
    r'<p><img alt="[^"]*" src="[^"]*cyclospora_heatmap\.svg"\s*/?></p>'
)

CHART = """<div class="chart-frame">
    <div id="figure-1-chart" role="img"
         aria-label="Heatmap of Cyclospora cayetanensis diagnostic reads per billion, by sewershed and fortnight"></div>
  </div>
  <p class="chart-hint">Hover any cell for its counts · click a cell with detections for its reads</p>
  <div id="reads-panel" hidden>
    <div class="reads-head">
      <div>
        <span class="reads-title" id="reads-title"></span>
        <span class="reads-sub" id="reads-sub"></span>
      </div>
      <div class="reads-acts">
        <button type="button" id="reads-copy">Copy FASTA</button>
        <button type="button" id="reads-save">Download</button>
        <button type="button" id="reads-close" aria-label="Close reads">Close</button>
      </div>
    </div>
    <pre id="reads-body"></pre>
  </div>"""


def collect_reads() -> dict[str, dict]:
    """Load the committed diagnostic reads, keyed by the sample id the chart uses.

    Reads are collapsed to distinct sequences within a run, which is the count the
    heatmap plots, with the number of observations kept alongside each one. The
    committed FASTAs hold every retained read, copies included; a cell that plots
    two distinct reads would otherwise open as sixteen near-identical records.

    Every file here is one public SRA run, so the payload the page ships carries
    nothing that is not already in the repository and in the BioProject.
    """
    payload: dict[str, dict] = {}
    for path in sorted(READS.glob("*.diagnostic_reads.fasta.gz")):
        sample, _, run = path.name.split(".", 1)[0].partition("__")
        distinct: dict[str, dict] = {}
        header = None
        for line in gzip.open(path, "rt"):
            line = line.strip()
            if line.startswith(">"):
                header = line[1:]
            elif line and header is not None:
                name, _, rest = header.partition(" ")
                kmers = re.search(r"diagnostic_kmers=(\d+)", rest)
                # A read and its reverse complement are one molecule, and that is
                # how the screen counts distinct reads, so collapse canonically.
                key = min(line, line.translate(COMPLEMENT)[::-1])
                seen = distinct.get(key)
                if seen is None:
                    distinct[key] = {
                        "id": name,
                        "k": int(kmers.group(1)) if kmers else None,
                        "s": line,
                        "n": 1,
                    }
                else:
                    seen["n"] += 1
                header = None
        records = list(distinct.values())
        if records:
            # A collection date can be deposited as more than one run, so the
            # entry accumulates rather than replacing what an earlier file left.
            entry = payload.setdefault(sample, {"runs": [], "reads": []})
            entry["runs"].append(run)
            entry["reads"].extend(records)
    if not payload:
        raise SystemExit(f"no diagnostic read files found under {READS}")
    return payload


def split_out_authors(text: str) -> tuple[str, str]:
    """Lift the authors and affiliations out of the body for the masthead."""
    match = re.search(
        r"^## Authors and affiliations\n\n(.+?)\n\n(.+?)\n\n(?=^## )", text, re.S | re.M
    )
    if not match:
        raise SystemExit("the Authors and affiliations section was not found")
    inline = markdown.Markdown(extensions=["nl2br"])
    byline = inline.convert(match.group(1)).removeprefix("<p>").removesuffix("</p>")
    inline.reset()
    affiliations = inline.convert(match.group(2)).removeprefix("<p>").removesuffix("</p>")
    block = (
        f'<p class="byline">{byline}</p>\n'
        f'  <p class="affiliations">{affiliations}</p>'
    )
    return text[: match.start()] + text[match.end():], block


def wrap_display_items(body: str) -> str:
    """Set tables and the figure in their own rules, with labeled captions."""

    def caption(text: str, kind: str, number: str) -> str:
        # "<strong>Table 1. What each filter removed.</strong> Each row counts…"
        text = re.sub(
            rf"^<strong>{kind} {number}\.\s*(.*?)</strong>",
            rf'<span class="fig-label">{kind} {number}</span><span class="fig-title">\1</span>',
            text,
            count=1,
        )
        return f"<figcaption>{text}</figcaption>"

    # Tables: the caption paragraph precedes its table in the source.
    def table_block(match: re.Match[str]) -> str:
        number, cap, table = match.group(2), match.group(1), match.group(3)
        return (
            f'<figure class="display" id="table-{number}">\n'
            f"  {caption(cap, 'Table', number)}\n"
            f'  <div class="table-scroll">{table}</div>\n'
            f"</figure>"
        )

    body = re.sub(
        r"<p>(<strong>Table (\d+)\..*?)</p>\s*(<table>.*?</table>)",
        table_block,
        body,
        flags=re.S,
    )

    # Figure: the caption paragraph follows the chart.
    def figure_block(match: re.Match[str]) -> str:
        number, cap = match.group(2), match.group(1)
        return (
            f'<figure class="display" id="figure-{number}">\n'
            f"  {CHART}\n"
            f"  {caption(cap, 'Figure', number)}\n"
            f"</figure>"
        )

    body, count = re.subn(
        r"@@CHART@@\s*<p>(<strong>Figure (\d+)\..*?)</p>",
        figure_block,
        body,
        flags=re.S,
    )
    if count != 1:
        raise SystemExit("the Figure 1 caption did not follow the figure")

    # Tables without a caption, such as Data availability, are not display items,
    # but they still need to scroll rather than push the page sideways.
    def scroll(match: re.Match[str]) -> str:
        table = match.group(0)
        return table if 'class="table-scroll"' in match.string[max(0, match.start() - 40): match.start()] \
            else f'<div class="table-scroll">{table}</div>'

    return re.sub(r"<table>.*?</table>", scroll, body, flags=re.S)


def link_callouts(body: str) -> str:
    """Turn every in-text mention of a display item into a link to it."""

    def replace(match: re.Match[str]) -> str:
        kind, number = match.group(1), match.group(2)
        return f'<a href="#{kind.lower()}-{number}">{kind} {number}</a>'

    # A mention opening a <strong> is the item's own caption, not a callout.
    return re.sub(r"(?<!<strong>)\b(Table|Figure) (\d+)\b(?!\.\s)", replace, body)


def rewrite_links(body: str) -> str:
    """Point repository-relative links at GitHub, leaving anchors and URLs alone."""

    def replace(match: re.Match[str]) -> str:
        attr, target = match.group(1), match.group(2)
        # Markdown obfuscates mailto: links into numeric character references, so
        # decode before deciding whether a target is already absolute.
        if html.unescape(target).startswith(
            ("http://", "https://", "#", "mailto:", "assets/")
        ):
            return match.group(0)
        path, _, fragment = target.partition("#")
        kind = "tree" if (ROOT / path).is_dir() else "blob"
        url = f"{REPO}/{kind}/{BRANCH}/{path.rstrip('/')}"
        if fragment:
            url = f"{url}#{fragment}"
        return f'{attr}="{url}"'

    return re.sub(r'(href|src)="([^"]+)"', replace, body)


def build() -> None:
    text = (ROOT / "README.md").read_text()

    # Strip the two lead-ins written for repository readers. One points here, which
    # is redundant on this page; the other is restated in the header below the title.
    text = re.sub(r"^\*\*The preferred way to read this manuscript.*?\*\*\n\n", "", text, flags=re.S | re.M)
    text = re.sub(r"^\*Note that supporting evidence.*?- DHO\*\n\n", "", text, flags=re.S | re.M)

    # The H1 becomes the page header rather than part of the flowing text.
    title_match = re.match(r"# (.+)\n", text)
    title_md = title_match.group(1) if title_match else "Manuscript"
    text = text[title_match.end():] if title_match else text

    text, authors = split_out_authors(text)

    converter = markdown.Markdown(
        extensions=["tables", "fenced_code", "attr_list", "sane_lists", "nl2br"]
    )
    body = converter.convert(text)
    title_html = markdown.Markdown().convert(title_md).removeprefix("<p>").removesuffix("</p>")

    body, count = FIGURE_IMAGE.subn("@@CHART@@", body, count=1)
    if count != 1:
        raise SystemExit("the Figure 1 image was not found in the rendered README")
    # Callouts first: while the captions are still plain <strong> lead-ins, the
    # lookbehind in link_callouts can tell a caption from a mention of it.
    body = link_callouts(body)
    body = wrap_display_items(body)
    body = rewrite_links(body)

    (DOCS / "assets").mkdir(exist_ok=True)
    shutil.copyfile(SVG, DOCS / "assets/cyclospora_heatmap.svg")
    spec = json.loads(SPEC.read_text())
    spec.pop("title", None)  # the figure caption in the text carries the title
    (DOCS / "assets/cyclospora_heatmap.vl.json").write_text(json.dumps(spec, separators=(",", ":")))
    reads = collect_reads()
    (DOCS / "assets/diagnostic_reads.json").write_text(json.dumps(reads, separators=(",", ":")))
    distinct = sum(len(entry["reads"]) for entry in reads.values())
    observed = sum(r["n"] for entry in reads.values() for r in entry["reads"])
    runs = sum(len(entry["runs"]) for entry in reads.values())
    print(
        f"  {distinct:,} distinct diagnostic reads ({observed:,} observations) "
        f"from {runs} runs available on click"
    )

    template = (DOCS / "template.html").read_text()
    page = (
        template.replace("{{TITLE_TEXT}}", html.escape(re.sub(r"<[^>]+>", "", title_html)))
        .replace("{{TITLE_HTML}}", title_html)
        .replace("{{AUTHORS}}", authors)
        .replace("{{REPO}}", REPO)
        .replace("{{BODY}}", body)
    )
    (DOCS / "index.html").write_text(page)
    print(f"wrote {DOCS / 'index.html'} ({len(page):,} bytes)")


if __name__ == "__main__":
    build()
