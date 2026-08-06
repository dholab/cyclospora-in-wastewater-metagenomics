#!/usr/bin/env python3
"""Build the GitHub Pages reading version of the manuscript from `README.md`.

The README is the single source of the manuscript text. This script renders it to
`docs/index.html` and makes the three changes the web version needs:

1. **The figure becomes interactive.** The static SVG in the README is replaced by
   the Vega-Lite specification behind it, rendered in place with vega-embed, so
   read counts, sequencing depth, and contributing runs appear on hover. The SVG
   stays as the `<noscript>` fallback.
2. **Repository links resolve.** Pages serves only `docs/`, so links to files in
   the repository are rewritten to `blob`/`tree` URLs on GitHub. Anchors and
   external links are left alone.
3. **The text is styled for reading**, in one column, light or dark to match the
   reader's system.

Run `pixi run build` in this directory. Regenerate and commit `index.html`
whenever the README changes.
"""

from __future__ import annotations

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

# The README line that places the static figure, replaced by the live chart.
FIGURE_IMAGE = re.compile(
    r'<p><img alt="[^"]*" src="[^"]*cyclospora_heatmap\.svg"\s*/?></p>'
)

FIGURE_BLOCK = """
<figure class="chart">
  <div id="heatmap" role="img"
       aria-label="Heatmap of Cyclospora cayetanensis diagnostic reads per billion, by sewershed and fortnight"></div>
  <noscript>
    <img src="assets/cyclospora_heatmap.svg"
         alt="Heatmap of Cyclospora cayetanensis diagnostic reads per billion, by sewershed and fortnight">
  </noscript>
  <p class="chart-note">Hover any cell for its read count, sequencing depth, and contributing runs.</p>
</figure>
"""


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

    # nl2br keeps the affiliation block one line per affiliation, as written. Every
    # paragraph in the README is a single source line, so nothing else is affected.
    converter = markdown.Markdown(
        extensions=["tables", "fenced_code", "attr_list", "sane_lists", "nl2br"]
    )
    body = converter.convert(text)
    title_html = markdown.Markdown().convert(title_md).removeprefix("<p>").removesuffix("</p>")

    body = FIGURE_IMAGE.sub(FIGURE_BLOCK, body, count=1)
    if "id=\"heatmap\"" not in body:
        raise SystemExit("the Figure 1 image was not found in the rendered README")
    body = rewrite_links(body)

    (DOCS / "assets").mkdir(exist_ok=True)
    shutil.copyfile(SVG, DOCS / "assets/cyclospora_heatmap.svg")
    spec = json.loads(SPEC.read_text())
    spec.pop("title", None)  # the figure caption in the text carries the title
    (DOCS / "assets/cyclospora_heatmap.vl.json").write_text(json.dumps(spec, separators=(",", ":")))

    template = (DOCS / "template.html").read_text()
    page = (
        template.replace("{{TITLE_TEXT}}", html.escape(re.sub(r"<[^>]+>", "", title_html)))
        .replace("{{TITLE_HTML}}", title_html)
        .replace("{{REPO}}", REPO)
        .replace("{{BODY}}", body)
    )
    (DOCS / "index.html").write_text(page)
    print(f"wrote {DOCS / 'index.html'} ({len(page):,} bytes)")


if __name__ == "__main__":
    build()
