# The reading version of the manuscript

This directory is the source of the GitHub Pages site at
<https://dholab.github.io/cyclospora-in-wastewater-metagenomics/>, the preferred way to read this
work. It is the same manuscript as the [repository README](../README.md), with Figure 1 rendered as
an interactive chart and every method, script, and data file linked to its place in the repository.

| File | What it is |
|---|---|
| [`build_site.py`](build_site.py) | Generates `index.html` from `../README.md`. Swaps the static figure for the Vega-Lite chart, renders the front matter as the masthead, sets the tables and figure off as captioned display items, turns every in-text mention of one into a link to it, and rewrites repository-relative links to GitHub URLs. |
| [`template.html`](template.html) | The page shell: house palette, typography, and the vega-embed call. |
| `index.html` | The generated page. Not committed. [`.github/workflows/pages.yml`](../.github/workflows/pages.yml) renders and publishes it on every push to `main`, so an edit made from any machine, or in the GitHub web editor, reaches the site without anyone running the build locally. |
| `assets/` | Also generated, and not committed. The Vega-Lite specification and the static SVG copied from [`02-screen-wastewater-metagenomes/results/figures/`](../02-screen-wastewater-metagenomes/results/figures/), plus the distinct diagnostic reads assembled from [`results/reads/`](../02-screen-wastewater-metagenomes/results/reads/) for the click-to-read panel. |
| `pixi.toml`, `pixi.lock` | Python 3.12 and Python-Markdown, pinned. |
| `.nojekyll` | Serves the generated HTML as it is, without Jekyll processing. |

The README is the single source of the manuscript text, and it is ordered so this build has as
little to rearrange as possible. Everything above the first `##` heading is front matter, which
becomes the masthead: the title, then the byline, then the affiliations. A front-matter paragraph
that should stay in the repository but not appear on the site is marked with an HTML comment, which
GitHub renders as nothing:

```markdown
<!-- repo-only -->
**The preferred way to read this manuscript is the interactive version…**
```

Editing the README and pushing is enough to update the site. To preview the page before pushing,
build it and serve this directory:

```bash
cd docs
pixi run build
python3 -m http.server 8000   # then open http://localhost:8000/
```

The chart needs `vega`, `vega-lite`, and `vega-embed` from jsDelivr. Without network access, or with
JavaScript disabled, the page falls back to the static SVG.
