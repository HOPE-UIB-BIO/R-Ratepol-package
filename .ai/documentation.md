# Documentation and website

## Sources and generated files

| Edit this source | Generated output (never edit by hand) | How to regenerate |
| --- | --- | --- |
| roxygen blocks in `R/*.R` | `man/*.Rd`, `NAMESPACE` | `devtools::document()` |
| `README.qmd` | `README.md`, `man/figures/README-*` | `quarto render README.qmd` |
| `vignettes/*.qmd` (pkgdown articles) | `vignettes/*_files/`, pages in `docs/articles/` | the website build |
| `_pkgdown.yml`, sources above | `docs/` | `Rscript tools/build_pkgdown_site.R` |

- The articles in `vignettes/` (`package-description.qmd`, `workflow-example.qmd`) are Quarto documents rendered only for the pkgdown website. `.Rbuildignore` excludes them from the package build, so they are not CRAN vignettes and `devtools::check()` does not run them.
- `docs/` is committed and served by GitHub Pages at <https://hope-uib-bio.github.io/R-Ratepol-package/>.
- `tools/build_pkgdown_site.R` fixes image alt attributes in `README.md`, builds the site, and then removes pages that are not part of it: agent instructions (`AGENTS.md`, `CLAUDE.md`) and the instruction files in `.github/`. When an adapter or a root-level Markdown file is added, check that it does not appear in `docs/`.
- `NEWS.md` is the changelog for people who use the package. Each release starts with `# RRatepol <version>`. Describe what changed for users and keep internal history out.

## Writing

- Never hard-wrap Markdown prose; keep each paragraph and list item on one source line. In `.qmd` files, write one sentence per line. Do not reflow unrelated existing text.
- Use British spelling (`Language: en-GB`).
- Describe the current state of the package, not the history of a change.
- Keep code in README and articles reproducible from `example_data` or public sources (the workflow article downloads data from Neotoma). R chunks follow `.ai/r-coding.md`.
- Never invent results, figures, or citations. Numbers in articles come from code chunks.

## Checks

Render only the affected source during normal documentation work and inspect warnings, links, and output. Before a release or after site-wide changes, rebuild the whole website and look at the rendered pages. `devtools::check()` remains the check for roxygen examples.
