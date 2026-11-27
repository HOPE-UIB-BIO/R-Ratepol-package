# RRatepol: agent guide

RRatepol is an R package for estimating rate of change (RoC) from community data in time series, mainly fossil pollen records. It implements the method of Mottl et al. (2021, Review of Palaeobotany and Palynology 293, <https://doi.org/10.1016/j.revpalbo.2021.104483>).

This file is the entry point for coding assistants. The main guidance lives in `.ai/`; this file routes each task to it. The guidance is self-contained: it does not rely on any user-level instruction files, so collaborators and other tools get the same rules.

## Precedence

1. The user's current request.
2. This file and the files in `.ai/`.
3. The assistant's user-level instructions, for anything this repository does not define.
4. Defaults.

## Required reading

| Task | Read first |
| --- | --- |
| Any work in this repository | `AGENTS.md` |
| R code in `R/`, `tests/`, or `tools/`, and R chunks in `.qmd` files | `.ai/r-coding.md` |
| Creating, editing, documenting, or testing a function | `.ai/r-functions.md` (follow it in order), `.ai/r-coding.md` |
| README, articles in `vignettes/`, `NEWS.md`, pkgdown website | `.ai/documentation.md` |
| Debugging or bug fixes | `.ai/debugging.md` |
| Git, branches, issues, pull requests, releases | `.ai/git-workflow.md` |
| Suggesting, writing, or reviewing a commit message | `.ai/git-workflow.md`, then reread `.ai/commit-messages.md` in the current turn |
| Reviews, plans for large changes, delegating to other agents | `.ai/review-checklist.md`, `.ai/task-routing.md`, `.ai/agents/` |

## Key commands

Run from the package root in a fresh R session (for example `Rscript -e "devtools::test()"`). On Windows, `Rscript` exits with code 1 whenever it prints a warning, so judge success from the output text, not the exit code.

| Purpose | Command |
| --- | --- |
| Load the current source | `devtools::load_all()` |
| Regenerate `man/` and `NAMESPACE` | `devtools::document()` |
| Run one test file | `devtools::load_all()`, then `testthat::test_file("tests/testthat/test-<function_name>.R")` |
| Run the full test suite | `devtools::test()` |
| Run the opt-in parallel timing test (about 25 s) | `Sys.setenv(RRATEPOL_TEST_TIMING = "true")`, then `testthat::test_file("tests/testthat/test-estimate_roc.R")` |
| Run R CMD check | `devtools::check()` |
| Render the README | `quarto render README.qmd` |
| Build the website into `docs/` | `Rscript tools/build_pkgdown_site.R` |

## Rules that always apply

- **Numerical results belong to the user.** RRatepol produces published scientific results. A suspected bug that could change RoC values, peak points, smoothing, binning, or any other numerical output is confirmed with a test or a reproducible example and reported with the evidence; it is not fixed until the user decides. The same applies to changing defaults, seeds, or algorithms.
- **Keep the public interface.** Exported function names, arguments, defaults, returned columns, and error behaviour change only when the task explicitly asks for it, together with roxygen, tests, and `NEWS.md`. The deprecated `fc_*` functions stay as `lifecycle` wrappers.
- **Tests first.** Behaviour changes follow the test-first sequence in `.ai/r-functions.md`.
- **One function per file, no exceptions.** A touched file that defines several functions is split in the same change (see `.ai/r-coding.md`).
- **Generated files are never edited by hand:** `NAMESPACE`, `man/`, `README.md`, `docs/`.
- **No new dependency** without the user's approval. Keep the package installable from CRAN: no GitHub-only packages, and code must run on the declared `R (>= 4.0)`.
- **Git is controlled by the user.** Inspect status and diffs freely, but never commit, push, merge, rebase, switch branches, tag, or create pull requests or issues unless the user asks for that operation by name.
- **Text formatting.** The 80-character limit applies only to R code, R comments, and roxygen lines. Never hard-wrap Markdown prose; keep each paragraph and list item on one source line. In `.qmd` files, write one sentence per line. Use ASCII only in R source.

## Tool adapters

| Tool | Entry file |
| --- | --- |
| Codex and other `AGENTS.md`-aware tools | `AGENTS.md` |
| Claude Code | `CLAUDE.md` -> `AGENTS.md` |
| GitHub Copilot | `.github/copilot-instructions.md`, `.github/instructions/*.instructions.md` (by `applyTo`), `.github/commit-instructions.md` |

Adapters only point to `AGENTS.md` and `.ai/`; never put rules in an adapter. When guidance changes, edit `AGENTS.md` or `.ai/` first, then check that every adapter still points to an existing file.
