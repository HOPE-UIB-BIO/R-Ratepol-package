# Commit messages

Read this file in the current turn before suggesting, writing, or reviewing a commit message, and inspect the actual diff or staged changes first. Do not infer a message from the task description alone.

When asked only for a commit message, return exactly one plain-text line: no body, bullets, quotes, code fences, labels, or explanation.

## Format

```text
<subject>: <short summary>
```

- One line, no body, no full stop at the end, at most 72 characters.
- The subject says what changed; the summary says what was done.
- No issue or pull-request numbers (`#115`, `issue 60`); issue relationships belong in the pull-request description.
- No temporary labels such as `phase 1` or `stage 3`; describe the lasting change.
- No `Co-Authored-By`, "Generated with", or other trailers.

## Subject

Use the narrowest meaningful subject:

| Change | Subject | Example |
| --- | --- | --- |
| One function | its name with parentheses | `estimate_roc(): add chi-squared dissimilarity method` |
| Several functions around one topic | a plain topic label | `Peak detection: switch default method to GAM-based approach` |
| Tests only | `tests` | `tests: add zero-row community matrix test for prepare_data()` |
| Roxygen, README, articles, website | `docs` | `docs: update package-description article for clarity` |
| Changelog | `NEWS` | `NEWS: update to 1.3.0` |
| Package metadata | `DESCRIPTION` | `DESCRIPTION: add gratia to Imports` |
| Website configuration or build script | `pkgdown` or the function name | `build_pkgdown_site(): drop agent instructions from the site` |
| Continuous integration | `ci` | `ci: update R-CMD-check workflow to r-lib/actions v2` |
| Agent guidance | `agents` | `agents: move guidance to AGENTS.md and .ai` |
| Editor configuration | `vscode` | `vscode: allow Rscript in chat terminal` |
| Other housekeeping | `chore` or `refactor` | `chore: add .ai folder to .Rbuildignore` |

## Wording

Start the summary with a specific verb: `add`, `adjust`, `correct`, `document`, `remove`, `replace`, `split`, `switch`, `update`.

Never use `enhance`, `feat`, `feature`, or `fix`; state what was done instead.
