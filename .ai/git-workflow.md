# Git, branches, and releases

## Git is controlled by the user

Inspect status, diffs, logs, and branches freely. Never run a state-changing operation unless the user asks for that operation by name: `git add`, `commit` (including `--amend`), `push`, `merge`, `rebase`, `reset`, `checkout -- <path>`, `restore`, `clean`, switching or deleting branches, creating or removing worktrees, tagging, or creating or editing issues, pull requests, or comments. Finishing a task is not a request to commit. Commit and push are separate permissions, and each approval covers only the stated operation.

Before an authorised commit, report the files in scope, the checks run and their results, and the proposed message, and run `git diff --check`. Before a push, report the branch, the commits, and the remote. Preserve unrelated tracked and untracked changes. Use `git mv` for tracked moves and renames.

## Commit messages

Before suggesting, writing, or reviewing a commit message, read `.ai/commit-messages.md` in the current turn and inspect the actual diff or staged changes.

## Branch model

- `master` holds released work.
- Work for the next release collects on a version branch named after the version (for example `1.3.0`).
- Branch features and fixes from the active version branch, never from another feature branch. Name branches after the problem (`add_chi2_method`, `issue60`).
- Feature branches return to the version branch through a pull request, squash-merged. When the release is ready, the version branch is merged into `master`.
- Delete merged branches only when the user asks.

## Releases

A release is prepared on its version branch:

1. Update `Version` in `DESCRIPTION` and the `# RRatepol <version>` section of `NEWS.md`.
2. Run `devtools::document()`, `devtools::test()`, and `devtools::check()`; for a CRAN submission also run `devtools::check(cran = TRUE)` and record the environments in `cran-comments.md`.
3. Render `README.qmd` and rebuild the website with `Rscript tools/build_pkgdown_site.R`; look at the result.
4. Merging, tagging, GitHub releases, and CRAN submission happen only on the user's explicit instruction.

Keep generated documentation, behaviour changes, and release metadata in separate commits unless the user combines them deliberately.

## Issue, pull-request, and comment text

Write for a reader who knows ecology and R but has never seen this work and has not read any plan, chat, or other issue. The text must explain itself; links support it and never replace it.

- Open with two or three plain sentences: what this is about, why it matters, and what the reader is asked to do.
- Give the background: which functions are affected, what happens now, and what is wrong or missing. Explain project-specific terms on first use.
- For an issue: the task and what "done" means, as checkable points.
- For a pull request: what changed and why, what was checked (command and result), what was not checked, and any decision left.
- For every reference to another issue or pull request, add a few words on what it is or decided. Write each relationship on its own list line (`- Closes #12`, `- Part of #10`).
- Never use plan codes, phase numbers, or session shorthand without a plain explanation.
- Do not hard-wrap the text; one paragraph or list item per line.
- Before posting, reread it as that newcomer and fill any gap.
