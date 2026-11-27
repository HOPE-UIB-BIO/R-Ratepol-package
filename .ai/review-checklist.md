# Review checklist

Report findings first, ordered by severity, with file and line references. If nothing is found, say so and name what was not checked.

Check:

- **Numerical results.** Does the change alter RoC values, peak points, or other outputs for `example_data`? Unintended changes are findings; intended ones need the user's decision, tests, and a `NEWS.md` entry.
- **Public interface.** Exported names, arguments, defaults, returned columns, and error behaviour are unchanged unless the task asked for it; deprecated `fc_*` wrappers still work.
- **One function per file** in every touched file under `R/` and `tests/testthat/`, helpers included, each function with its matching test file. A file with several functions is a finding, never a deferral.
- **Tests** check the intended behaviour from the specification, cover invalid inputs and `silent`, and would fail without the change.
- **Package rules** from `.ai/r-coding.md`: `%>%` rather than `|>` or `\(x)` (R >= 4.0), no new dependencies, `assertthat` and `util_check_*()` checks with messages, output only through `util_output_*()` under `silent`.
- **Style** of new and rewritten code: assignment layout, multi-line conditions, named arguments, namespaces, 80 characters, ASCII only.
- **Documentation in step:** roxygen, `man/`, `NAMESPACE`, `NEWS.md`, README, and articles.
- **Generated files** changed only by their generator.
