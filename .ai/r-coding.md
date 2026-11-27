# RRatepol R coding

Main R guidance for package code in `R/`, tests in `tests/testthat/`, R chunks in `README.qmd` and `vignettes/*.qmd`, and scripts in `tools/`. For creating, documenting, and testing functions, also read `.ai/r-functions.md`.

## New code and legacy code

Most of the package predates these conventions. It uses `%>%`, `$`, `paste()`, `lapply()`, bare base R calls, and single-line function signatures.

- Apply the conventions below to new code and to code that the task materially rewrites.
- Leave untouched lines in their existing style. Do not combine a behaviour change with restyling unrelated code, and do not propose a package-wide restyle.
- Keep one function or pipeline internally consistent; do not mix styles inside it.
- The one-function-per-file rule is the exception: it applies to every file a change touches (see "Files and functions").

## Package-specific rules

These override the general conventions further down.

- **Pipe.** Use `%>%`, which the package re-exports in `R/utils-pipe.R`. The native pipe `|>` and the `\(x)` lambda need R 4.1, but `DESCRIPTION` declares `R (>= 4.0)`. Use them only after the user decides to raise the R requirement.
- **Dependencies.** Imports are limited to those in `DESCRIPTION`. `stringr`, `readr`, `cli`, and `glue` are not imported, so `paste()` and `paste0()` stay allowed in package code; do not add a package to avoid them. `dplyr::join_by()`, `dplyr::reframe()`, and `purrr::list_rbind()` need `dplyr (>= 1.1.0)` and `purrr (>= 1.0.0)`, which `DESCRIPTION` does not declare; ask before relying on them. `RUtilpol` was removed in 1.3.0 and must not come back.
- **Argument checks.** Use `assertthat::assert_that()` with an explicit `msg`, and the internal helpers `util_check_class()`, `util_check_col_names()`, `util_check_vector_values()`, and `util_check_if_integer()` (in `R/util_internal.R`). Do not use plain `stop()` or `stopifnot()` in new code.
- **Console output.** Print only through `util_output_comment()`, `util_output_heading()`, and `util_output_warning()`, wrapped in `if (isFALSE(silent))`. Every function that prints takes `silent = FALSE` and forwards it to the helpers it calls. `silent` never suppresses errors. Several functions (`estimate_roc()` and the internal functions it calls) also take `verbose`, which switches detailed progress messages on; `silent = TRUE` overrides it. Keep both arguments and their names; do not merge or rename them.
- **Internal functions.** They have no `@export` and carry `@keywords internal`. Existing helpers use the `util_` prefix without a leading dot; follow that pattern for new helpers of the same kind.
- **Parallel iteration.** `estimate_roc()` uses `pbapply::pblapply()` with a `parallel` cluster and progress-bar options. Follow that pattern for new parallel work, and only for CPU-heavy, independent iterations.
- **Randomness.** New tests and examples use `set.seed(900723)`. Existing seeds (for example `set.seed(123)` in `tests/testthat/helper-fixtures.R`) stay unless the user decides otherwise, because changing them changes expected values.
- **Plots.** `plot_roc()` returns a `ggplot2` object with `ggplot2::theme_classic()` and `ggplot2::labs()`; it does not save files. The package does not use `ggview`.
- **Spelling.** `DESCRIPTION` sets `Language: en-GB`. Use British spelling in documentation and messages (standardise, colour). Keep existing argument names that use other spellings.

## Files and functions

One function per file. This rule has no exceptions.

- Every function, exported or internal, lives alone in `R/<function_name>.R`, and the file name equals the function name.
- Every function has a matching `tests/testthat/test-<function_name>.R` that tests that function directly.
- Test helpers follow the same rule: one function per `tests/testthat/helper-<function_name>.R`.
- Do not define named functions inside another function.
- When a change edits, renames, or moves a function that shares its file with other functions, split that file in the same change, together with matching test files. Never defer the split to a later task, and never offer keeping functions together as an option in a plan.
- Files that currently define several functions: `R/util_internal.R` (nine `util_*` helpers) and `tests/testthat/helper-fixtures.R`. List them before proposing a commit with `grep -cE '^[A-Za-z_.][A-Za-z0-9_.]*[[:space:]]*<-[[:space:]]*function' R/*.R R/*.r tests/testthat/*.R | grep -v ':[01]$'`. A change is not ready while a file it touches is on that list.
- `DESCRIPTION` lists test files under `Config/testthat/start-first`; keep that list in step when a test file is renamed or split.

## Naming

Use lower `snake_case` and full descriptive words for objects, arguments, functions, and columns. No dots in names.

Object prefixes, when the type is stable:

| Prefix | Use |
| --- | --- |
| `data_` | data frame or tibble |
| `table_` | summary or reporting table |
| `list_` | list |
| `vec_` | vector |
| `mat_` | matrix |
| `mod_` | model |
| `plot_` | plot object |
| `path_` | file path |
| `flag_` | logical control variable |
| `sel_` | the selected item, such as the current loop value |
| `n_`, `index_` | counts, indices |
| `res_` | only the object a function returns |

Function names are verbs describing the returned value or side effect (`estimate_roc()`, `detect_peak_points()`, `subset_samples()`). Avoid vague `get_` and `make_` for new functions, and vague variants such as `_simple`, `_new`, or `_final`. Exported names and argument names (for example `data_source_community`, `smooth_method`) are part of the public interface; do not rename them.

Give each transformation stage a new name instead of overwriting the earlier object:

```r
data_community_filtered <-
  data_community_raw %>%
  dplyr::filter(sample_id %in% vec_sample_ids)
```

## Layout

- Assignment is `<-`, never `=` or `->`. When the right-hand side is a function call, put it on the next line, indented two spaces. Only a scalar literal or `NULL` stays on the assignment line.
- Function definitions keep `<-` on the same line as `function`. With two or more arguments, put each on its own line indented two spaces and close with `) {` on its own line. Do not align arguments under the opening bracket.
- Name arguments in calls (exception: `...`). With two or more arguments, put each on its own line.
- `if`, `for`, and `while` conditions always go on their own indented line, with `) {` on its own line, even when short. Keep `} else {` together.
- Indent with two spaces. No semicolons. `TRUE` and `FALSE`, never `T` and `F`.
- A space after commas and around infix operators; none around `::`, `$`, `[`, `[[`, `^`, `:`, or unary `-`.
- No R line over 80 characters, including comments, roxygen, and `test_that()` descriptions.
- Separate top-level statements in a block with one blank line; no blank lines inside one continued expression.

```r
estimate_example <- function(
  data_source,
  silent = FALSE
) {
  util_check_class(
    x = data_source,
    sel_class = "data.frame"
  )

  if (
    isFALSE(silent)
  ) {
    util_output_comment(msg = "Estimating example values")
  }

  res_example <-
    data_source %>%
    dplyr::mutate(
      value_scaled = value / max(value)
    )

  return(res_example)
}
```

## Comments

- Start with a capital letter and one space after `#`. Explain why a non-standard approach was chosen; do not narrate obvious lines.
- A multi-line comment continues on an indented line (`#    continued`).
- ASCII only: write `->`, never an arrow symbol, and no typographic quotes or dashes.
- Comments describe the current code; never add `# NEW:`, `# FIX:`, issue or pull-request numbers.
- Long functions may use section headers ending in `-----`, matching the style already in the file.

## Namespaces and tidyverse

- Use `pkg::function()` for every call to a non-base package. Base R calls may stay bare in package code; tests namespace every call, base R included (see `.ai/r-functions.md`). Never call `library()` or `require()` in package code.
- In new code, avoid `$` for element access: use `dplyr::pull()` for columns and `purrr::chuck()` for list elements.
- In new code, avoid the `apply` family; use `purrr::map()`, `purrr::map_dbl()`, `purrr::map2()`, and so on. Pipe the input into `purrr::map()`, with `rlang::set_names()` first when names should carry over.
- Anonymous functions in `purrr` use the tilde form (`~ .x`, `..1`). With nested maps, capture the outer `.x` in a named variable first.
- Use `purrr::walk()` for side effects when the return value is not needed.
- Do not use `base::attr()` for metadata, `eval(parse())`, `get()` inside data masks, `dplyr::transmute()`, or `purrr::map_dfr()`.
- Follow `dplyr::group_by()` with `dplyr::ungroup()` or `.groups = "drop"` unless grouped output is intended.
- Data masking: forward bare column arguments with `{{ }}`; access a column named by a string with `.data[[col_name]]`.

## Performance

Profile before optimising (`profvis::profvis()`, `bench::mark()`). Never grow an object in a loop; preallocate or use `purrr::map()`. Parallelise only CPU-heavy, independent work where the computation clearly outweighs the overhead.
