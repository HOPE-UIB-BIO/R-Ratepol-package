# Functions, documentation, and tests

Read `.ai/r-coding.md` first; it holds the layout, naming, and package-specific rules that also apply here.

## Test-first sequence

Every new function and every behaviour change follows these steps in order:

1. **Specification.** Write or update the roxygen block first: purpose, arguments, return value, errors, and examples. For a new function, the body is a stub (`return(NULL)` or `stop("not implemented")`).
2. **Tests.** Write or update `tests/testthat/test-<function_name>.R` from the specification, not from the implementation.
3. **Confirm failure.** Run the test file and confirm that the new tests fail, for the intended reason and not because of broken setup. A test that passes against a stub checks nothing; revise it.
4. **Implement** the smallest coherent change until the tests pass. Keep the roxygen block in step with any change to arguments or behaviour.
5. **Document.** Run `devtools::document()` and inspect the changes to `NAMESPACE` and `man/`.
6. **Full suite.** Run `devtools::test()`.
7. **Check.** Run `devtools::check()`; it must pass with no new errors, warnings, or notes.

A task is complete only when steps 6 and 7 pass. Report which steps were run. If a pure documentation correction cannot fail first, say why and check the generated output instead. Never commit at the end of the sequence; stop and report (see `.ai/git-workflow.md`).

Focused run:

```r
devtools::load_all()

testthat::test_file(
  "tests/testthat/test-<function_name>.R"
)
```

## Writing the function

- One function per file, with a matching test file (see `.ai/r-coding.md`).
- Validate arguments at the start with `assertthat::assert_that()` (always with `msg`) and the `util_check_*()` helpers.
- Take `silent = FALSE` when the function prints anything, and forward it.
- Build the result under a name (`res` or `res_*`) and end with `return(res)`. Early `return(NULL)` and `return(invisible(res))` are fine.
- Before writing a new helper, search `R/` for an existing one that does nearly the same thing, and extend it with its tests where reasonable.

## Roxygen

Every function has a roxygen block directly above it, within 80 characters per line. Put tag text on the line after the tag. State argument types in backticks (`` `numeric` ``, `` `data.frame` ``, `` `NULL` ``), valid values, and what the default means.

```r
#' @title Estimate something per working unit
#' @description
#' One or two sentences on what the function computes.
#' @param data_source
#' `data.frame` with one row per sample and a column `age`.
#' @param silent
#' `logical`. If `TRUE`, suppress all console output.
#' @return
#' A `tibble` with columns `age` and `value`.
#' @details
#' Assumptions, units, and edge-case behaviour.
#' @seealso [estimate_roc()]
#' @export
#' @examples
#' \donttest{
#' data("example_data", package = "RRatepol")
#'
#' res <-
#'   estimate_something(
#'     data_source = example_data$sample_age[[1]]
#'   )
#' }
```

- `@title` is one short sentence in sentence case with no final full stop.
- Exported functions have `@export` and runnable `@examples` that use `example_data`. Internal functions have `@keywords internal` and no `@export`.
- New examples run directly when fast; wrap slow ones in `\donttest{}`, never `\dontrun{}`. The existing exported examples use `\dontrun{}`; changing them is a separate decision for the user.
- Deprecated functions use `lifecycle::badge("deprecated")` in the roxygen and `lifecycle` warnings in the body, as in the `fc_*` files.
- Cite methods with `@references`.

## Tests

Tests live in `tests/testthat/`, one file per function, `test-<function_name>.R`. `devtools::load_all()` or `tests/testthat.R` loads the package; never call `library()` or `source()` in a test file. Shared fixtures come from `tests/testthat/helper-*.R`, preferably built from `example_data`.

Style rules most often missed:

1. The 80-character limit applies to every R line, including `test_that()` descriptions. Shorten a long description; do not split it with `paste0()`.
2. Never use `$`; use `dplyr::pull(res, column)` or `purrr::chuck(list, "element")`.
3. After `<-`, put a function call on the next line for every object. Only scalar literals and `NULL` stay on the line.
4. Namespace every call, base R included: `testthat::test_that()`, `testthat::expect_equal()`, `base::nrow()`, `base::colnames()`, `base::length()`, `base::set.seed()`. Note `stats::sd()`, not `base::sd()`.
5. Use only packages in `Imports` or `Suggests`.
6. Expectations may pass `object` and `expected` by position (`testthat::expect_equal(base::nrow(res), 3L)`); this is the only exception to named arguments.
7. No project banner at the top of the file. Group tests into sections in the style already used in `tests/testthat/` and name each test `"function_name() <behaviour>"`.

What to test:

- **Intended behaviour,** inferred from the function name, the roxygen specification, and argument names and defaults, before the implementation. Where they disagree with the code, test the specification and report the disagreement.
- **Input checks:** a valid call; a wrong type or class for each argument; missing columns; invalid values (unsupported `method` strings, negative sizes); `NA`, `NULL`, empty vectors, and zero-row data where relevant. Match error messages with a regex.
- **Output structure:** class, column names, dimensions, and every output mode an argument switches.
- **Correctness:** small inputs whose expected result is computed by hand in the test, compared with `testthat::expect_equal()` (`tolerance = 1e-8` for floating point); invariants such as monotonic ages or proportions that sum to one.
- **Argument combinations:** each `method` value and logical switch, optional arguments given and omitted.
- **Messages:** output appears by default and `silent = TRUE` removes it; errors still occur with `silent = TRUE`.
- **Randomness:** wrap calls in `base::set.seed(900723)`.
