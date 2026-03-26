# RRatepol 1.3.0

## Complete unit test suite

- `tests/testthat/` created from scratch (~3,000 lines across 19 test files, one per function)
- shared `helper-fixtures.R` added for reusable test data
- every function now has tests for valid input, invalid input with informative errors, and function-specific edge cases

## CI and test coverage

- new `test-coverage.yaml` workflow uploads coverage results to Codecov on every push
- `R-CMD-check.yaml` updated to current `r-lib/actions` conventions
- `testthat` (≥ 3.0.0) added to `Suggests`; test execution order declared via `Config/testthat/start-first`

## Removed runtime dependencies

- `RUtilpol` (GitHub-only) and `usethis` removed from `Imports`
- `Remotes:` field removed; package is now installable from CRAN without any GitHub dependency

## New internal utilities

- `R/util_internal.R`: `util_check_class()`, `util_check_col_names()`, `util_check_vector_values()`, `util_check_if_integer()`, `util_output_comment()`, `util_output_heading()`, `util_output_warning()`, `util_flatten_list_by_one()`
- `R/util_search_parameter.R`: Grimm-smoothing window-growth logic extracted into `util_search_parameter()`

## New `silent` argument

- `silent = FALSE` added to `estimate_roc()`, `run_iteration()`, `extract_data()`, `transform_into_proportions()`, `plot_roc()`, and all functions that produce console output
- when `TRUE`, suppresses all messages and warnings without affecting `verbose`

## Input validation — new assertions

- `estimate_roc()`: `standardise`, `tranform_to_proportions`, `verbose` must each be a single `TRUE`/`FALSE`; `interest_threshold` must be a single value when supplied; `use_parallel`, when numeric, must not be `0` or `NA` (closes #95)
- `plot_roc()`: `peaks` must be a single non-`NA` value; warning raised when `ROC` column contains `NA`s (closes #94)
- `detect_sni()`: columns 2 and 3 of `CharData` must not contain `NA`s (closes #93)
- `make_trend()`: `sel_method` must be a single value; return value always coerced to numeric vector, fixing `"non_linear"` returning an array (closes #92, closes #78)
- `transform_into_proportions()`: input must not be an empty data frame; `sel_method` must be a single value (closes #91)
- `subset_community()`: input must not be an empty data frame (closes #89)
- `reduce_data_simple()`: input must not be empty; `check_taxa` and `check_levels` must each be a single `TRUE`/`FALSE` (closes #88)
- `subset_samples()`: `bin_selection` must be `"first"`, `"random"`, or `NULL` (closes #87)
- `run_iteration()`: `standardise`, `tranform_to_proportions`, `verbose` must each be a single `TRUE`/`FALSE`; `time_standardisation` must not be `0` or `NA`; standardisation-failure error is now unconditional, not gated on `verbose` (closes #86, closes #62)
- `make_bins()`: `working_units` must be a single value (closes #85)
- `prepare_data()`: community and age inputs must not be `NULL` or empty; `bin_size` assertion runs only when `working_units != "levels"`, fixing spurious error with `working_units = "levels"` and `bin_size = NULL` (closes #84, closes #79)
- `reduce_data()`: community and age inputs must not be `NULL`; sample alignment uses `intersect()` across community, age, and age_un so all-zero rows and mismatched samples are dropped consistently (closes #82, closes #37)
- `extract_data()`: `age` values must not all be identical; `age_uncertainty` columns must not all be identical (closes #81)

## Documentation and site

- Roxygen2 updated from 7.2.3 → 7.3.3; all `man/*.Rd` files regenerated
- `pkgdown` site rebuilt; reference index and article pages updated
- `quarto` added to `Suggests` to support vignette rendering

## New contributor

- Friederike Wolke (@FriedaRosa) added as contributor (`ctb`) for authoring the test suite and filing the issues that drove this release's validation improvements

# RRatepol 1.2.3

- fix an issue with subsetting the uncertainty matrix to correctly align with the rest of the data (thanks to Giacomo Galli for reporting the bug)
- replace all not a basic ASCII character in documentation

# RRatepol 1.2.2

- add more information about additional resources and workshops
- use {neotoma2} from CRAN
- fix an issue with a case, where the age difference between samples is zero. The zeros are replaced by arbitrary 0.1.

# RRatepol 1.2.1

## Uncertainty sampling

- fix an issue with an error message for a large number of iterations in sampling age-sequence

# RRatepol 1.2.0

## Name uniformity

- In order to make all (exported) functions easier to use (use verbs, snake_style,...), rename all functions and their argument. For back-compatibility, the original functions and arguments are kept but flagged as deprecated. More details [here](https://github.com/HOPE-UIB-BIO/R-Ratepol-package/milestone/1?closed=1).
- main exported function name change:
  - rename `fc_estimate_RoC()` to `estimate_roc()` ([#23](https://github.com/HOPE-UIB-BIO/R-Ratepol-package/pull/23))
  - rename `fc_plot_RoC_sequence()` to `plot_roc()` ([#24](https://github.com/HOPE-UIB-BIO/R-Ratepol-package/pull/24))
  - rename `fc_detect_peak_points()` to `detect_peak_points()` ([#25](https://github.com/HOPE-UIB-BIO/R-Ratepol-package/pull/25))
  - rename all internal functions ([#26](https://github.com/HOPE-UIB-BIO/R-Ratepol-package/pull/26))
- Update README ([#27](https://github.com/HOPE-UIB-BIO/R-Ratepol-package/pull/27))

# RRatepol 1.1.0

## `neotoma` -> `neotoma2`

- update the vignette to use `neotoma2` package (`neotoma` was deprecated)

## utility functions  Internal checks

- change `util_*` functions to use [`RUtilpol` package](https://github.com/HOPE-UIB-BIO/R-Utilpol-package)

## Function references

- remove several functions references by flagging them as `internal`

## Bug fixes

- correctly not transform to proportions when selected by a user.

## Other

- within `select()` replace `.data$foo` with `"foo"` (suggested in `tidyselect`)
- replace `1:length(...` with `seq_along`

# RRatepol 1.0.0

## General

- package has been re-factored (most internal processes have been updated)
- several arguments of functions have been changed and marked with `lifecycle`
- add several more checks while data sourcing

## Dissimilarity calculation

- use `vegan` package to calculate all dissimilarities except of standardised euclidean.

## parallel computation

- use the `pbapply` package
- data for iteration is prepared separately
- added progress bar

## example data

- names for example columns as been adjusted to avoid using `.`

# RRatepol 0.6.1

## Dissimilarity calculation

- use `vegan` package to calculate Gower and Bray-Curtis (new) dissimilarity

## Univariate RoC estimation

- fix issues when using only single-taxon dataset

## Code cosmetics changes

- TRUE/FALSE written as full
- all arguments are written for each function
- spelling fixes

# RRatepol 0.6.0

## Overhaul change of the estimation of RoC

- RoC score is now returned in units selected by user (time_standardisation)
- Method of sample selection from bin is now optional (bin_selection)
- RoC is estimated only for subsequent bins (only_subsequent == TRUE)
- Bins are now created form the beginning of the core (instead of 0)

## Console outputs

- General overhaul of console outputs presented to user while running RRatepol
- more information and warning messages are presented with each run

## doParallel

have been curated by you during the last years packages `doSNOW` and `snow` are replaced by `doParallel`

- progress bar for randomisation is currently not present

## Other

- Added various checks for correct argument selection
- all hyphens are replaced with dashes (`–`) in the whole package
- examples are wrapped into `dontrun`

# RRatepol 0.5.6

## README

- Figures are now saved in the new folder (`man` folder) with the new names.
- cosmetics changes in the code (spaces, new lines, etc)

## Vignette

- stop using pre-saved data
- decrease the `i_multiplier` to 0.5 to speed up the vignette building
- change the `treads` arguments in `fc_estimate_RoC` to FALSE so the vignette can be build in machines without multiple cores.
- rename r chunks to remove empty spaces in names

## Other

- remove `RRatepol` namespace inside package for internal functions
- fix couple of typos in function descriptions
- add cran-comments

# RRatepol 0.5.5

- `sample.id` or `sample_id` can be used in as the sample identification
- RRatepol::fc_plot_RoC_sequence: fix typo in argument `age_threshold`
- citation(package = "RRatepol") now show correct citation
- update the Figure in README file
