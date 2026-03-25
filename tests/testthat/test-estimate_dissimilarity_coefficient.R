# ================================================================ #
# 1. data_source_dc: input validation                              #
# ================================================================ #

testthat::test_that(
  "estimate_dissimilarity_coefficient() errors on missing data_source_dc",
  {
    testthat::expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = ,
        dissimilarity_coefficient = "euc"
      ),
      'argument "data_source_dc" is missing, with no default'
    )
  }
)

testthat::test_that(
  "estimate_dissimilarity_coefficient() errors on invalid data_source_dc",
  {
    purrr::walk(
      .x = list(NULL, "my_data", 123, list(), data.frame(), 0, NA),
      .f = function(bad_input) {
        testthat::expect_error(
          estimate_dissimilarity_coefficient(
            data_source_dc = bad_input,
            dissimilarity_coefficient = "euc"
          ),
          "invalid 'length' argument"
        )
      }
    )
  }
)

# ================================================================ #
# 2. dissimilarity_coefficient: input validation                   #
# ================================================================ #

testthat::test_that(
  "estimate_dissimilarity_coefficient() errors on NULL coefficient",
  {
    data_dc <-
      make_dc_data()

    testthat::expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = data_dc,
        dissimilarity_coefficient = NULL
      ),
      "argument is of length zero"
    )
  }
)

testthat::test_that(
  "estimate_dissimilarity_coefficient() errors on unknown coefficient",
  {
    data_dc <-
      make_dc_data()

    purrr::walk(
      .x = list("my_coefficient", 123, 0),
      .f = function(bad_coeff) {
        testthat::expect_error(
          estimate_dissimilarity_coefficient(
            data_source_dc = data_dc,
            dissimilarity_coefficient = bad_coeff
          ),
          "object 'corrmat' not found"
        )
      }
    )
  }
)

testthat::test_that(
  "estimate_dissimilarity_coefficient() errors on multiple coefficients",
  {
    data_dc <-
      make_dc_data()

    testthat::expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = data_dc,
        dissimilarity_coefficient = c("euc", "chord")
      ),
      "the condition has length > 1"
    )
  }
)

# ================================================================ #
# 3. Valid inputs: all coefficients x both data types              #
# ================================================================ #

testthat::test_that(
  "estimate_dissimilarity_coefficient() returns numeric for all coefficients",
  {
    vec_coefficients <-
      c("euc", "euc.sd", "chord", "chisq", "gower", "bray")

    purrr::walk(
      .x = c("proportions", "percentages"),
      .f = function(sel_method) {
        data_dc <-
          make_dc_data(sel_method = sel_method)

        purrr::walk(
          .x = vec_coefficients,
          .f = function(coeff) {
            result <-
              estimate_dissimilarity_coefficient(
                data_source_dc = data_dc,
                dissimilarity_coefficient = coeff
              )
            testthat::expect_type(result, "double")
          }
        )
      }
    )
  }
)
