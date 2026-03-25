# ================================================================ #
# 1. Required arguments: input validation                          #
# ================================================================ #

testthat::test_that(
  "extract_data() errors on missing data_community_extract",
  {
    testthat::expect_error(
      extract_data(
        data_community_extract = ,
        data_age_extract =
          RRatepol::example_data$sample_age[[1]],
        silent = TRUE
      ),
      'argument "data_community_extract" is missing, with no default'
    )
  }
)

testthat::test_that(
  "extract_data() errors on invalid data types",
  {
    purrr::walk(
      .x = list(NULL, "my_data", 123, NA),
      .f = function(bad_input) {
        testthat::expect_error(
          extract_data(
            data_community_extract = bad_input,
            data_age_extract =
              RRatepol::example_data$sample_age[[1]],
            silent = TRUE
          )
        )
      }
    )
  }
)

# ================================================================ #
# 2. Valid output: structure and column names                      #
# ================================================================ #

testthat::test_that(
  "extract_data() returns list with expected structure",
  {
    result <-
      suppressWarnings(
        extract_data(
          data_community_extract =
            RRatepol::example_data$pollen_data[[1]],
          data_age_extract =
            RRatepol::example_data$sample_age[[1]],
          silent = TRUE
        )
      )

    testthat::expect_type(result, "list")
    testthat::expect_true(
      base::all(
        c("community", "age", "age_un") %in% base::names(result)
      )
    )
    testthat::expect_s3_class(
      purrr::chuck(result, "community"), "data.frame"
    )
    testthat::expect_true(
      base::nrow(purrr::chuck(result, "community")) > 0
    )
  }
)

testthat::test_that(
  "extract_data() errors when age is all identical",
  {
    data_community <-
      RRatepol::example_data$pollen_data[[1]]

    data_age <-
      RRatepol::example_data$sample_age[[1]]

    data_age$age <- base::rep(1000, base::nrow(data_age))

    testthat::expect_error(
      extract_data(
        data_community_extract = data_community,
        data_age_extract = data_age,
        silent = TRUE
      ),
      "'age' values must not all be identical across samples"
    )
  }
)

testthat::test_that(
  "extract_data() errors when age_uncertainty columns are all identical",
  {
    data_community <-
      RRatepol::example_data$pollen_data[[1]]

    data_age <-
      RRatepol::example_data$sample_age[[1]]

    n_samples <-
      base::nrow(data_age)

    # each column has a single repeated value (all-identical)
    age_un_identical <-
      base::matrix(
        base::rep(1000, 10 * n_samples),
        nrow = 10,
        ncol = n_samples
      )

    testthat::expect_error(
      extract_data(
        data_community_extract = data_community,
        data_age_extract = data_age,
        age_uncertainty = age_un_identical,
        silent = TRUE
      ),
      "'age_uncertainty' columns must not all be identical"
    )
  }
)
