# ================================================================ #
# 1. Invalid CharData types                                        #
# ================================================================ #

testthat::test_that(
  "detect_sni() errors on invalid CharData types",
  {
    valid_bw <- 500

    purrr::walk(
      .x = list(NULL, "my_data", 123, 0, NA, list()),
      .f = function(bad_input) {
        testthat::expect_error(
          detect_sni(
            CharData = bad_input,
            BandWidth = valid_bw
          )
        )
      }
    )
  }
)

# ================================================================ #
# 2. Invalid BandWidth values                                      #
# ================================================================ #

testthat::test_that(
  "detect_sni() errors on NULL BandWidth",
  {
    data_source <-
      make_roc_data()

    pred_gam <-
      make_trend(
        data_source = data_source,
        sel_method = "non_linear"
      )

    char_data <-
      data.frame(
        dplyr::pull(data_source, Age),
        dplyr::pull(data_source, ROC),
        pred_gam
      )

    testthat::expect_error(
      detect_sni(
        CharData = char_data,
        BandWidth = NULL
      ),
      "argument is of length zero"
    )
  }
)

testthat::test_that(
  "detect_sni() errors on non-numeric or non-positive BandWidth",
  {
    data_source <-
      make_roc_data()

    pred_gam <-
      make_trend(
        data_source = data_source,
        sel_method = "non_linear"
      )

    char_data <-
      data.frame(
        dplyr::pull(data_source, Age),
        dplyr::pull(data_source, ROC),
        pred_gam
      )

    testthat::expect_error(
      detect_sni(CharData = char_data, BandWidth = "5"),
      "non-numeric argument to binary operator"
    )

    testthat::expect_error(
      detect_sni(CharData = char_data, BandWidth = 0),
      "'f' must be finite and > 0"
    )

    testthat::expect_error(
      detect_sni(CharData = char_data, BandWidth = -5),
      "'f' must be finite and > 0"
    )
  }
)

# ================================================================ #
# 3. Valid input: output structure                                  #
# ================================================================ #

testthat::test_that(
  "detect_sni() returns list with expected names and types",
  {
    data_source <-
      make_roc_data()

    pred_gam <-
      make_trend(
        data_source = data_source,
        sel_method = "non_linear"
      )

    char_data <-
      data.frame(
        dplyr::pull(data_source, Age),
        dplyr::pull(data_source, ROC),
        pred_gam
      )

    band_width <-
      5 * base::mean(
        base::diff(
          dplyr::pull(data_source, Age)
        )
      )

    result <-
      detect_sni(
        CharData = char_data,
        BandWidth = band_width
      )

    testthat::expect_type(result, "list")

    testthat::expect_true(
      base::all(
        c("SNI", "winInd", "popN", "popS", "meanN", "stdN", "CF") %in%
          base::names(result)
      )
    )

    testthat::expect_type(
      purrr::chuck(result, "SNI"),
      "double"
    )

    testthat::expect_equal(
      base::length(purrr::chuck(result, "SNI")),
      base::nrow(char_data)
    )
  }
)
