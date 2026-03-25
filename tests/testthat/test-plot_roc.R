# ================================================================ #
# 1. data_source: input validation                                  #
# ================================================================ #

testthat::test_that(
  "plot_roc() errors on missing data_source",
  {
    testthat::expect_error(
      plot_roc(
        data_source = ,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = NULL
      ),
      'argument "data_source" is missing, with no default'
    )
  }
)

testthat::test_that(
  "plot_roc() errors on invalid data_source type",
  {
    purrr::walk(
      .x = list(NULL, "data_source", 123, 0, NA, list()),
      .f = function(bad_input) {
        testthat::expect_error(
          plot_roc(
            data_source = bad_input,
            age_threshold = NULL,
            roc_threshold = NULL,
            peaks = FALSE,
            trend = NULL
          ),
          "'data_source' must be one of the following: 'data.frame'"
        )
      }
    )
  }
)

testthat::test_that(
  "plot_roc() errors when data_source lacks required columns",
  {
    testthat::expect_error(
      plot_roc(
        data_source = data.frame(),
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = NULL
      ),
      paste0(
        "'data_source' must contains following columns:",
        " 'Age', 'ROC', 'ROC_up', 'ROC_dw'"
      )
    )
  }
)

testthat::test_that(
  "plot_roc() errors on zero-row data_source",
  {
    testthat::expect_error(
      suppressWarnings(
        plot_roc(
          data_source = data.frame(
            Age = numeric(0),
            ROC = numeric(0),
            ROC_up = numeric(0),
            ROC_dw = numeric(0)
          ),
          age_threshold = NULL,
          roc_threshold = NULL,
          peaks = FALSE,
          trend = NULL
        )
      ),
      # with age_threshold == NULL it uses max(Age) -> NULL
      "'to' must be a finite number"
    )
  }
)

# ================================================================ #
# 2. age_threshold: input validation                                #
# ================================================================ #

testthat::test_that(
  "plot_roc() accepts missing or NULL age_threshold",
  {
    data_source <-
      make_plot_roc_data()

    testthat::expect_no_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = NULL
      )
    )
  }
)

testthat::test_that(
  "plot_roc() errors on invalid age_threshold type",
  {
    data_source <-
      make_plot_roc_data()

    purrr::walk(
      .x = list("8000", NA, list(), data.frame()),
      .f = function(bad_val) {
        testthat::expect_error(
          plot_roc(
            data_source = data_source,
            age_threshold = bad_val,
            roc_threshold = NULL,
            peaks = FALSE,
            trend = NULL
          ),
          "'age_threshold' must be one of the following: 'NULL', 'numeric'"
        )
      }
    )
  }
)

testthat::test_that(
  "plot_roc() errors when age_threshold has length > 1",
  {
    data_source <-
      make_plot_roc_data()

    testthat::expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = c(100, 8000),
        roc_threshold = NULL,
        peaks = FALSE,
        trend = NULL
      ),
      "'to' must be of length 1"
    )
  }
)

testthat::test_that(
  "plot_roc() errors on negative age_threshold",
  {
    data_source <-
      make_plot_roc_data()

    testthat::expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = -8000,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = NULL
      ),
      "wrong sign in 'by' argument"
    )
  }
)

# ================================================================ #
# 3. roc_threshold: input validation                                #
# ================================================================ #

testthat::test_that(
  "plot_roc() accepts missing or NULL roc_threshold",
  {
    data_source <-
      make_plot_roc_data()

    testthat::expect_no_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = NULL
      )
    )
  }
)

testthat::test_that(
  "plot_roc() errors on invalid roc_threshold type",
  {
    data_source <-
      make_plot_roc_data()

    purrr::walk(
      .x = list("1", NA, list(), data.frame()),
      .f = function(bad_val) {
        testthat::expect_error(
          plot_roc(
            data_source = data_source,
            age_threshold = NULL,
            roc_threshold = bad_val,
            peaks = FALSE,
            trend = NULL
          ),
          "'roc_threshold' must be one of the following: 'NULL', 'numeric'"
        )
      }
    )
  }
)

testthat::test_that(
  "plot_roc() errors when roc_threshold has length > 1",
  {
    data_source <-
      make_plot_roc_data()

    testthat::expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = c(1, 2),
        peaks = FALSE,
        trend = NULL
      ),
      "`ylim` must be a vector of length 2"
    )
  }
)

# ================================================================ #
# 4. peaks: input validation                                        #
# ================================================================ #

testthat::test_that(
  "plot_roc() accepts missing peaks and uses FALSE as default",
  {
    data_source <-
      make_plot_roc_data()

    testthat::expect_no_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = , # uses FALSE as default
        trend = NULL
      )
    )
  }
)

testthat::test_that(
  "plot_roc() errors on invalid non-logical peaks type",
  {
    data_source <-
      make_plot_roc_data()

    purrr::walk(
      .x = list(NULL, "FALSE", 1, list(), data.frame()),
      .f = function(bad_val) {
        testthat::expect_error(
          plot_roc(
            data_source = data_source,
            age_threshold = NULL,
            roc_threshold = NULL,
            peaks = bad_val,
            trend = NULL
          ),
          "'peaks' must be one of the following: 'logical'"
        )
      }
    )
  }
)

# ================================================================ #
# 5. trend: input validation                                        #
# ================================================================ #

testthat::test_that(
  "plot_roc() accepts missing or NULL trend",
  {
    data_source <-
      make_plot_roc_data()

    testthat::expect_no_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = NULL
      )
    )
  }
)

testthat::test_that(
  "plot_roc() errors on invalid non-character trend type",
  {
    data_source <-
      make_plot_roc_data()

    purrr::walk(
      .x = list(FALSE, 1, NA, list(), data.frame()),
      .f = function(bad_val) {
        testthat::expect_error(
          plot_roc(
            data_source = data_source,
            age_threshold = NULL,
            roc_threshold = NULL,
            peaks = FALSE,
            trend = bad_val
          ),
          "'trend' must be one of the following: 'NULL', 'character'"
        )
      }
    )
  }
)

testthat::test_that(
  "plot_roc() errors on invalid trend method name",
  {
    data_source <-
      make_plot_roc_data()

    testthat::expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = "my_trend"
      ),
      paste0(
        "'trend' must contains one of the following values:",
        " 'threshold', 'trend_linear', 'trend_non_linear'"
      )
    )
  }
)

testthat::test_that(
  "plot_roc() errors when trend has length > 1",
  {
    data_source <-
      make_plot_roc_data()

    testthat::expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = c("threshold", "trend_linear")
      ),
      "the condition has length > 1"
    )
  }
)

# ================================================================ #
# 6. Functionality                                                   #
# ================================================================ #

testthat::test_that(
  "plot_roc() returns a ggplot object",
  {
    data_source <-
      make_plot_roc_data()

    res <-
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = NULL
      )

    testthat::expect_s3_class(res, "ggplot")
  }
)

testthat::test_that(
  "plot_roc() correctly displays peak points",
  {
    data_source <-
      make_plot_roc_data()

    p <-
      plot_roc(
        data_source = data_source,
        age_threshold = 8000,
        roc_threshold = NULL,
        peaks = TRUE,
        trend = NULL
      )

    built <-
      ggplot2::ggplot_build(p)

    point_data <-
      built$data[[which(
        sapply(
          built$plot$layers, function(l) {
            inherits(l$geom, "GeomPoint")
          }
        )
      )]]

    expected_peaks <-
      dplyr::filter(data_source, Peak == TRUE)

    testthat::expect_equal(
      base::nrow(point_data),
      base::nrow(expected_peaks)
    )

    testthat::expect_true(
      base::all(point_data$colour == "green")
    )
  }
)

testthat::test_that(
  "plot_roc() sets proper axis limits from threshold parameters",
  {
    data_source <-
      estimate_roc(
        data_source_community =
          RRatepol::example_data$pollen_data[[1]],
        data_source_age =
          RRatepol::example_data$sample_age[[1]],
        smooth_method = "shep",
        smooth_n_points = 5,
        working_units = "levels",
        standardise = FALSE,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        use_parallel = FALSE,
        silent = TRUE
      )

    age_threshold <- 8000
    roc_threshold <- 2

    p <-
      plot_roc(
        data_source = data_source,
        age_threshold = age_threshold,
        roc_threshold = roc_threshold,
        peaks = FALSE,
        trend = NULL
      )

    testthat::expect_s3_class(p, "ggplot")

    # Check the limits stored in the coord_flip object directly —
    # panel_params includes ggplot2 axis expansion so those values differ.
    # In coord_flip, $limits$x corresponds to Age (xlim) and
    # $limits$y corresponds to ROC (ylim).
    coord_limits <-
      p$coordinates$limits

    testthat::expect_equal(
      base::sort(coord_limits$x),
      c(0, age_threshold)
    )
    testthat::expect_equal(
      coord_limits$y,
      c(0, roc_threshold)
    )
  }
)

testthat::test_that(
  "plot_roc() adds trend line when trend is specified",
  {
    data_source <-
      estimate_roc(
        data_source_community =
          RRatepol::example_data$pollen_data[[1]],
        data_source_age =
          RRatepol::example_data$sample_age[[1]],
        smooth_method = "shep",
        smooth_n_points = 5,
        working_units = "levels",
        standardise = FALSE,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        use_parallel = FALSE,
        silent = TRUE
      ) |>
      detect_peak_points(
        sel_method = "trend_linear",
        sd_threshold = 2
      )

    p <-
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = TRUE,
        trend = "trend_linear"
      )

    # base layers + point layer + trend layer
    testthat::expect_true(base::length(p$layers) > 4)
  }
)

testthat::test_that(
  "plot_roc() produces different plots for different trend methods",
  {
    data_source <-
      make_plot_roc_data()

    p_linear <-
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = TRUE,
        trend = "trend_linear"
      )
    p_nonlinear <-
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = TRUE,
        trend = "trend_non_linear"
      )
    p_threshold <-
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = TRUE,
        trend = "threshold"
      )

    testthat::expect_false(identical(p_linear, p_nonlinear))
    testthat::expect_false(identical(p_linear, p_threshold))
    testthat::expect_false(identical(p_threshold, p_nonlinear))
  }
)
