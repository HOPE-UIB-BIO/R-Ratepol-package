# Parameters:
## data_source (estimate_roc output)
## age_threshold (numeric; max age)
## roc_threshold (numeric;max roc)
## peaks (logical)
## trend (threshold, trend_linear, trend_non_linear)

data_source <-
  estimate_roc(
    data_source_community = RRatepol::example_data$pollen_data[[1]],
    data_source_age = RRatepol::example_data$sample_age[[1]],
    age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
    smooth_method = "grim",
    smooth_n_points = 5,
    smooth_age_range = 500,
    smooth_n_max = 9,
    working_units = "levels",
    bin_size = 500,
    number_of_shifts = 1,
    bin_selection = "first",
    standardise = TRUE,
    n_individuals = 150,
    dissimilarity_coefficient = "euc",
    tranform_to_proportions = TRUE,
    rand = 10,
    use_parallel = FALSE,
    interest_threshold = NULL,
    time_standardisation = NULL,
    verbose = FALSE
  ) %>%
  detect_peak_points(
    sel_method = "trend_linear"
  )

plot_roc(
  data_source = data_source,
  age_threshold = NULL,
  roc_threshold = NULL,
  peaks = FALSE,
  trend = NULL
)


test_that(
  "plot_roc throws error with empty data_source",
  {
    expect_error(
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

# NULL
test_that(
  "plot_roc throws error with NULL data_source",
  {
    expect_error(
      plot_roc(
        data_source = NULL,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = NULL
      ),
      "'data_source' must be one of the following: 'data.frame'"
    )
  }
)

# character
test_that(
  "plot_roc() rejects character as data_source",
  {
    expect_error(
      plot_roc(
        data_source = "data_source",
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = NULL
      ),
      "'data_source' must be one of the following: 'data.frame'"
    )
  }
)

# Numeric
test_that(
  "plot_roc() rejects numeric as data_source",
  {
    expect_error(
      plot_roc(
        data_source = 123,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = NULL
      ),
      "'data_source' must be one of the following: 'data.frame'"
    )
  }
)

# zero
test_that(
  "plot_roc() rejects zero value as data_source",
  {
    expect_error(
      plot_roc(
        data_source = 0,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = NULL
      ),
      "'data_source' must be one of the following: 'data.frame'"
    )
  }
)

# NA
test_that(
  "plot_roc() rejects NA as data_source",
  {
    expect_error(
      plot_roc(
        data_source = NA,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = NULL
      ),
      "'data_source' must be one of the following: 'data.frame'"
    )
  }
)

# empty list
test_that(
  "plot_roc() rejects empty list as data_source",
  {
    expect_error(
      plot_roc(
        data_source = list(),
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = NULL
      ),
      "'data_source' must be one of the following: 'data.frame'"
    )
  }
)
# empty data.frame
test_that(
  "plot_roc() rejects empty data.frame as data_source",
  {
    expect_error(
      plot_roc(
        data_source = data.frame(),
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = NULL
      ),
      "'data_source' must contains following columns: 'Age', 'ROC', 'ROC_up', 'ROC_dw'"
    )
  }
)

# 0 row dataframe
test_that(
  "plot_roc() rejects zero-row data.frame as data_source",
  {
    expect_error(
      suppressWarnings(
        plot_roc(
          data_source = data.frame(
            Age = numeric(
              0
            ), ROC = numeric(
              0
            ), ROC_up = numeric(
              0
            ), ROC_dw = numeric(
              0
            )
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

## age_threshold
# empty
test_that(
  "plot_roc() accepts missing age_threshold parameter and uses NULL as default",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_no_error(
      plot_roc(
        data_source = data_source,
        age_threshold = , # uses NULL as default
        roc_threshold = NULL,
        peaks = FALSE,
        trend = NULL
      )
    )
  }
)

# NULL
test_that(
  "plot_roc() uses max age when age_threshold is NULL",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_no_error(
      res <-
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

# character
test_that(
  "plot_roc() rejects character as age_threshold",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = "8000",
        roc_threshold = NULL,
        peaks = FALSE,
        trend = NULL
      ),
      "'age_threshold' must be one of the following: 'NULL', 'numeric'"
    )
  }
)

# NA
test_that(
  "plot_roc() rejects NA as age_threshold",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NA,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = NULL
      ),
      "'age_threshold' must be one of the following: 'NULL', 'numeric'"
    )
  }
)
# multiple
test_that(
  "plot_roc() rejects vector as age_threshold",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = c(
          100, 8000
        ),
        roc_threshold = NULL,
        peaks = FALSE,
        trend = NULL
      ),
      "'to' must be of length 1"
    )
  }
)

# negative
test_that(
  "plot_roc() rejects negative age_threshold",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
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

# empty list
test_that(
  "plot_roc() rejects list as age_threshold",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = list(),
        roc_threshold = NULL,
        peaks = FALSE,
        trend = NULL
      ),
      "'age_threshold' must be one of the following: 'NULL', 'numeric'"
    )
  }
)

# empty data.frame
test_that(
  "plot_roc() rejects data.frame as age_threshold",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = data.frame(),
        roc_threshold = NULL,
        peaks = FALSE,
        trend = NULL
      ),
      "'age_threshold' must be one of the following: 'NULL', 'numeric'"
    )
  }
)

# roc_threshold
# empty
test_that(
  "plot_roc() accepts missing roc_threshold parameter",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_no_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = , # uses NULL as default
        peaks = FALSE,
        trend = NULL
      )
    )
  }
)

# NULL (no error)
test_that(
  "plot_roc() uses max ROC when roc_threshold is NULL",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_no_error(
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

# character
test_that(
  "plot_roc() rejects character as roc_threshold",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = "1",
        peaks = FALSE,
        trend = NULL
      ),
      "'roc_threshold' must be one of the following: 'NULL', 'numeric'"
    )
  }
)

# NA
test_that(
  "plot_roc() rejects NA as roc_threshold",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NA,
        peaks = FALSE,
        trend = NULL
      ),
      "'roc_threshold' must be one of the following: 'NULL', 'numeric'"
    )
  }
)

# multiple
test_that(
  "plot_roc() rejects vector as roc_threshold",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = c(
          1, 2
        ),
        peaks = FALSE,
        trend = NULL
      ),
      "`ylim` must be a vector of length 2, not a double vector of length 3."
    )
  }
)

# empty list
test_that(
  "plot_roc() rejects list as roc_threshold",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = list(),
        peaks = FALSE,
        trend = NULL
      ),
      "'roc_threshold' must be one of the following: 'NULL', 'numeric'"
    )
  }
)

# empty data.frame
test_that(
  "plot_roc() rejects data.frame as roc_threshold",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = data.frame(),
        peaks = FALSE,
        trend = NULL
      ),
      "'roc_threshold' must be one of the following: 'NULL', 'numeric'"
    )
  }
)


# peaks tests
# empty
test_that(
  "plot_roc() accepts missing peaks parameter and uses FALSE as default",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_no_error(
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
# NULL
test_that(
  "plot_roc() rejects NULL as peaks",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = NULL,
        trend = NULL
      ),
      "'peaks' must be one of the following: 'logical'"
    )
  }
)


# character
test_that(
  "plot_roc() rejects character as peaks",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = "FALSE",
        trend = NULL
      ),
      "'peaks' must be one of the following: 'logical'"
    )
  }
)

# numeric
test_that(
  "plot_roc() rejects numeric as peaks",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = 1,
        trend = NULL
      ),
      "'peaks' must be one of the following: 'logical'"
    )
  }
)

# NA
test_that(
  "plot_roc() rejects NA as peaks",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = NA,
        trend = NULL
      ),
      # none programmed into the function - it just checks for logical and isFalse()
      "'peaks' must be one of the following: 'TRUE', 'FALSE'"
    )
  }
)

# multiple
test_that(
  "plot_roc() rejects vector as peaks",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = c(
          TRUE, FALSE
        ),
        trend = NULL
      ),
      # none programmed into the function
      # e.g.,
      "peaks argument must be of length 1"
    )
  }
)

# empty list
test_that(
  "plot_roc() rejects list as peaks",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = list(),
        trend = NULL
      ),
      "'peaks' must be one of the following: 'logical'"
    )
  }
)

# empty data.frame
test_that(
  "plot_roc() rejects data.frame as peaks",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = data.frame(),
        trend = NULL
      ),
      "'peaks' must be one of the following: 'logical'"
    )
  }
)


# Trend tests
# empty
test_that(
  "plot_roc() accepts missing trend parameter and uses NULL as default",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_no_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = # uses NULL as default
        )
    )
  }
)

# NULL (no error)
test_that(
  "plot_roc() accepts NULL as trend",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_no_error(
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

# FALSE
test_that(
  "plot_roc() rejects logical as trend",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = FALSE
      ),
      "'trend' must be one of the following: 'NULL', 'character'"
    )
  }
)

# invalid character
test_that(
  "plot_roc() rejects invalid trend method name",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = "my_trend"
      ),
      "'trend' must contains one of the following values: 'threshold', 'trend_linear', 'trend_non_linear'"
    )
  }
)

# multiple
test_that(
  "plot_roc() rejects multiple trend methods",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
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

# numeric
test_that(
  "plot_roc() rejects numeric as trend",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = 1
      ),
      "'trend' must be one of the following: 'NULL', 'character'"
    )
  }
)
# NA
test_that(
  "plot_roc() rejects NA as trend",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = NA
      ),
      "'trend' must be one of the following: 'NULL', 'character'"
    )
  }
)

# empty list
test_that(
  "plot_roc() rejects list as trend",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = list()
      ),
      "'trend' must be one of the following: 'NULL', 'character'"
    )
  }
)
# empty data.frame
test_that(
  "plot_roc() rejects data.frame as trend",
  {
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
    expect_error(
      plot_roc(
        data_source = data_source,
        age_threshold = NULL,
        roc_threshold = NULL,
        peaks = FALSE,
        trend = data.frame()
      ),
      "'trend' must be one of the following: 'NULL', 'character'"
    )
  }
)


# Functionality tests
test_that(
  "plot_roc() correctly displays peak points",
  {
    # Setup test data with known peaks
    set.seed(123)
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )

    # Generate plot with peaks
    p <-
      plot_roc(
        data_source,
        age_threshold = 8000, peaks = TRUE
      )

    # Build the plot
    built <-
      ggplot2::ggplot_build(p)

    # Check the point data - green points should correspond to peaks
    point_data <-
      built$data[[which(
        sapply(
          built$plot$layers, function(l) {
            inherits(
              l$geom, "GeomPoint"
            )
          }
        )
      )]]

    # Filter the original data for peaks
    expected_peaks <-
      dplyr::filter(data_source, Peak == TRUE)

    # Test that the number of points matches the number of peaks
    expect_equal(
      nrow(
        point_data
      ),
      nrow(expected_peaks)
    )

    # Test that colors are correct
    expect_true(
      all(
        point_data$colour == "green"
      )
    )
  }
)

test_that(
  "plot_roc() sets proper axis limits from threshold parameters",
  {
    # Create test data
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        smooth_method = "shep",
        smooth_n_points = 5,
        working_units = "levels",
        standardise = FALSE,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        use_parallel = FALSE,
        verbose = FALSE
      )

    # Define test parameters
    age_threshold <-
      8000
    roc_threshold <-
      2

    # Create plot
    p <-
      plot_roc(
        data_source,
        age_threshold = age_threshold,
        roc_threshold = roc_threshold
      )

    # Test plot structure
    expect_s3_class(
      p, "ggplot"
    )

    # Extract and test coordinate system
    built <-
      ggplot2::ggplot_build(p)
    panel <-
      built$layout$panel_params[[1]]

    # Test x-axis (age) limits - note that in coord_flip, x and y are switched
    expect_equal(
      panel$y.range[1], 0
    )
    expect_equal(
      panel$y.range[2], age_threshold
    )

    # Test y-axis (RoC) limits
    expect_equal(
      panel$x.range[1], 0
    )
    expect_equal(
      panel$x.range[2], roc_threshold
    )
  }
)

test_that(
  "plot_roc() adds trend line when trend parameter is specified",
  {
    # Create test data
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        smooth_method = "shep",
        smooth_n_points = 5,
        working_units = "levels",
        standardise = FALSE,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        use_parallel = FALSE,
        verbose = FALSE
      )

    # Add peaks
    data_with_peaks <-
      detect_peak_points(
        data_source = data_source,
        sel_method = "trend_linear",
        sd_threshold = 2
      )

    # Create plot with trend
    p <-
      plot_roc(
        data_with_peaks,
        peaks = TRUE,
        trend = "trend_linear"
      )

    # Count number of layers - should have extra line for trend
    expect_true(
      length(
        p$layers
      ) > 4
    ) # base layers + point layer + trend layer
  }
)

# different trend methods produce different results
test_that(
  "plot_roc() produces different plots for different trend methods",
  {
    set.seed(123)
    data_source <-
      estimate_roc(
        data_source_community = RRatepol::example_data$pollen_data[[1]],
        data_source_age = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = NULL,
        verbose = FALSE
      ) %>%
      detect_peak_points(
        sel_method = "trend_linear"
      )
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

    expect_false(
      identical(
        p_linear,
        p_nonlinear
      )
    )

    expect_false(
      identical(
        p_linear,
        p_threshold
      )
    )

    expect_false(
      identical(
        p_threshold,
        p_nonlinear
      )
    )
  }
)
