# ---------------------------------------------------------- #
#               Data_source_prep Input Tests                 #
#                            Errors                          #
# ---------------------------------------------------------- #
# 1.1 data_source_prep = NULL
# i) levels
test_that(
  "prepare_data, working_units = 'levels' throws error for NULL input for data_source",
  {
    expect_error(
      prepare_data(
        data_source_prep = NULL,
        working_units = "levels",
        rand = NULL
      ), "'data_source_prep' must be one of the following: 'list'"
    )
  }
)

# ii) bins
test_that(
  "prepare_data, working_units = 'bins' throws error for NULL input for data_source",
  {
    expect_error(
      prepare_data(
        data_source_prep = NULL,
        working_units = "bins",
        bin_size = 500,
        rand = NULL
      ), "'data_source_prep' must be one of the following: 'list'"
    )
  }
)

# iii) MW
test_that(
  "prepare_data, working_units = 'MW' throws error for NULL input for data_source",
  {
    expect_error(
      prepare_data(
        data_source_prep = NULL,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      ), "'data_source_prep' must be one of the following: 'list'"
    )
  }
)

# --------------------------------------------------- #

# 1.2 data_source_prep = empty list with data = zero
# (these fail because there is no error programmed into the function yet)
# i) levels
test_that(
  "prepare_data, working_units = 'levels' throws error for data_source_prep = empty named list",
  {
    expect_error(
      prepare_data(
        data_source_prep =
          list(
            community = data.frame(
              0
            ),
            age = data.frame(
              age = 0
            ),
            age_un = data.frame(0)
          ),
        working_units = "levels",
        rand = NULL
      ),
      # none programmed into the function yet.
      # e.g., "'data_source_prep' must be non-empty"
    )
  }
)

# ii) bins
test_that(
  "prepare_data, working_units = 'bins' throws error for data_source_prep = empty named list",
  {
    expect_error(
      prepare_data(
        data_source_prep =
          list(
            community = data.frame(
              0
            ),
            age = data.frame(
              age = 0
            ),
            age_un = data.frame(0)
          ),
        working_units = "bins",
        bin_size = 500,
        rand = NULL
      ),
      # none programmed into the function yet.
      # e.g., "'data_source_prep' must be non-empty"
    )
  }
)

# iii) MW
test_that(
  "prepare_data, working_units = 'MW' throws error for data_source_prep = empty named list",
  {
    expect_error(
      prepare_data(
        data_source_prep =
          list(
            community = data.frame(
              0
            ),
            age = data.frame(
              age = 0
            ),
            age_un = data.frame(0)
          ),
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      ),
      # none programmed into the function yet.
      # e.g., "'data_source_prep' must be non-empty"
    )
  }
)

# ----------------------------------------------------------------- #
# ----------------------------------------------------------------- #


test_that(
  "prepare_data rejects NULL data_source",
  {
    expect_error(
      prepare_data(
        data_source_prep = NULL,
        working_units = c("levels", "bins", "MW"),
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      ), "'data_source_prep' must be one of the following: 'list'"
    )
  }
)

test_that(
  "prepare_data rejects incomplete list(
    data.frame) in data_source",
  {
    expect_error(
      prepare_data(
        data_source_prep = RRatepol::example_data$pollen_data[[1]],
        working_units = c("levels", "bins", "MW"),
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      ), "'data_source_prep' must be one of the following: 'list'"
    )
  }
)

test_that(
  "prepare_data throws error with empty data_source_prep list",
  {
    expect_error(
      prepare_data(
        data_source_prep = list(),
        working_units = "levels"
      ),
      "argument of length 0"
    )
  }
)

test_that(
  "prepare_data  throws error with data_source_prep with wrong structure",
  {
    wrong_structure <-
      list(
        wrong_name = data.frame(
          x = 1:5
        ),
        another_wrong = data.frame(y = 1:5)
      )

    expect_error(
      prepare_data(
        data_source_prep = wrong_structure,
        working_units = "levels"
      ),
      "argument of length 0"
    )
  }
)

# ----------------------------------------------------------------- #
#                 input data internal structure validation
# ----------------------------------------------------------------- #

test_that(
  "prepare_data fails if no community in data_source",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )
    example_data$community <-
      NULL
    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = c("levels", "bins", "MW"),
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      ), "is not TRUE"
    )
  }
)

test_that(
  "prepare_data fails if no age in data_source",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )
    example_data$age <-
      NULL
    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = c("levels", "bins", "MW"),
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      ), "`names` must be a character vector"
    )
  }
)

# no error without age_uncertainty
test_that(
  "prepare_data works if no age_uncertainty in data_source",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )
    example_data$age_un <-
      NULL
    expect_no_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = c("levels", "bins", "MW"),
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      )
    )
  }
)

# ----------------------------------------------------------------- #
#   More specific tests for 1) community 2) age and 3) age_un       #
# ----------------------------------------------------------------- #

## Community data
test_that(
  "prepare_data handles single-row community data",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]][1, , drop = FALSE],
        data_age_extract = RRatepol::example_data$sample_age[[1]][1, , drop = FALSE],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]][, 1, drop = FALSE],
        verbose = FALSE
      )

    expect_no_error(
      result <-
        prepare_data(
          data_source_prep = example_data,
          working_units = "levels"
        )
    )

    expect_equal(nrow(result[[1]][[1]]$bins), 1)
  }
)

test_that(
  "prepare_data handles single-column community data",
  {
    community_single <-
      RRatepol::example_data$pollen_data[[1]][, 1:2] # sample_id + 1 species

    example_data <-
      extract_data(
        data_community_extract = community_single,
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = FALSE
      )

    expect_no_error(
      result <-
        prepare_data(
          data_source_prep = example_data,
          working_units = "levels"
        )
    )

    expect_equal(
      ncol(
        result[[1]][[1]]$data
      ) - 1, 1
    ) # -1 for age column
  }
)

# Zeros
test_that(
  "prepare_data with working_units='levels' works if there are 0s in community data",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    example_data$community[1:3] <-
      0

    expect_no_error(
      prepare_data(
        example_data,
        working_units = "levels"
      )
    )
  }
)

test_that(
  "prepare_data with working_units='bins' works if there are 0s in community data",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    example_data$community[1:3] <-
      0

    expect_no_error(
      prepare_data(
        example_data,
        working_units = "bins",
        bin_size = 500
      )
    )
  }
)

test_that(
  "prepare_data with working_units='MW' works if there are 0s in community data",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    example_data$community[1:3] <-
      0

    expect_no_error(
      prepare_data(
        example_data,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5
      )
    )
  }
)

# NAs in Community
test_that(
  "prepare_data with working_units='levels' works if there are NAs in community data",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    example_data$community[1:3] <-
      NA

    expect_no_error(
      prepare_data(
        example_data,
        working_units = "levels"
      )
    )
  }
)

test_that(
  "prepare_data with working_units='bins' works if there are NAs in community data",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    example_data$community[1:3] <-
      NA

    expect_no_error(
      prepare_data(
        example_data,
        working_units = "bins",
        bin_size = 500
      )
    )
  }
)

test_that(
  "prepare_data with working_units='MW' works if there are NAs in community data",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    example_data$community[1:3] <-
      NA

    expect_no_error(
      prepare_data(
        example_data,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5
      )
    )
  }
)

## Age data
# Zeros
test_that(
  "prepare_data with working_units='levels' works if there are 0s in age data",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    example_data$age$age[1:3] <-
      0

    expect_no_error(
      prepare_data(
        example_data,
        working_units = "levels"
      )
    )
  }
)

test_that(
  "prepare_data with working_units='bins' works if there are 0s in age data",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    example_data$age$age[1:3] <-
      0

    expect_no_error(
      res <-
        prepare_data(
          example_data,
          working_units = "bins",
          bin_size = 500
        )
    )
  }
)

test_that(
  "prepare_data with working_units='MW' works if there are 0s in age data",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    example_data$age$age[1:3] <-
      0

    expect_no_error(
      prepare_data(
        example_data,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5
      )
    )
  }
)

# ---------------------------------------------------------- #
#                  Working Units Input Tests                 #
# ---------------------------------------------------------- #

test_that(
  "prepare_data validates working_units is one of the three valid options",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "invalid",
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      ),
      "'working_units' must contains one of the following values: 'levels', 'bins', 'MW'"
    )
  }
)

test_that(
  "prepare_data validates working_units parameter is not numeric",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = 123,
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      ),
      "'working_units' must be one of the following: 'character'"
    )
  }
)

test_that(
  "prepare_data validates working_units parameter is not numeric",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = NA,
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      ),
      "'working_units' must be one of the following: 'character'"
    )
  }
)

test_that(
  "prepare_data validates working_units parameter is not a vector of parameters",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = c("MW", "bins", "levels"),
        rand = NULL
      ),
      "'arg' must be of length 1"
    )
  }
)

test_that(
  "prepare_data validates working_units parameter is not a vector of parameters",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = c("MW", "bins")
      ),
      "'arg' must be of length 1"
    )
  }
)

# --------------------------------------------------- #

# 1.3 working_units = NULL
# 1.3.1 m.avg
test_that(
  "prepare_data, smooth_method = 'm.avg' throws error for working_units = NULL",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_work_mavg <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "m.avg",
        smooth_n_points = 5
      ) %>%
      reduce_data()

    expect_error(
      prepare_data(
        data_source_prep = data_work_mavg,
        working_units = NULL,
        bin_size = NULL,
        number_of_shifts = NULL,
        rand = NULL
      ), "'working_units' must be one of the following: 'character'"
    )
  }
)

# 1.3.2. shep
test_that(
  "prepare_data, smooth_method = 'shep' throws error for working_units = NULL",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_work_shep <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "shep"
      ) %>%
      reduce_data()

    expect_error(
      prepare_data(
        data_source_prep = data_work_shep,
        working_units = NULL,
        bin_size = NULL,
        number_of_shifts = NULL,
        rand = NULL
      ), "'working_units' must be one of the following: 'character'"
    )
  }
)

# 1.3.3. age.w
test_that(
  "prepare_data, smooth_method = 'age.w' throws error for working_units = NULL",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_work_agew <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "age.w",
        smooth_n_points = 5,
        smooth_age_range = 500
      ) %>%
      reduce_data()


    expect_error(
      prepare_data(
        data_source_prep = data_work_agew,
        working_units = NULL,
        bin_size = NULL,
        number_of_shifts = NULL,
        rand = NULL
      ), "'working_units' must be one of the following: 'character'"
    )
  }
)

# 1.3.4. grim
test_that(
  "prepare_data, smooth_method = 'grim' throws error for working_units = NULL",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_work_grim <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9
      ) %>%
      reduce_data()

    expect_error(
      prepare_data(
        data_source_prep = data_work_grim,
        working_units = NULL,
        bin_size = NULL,
        number_of_shifts = NULL,
        rand = NULL
      ), "'working_units' must be one of the following: 'character'"
    )
  }
)

# --------------------------------------------------- #
# 1.4 working_units = invalid character
# 1.4.1 m.avg
test_that(
  "prepare_data, smooth_method = 'm.avg' throws error for working_units = invalid character",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_work_mavg <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "m.avg",
        smooth_n_points = 5
      ) %>%
      reduce_data()

    expect_error(
      prepare_data(
        data_source_prep = data_work_mavg,
        working_units = "invalid_unit",
        bin_size = NULL,
        number_of_shifts = NULL,
        rand = NULL
      ), "'working_units' must contains one of the following values: 'levels', 'bins', 'MW'"
    )
  }
)

# 1.4.2. shep
test_that(
  "prepare_data, smooth_method = 'shep' throws error for working_units = invalid character",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_work_shep <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "shep"
      ) %>%
      reduce_data()

    expect_error(
      prepare_data(
        data_source_prep = data_work_shep,
        working_units = "invalid_unit",
        bin_size = NULL,
        number_of_shifts = NULL,
        rand = NULL
      ), "'working_units' must contains one of the following values: 'levels', 'bins', 'MW'"
    )
  }
)

# 1.4.3. age.w
test_that(
  "prepare_data, smooth_method = 'age.w' throws error for working_units = invalid character",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_work_agew <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "age.w",
        smooth_n_points = 5,
        smooth_age_range = 500
      ) %>%
      reduce_data()


    expect_error(
      prepare_data(
        data_source_prep = data_work_agew,
        working_units = "invalid_unit",
        bin_size = NULL,
        number_of_shifts = NULL,
        rand = NULL
      ), "'working_units' must contains one of the following values: 'levels', 'bins', 'MW'"
    )
  }
)

# 1.4.4. grim
test_that(
  "prepare_data, smooth_method = 'grim' throws error for working_units = invalid character",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_work_grim <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9
      ) %>%
      reduce_data()

    expect_error(
      prepare_data(
        data_source_prep = data_work_grim,
        working_units = "invalid_unit",
        bin_size = NULL,
        number_of_shifts = NULL,
        rand = NULL
      ), "'working_units' must contains one of the following values: 'levels', 'bins', 'MW'"
    )
  }
)

# --------------------------------------------------- #
# 1.5 working_units = multiple
# 1.5.1 m.avg
test_that(
  "prepare_data, smooth_method = 'm.avg' throws error for working_units = multiple values",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_work_mavg <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "m.avg",
        smooth_n_points = 5
      ) %>%
      reduce_data()

    expect_error(
      prepare_data(
        data_source_prep = data_work_mavg,
        working_units = c("levels", "bins"),
        bin_size = 500,
        rand = NULL
      ), "'arg' must be of length 1"
    )
  }
)

# 1.4.2. shep
test_that(
  "prepare_data, smooth_method = 'shep' throws error for working_units = multiple values",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_work_shep <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "shep"
      ) %>%
      reduce_data()

    expect_error(
      prepare_data(
        data_source_prep = data_work_shep,
        working_units = c("levels", "bins"),
        bin_size = 500,
        rand = NULL
      ), "'arg' must be of length 1"
    )
  }
)

# 1.4.3. age.w
test_that(
  "prepare_data, smooth_method = 'age.w' throws error for working_units = multiple values",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_work_agew <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "age.w",
        smooth_n_points = 5,
        smooth_age_range = 500
      ) %>%
      reduce_data()


    expect_error(
      prepare_data(
        data_source_prep = data_work_agew,
        working_units = c("levels", "bins"),
        bin_size = 500,
        rand = NULL
      ), "'arg' must be of length 1"
    )
  }
)

# 1.4.4. grim
test_that(
  "prepare_data, smooth_method = 'grim' throws error for working_units = multiple values",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_work_grim <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9
      ) %>%
      reduce_data()

    expect_error(
      prepare_data(
        data_source_prep = data_work_grim,
        working_units = c("levels", "bins"),
        bin_size = 500,
        rand = NULL
      ), "'arg' must be of length 1"
    )
  }
)

test_that(
  "prepare_data validates working_units parameter is not 0",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = 0
      ),
      "'working_units' must be one of the following: 'character'"
    )
  }
)

test_that(
  "prepare_data validates working_units parameter is not NULL",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = NULL
      ),
      "'working_units' must be one of the following: 'character'"
    )
  }
)

# 1.5 no user-defined input for working_units
## (these tests fail because function will use defaults silently)
test_that(
  "prepare_data throws warning if no user-defined input for working_units is supplied",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_work_mavg <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data() %>%
      reduce_data()

    expect_warning(
      res <-
        prepare_data(
          data_source_prep = data_work_mavg,
          rand = NULL
        ),
      # none programmed into the function yet.
      # e.g., "No user-defined input for working_units. Defaulting to 'levels'."
    )
  }
)

# 1.5.1  m.avg
test_that(
  "prepare_data, smooth_method = 'm.avg' uses 'levels' if no user-defined input for working_units",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_work_mavg <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "m.avg",
        smooth_n_points = 5
      ) %>%
      reduce_data()

    expect_no_error(
      res <-
        prepare_data(
          data_source_prep = data_work_mavg,
          rand = NULL
        )
    )

    # Control: if nothing is user-supplied, uses working_units = "levels"
    res_levels <-
      prepare_data(
        data_source_prep = data_work_mavg,
        working_units = "levels"
      )

    expect_identical(res, res_levels)
  }
)

# 1.5.2. shep
test_that(
  "prepare_data, smooth_method = 'shep' uses 'levels' if no user-defined input for working_units",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_work_shep <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "shep"
      ) %>%
      reduce_data()

    expect_no_error(
      res <-
        prepare_data(
          data_source_prep = data_work_shep,
          rand = NULL
        )
    )
    # Control: if nothing is user-supplied, uses working_units = "levels"
    res_levels <-
      prepare_data(
        data_source_prep = data_work_shep,
        working_units = "levels"
      )

    expect_identical(res, res_levels)
  }
)
# 1.5.3. age.w
test_that(
  "prepare_data, smooth_method = 'age.w' uses 'levels' if no user-defined input for working_units",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_work_agew <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "age.w",
        smooth_n_points = 5,
        smooth_age_range = 500
      ) %>%
      reduce_data()


    expect_no_error(
      res <-
        prepare_data(
          data_source_prep = data_work_agew,
          rand = NULL
        )
    )

    # Control: if nothing is user-supplied, uses working_units = "levels"
    res_levels <-
      prepare_data(
        data_source_prep = data_work_agew,
        working_units = "levels"
      )

    expect_identical(res, res_levels)
  }
)
# 1.5.4. grim
test_that(
  "prepare_data, smooth_method = 'grim' uses 'levels' if no user-defined input for working_units",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_work_grim <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9
      ) %>%
      reduce_data()


    expect_no_error(
      res <-
        prepare_data(
          data_source_prep = data_work_grim,
          rand = NULL
        )
    )

    # Control: if nothing is user-supplied, uses working_units = "levels"
    res_levels <-
      prepare_data(
        data_source_prep = data_work_grim,
        working_units = "levels"
      )

    expect_identical(res, res_levels)
  }
)


### Bin Size ###

# ---------------------------------------------------------- #
#                   Bin_size Input Tests                     #
# ---------------------------------------------------------- #

# 1.1 "levels"
## skipped because bin_size is not required for 'levels'
# 1.2 "bins"
test_that(
  "prepare_data validates bin_size is numeric",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "bins",
        bin_size = "500",
        number_of_shifts = NULL,
        rand = NULL
      ),
      "'bin_size' must be one of the following: 'numeric'"
    )
  }
)

test_that(
  "prepare_data validates bin_size is not Inf",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "bins",
        bin_size = Inf,
        number_of_shifts = NULL,
        rand = NULL
      ),
      "'to' must be a finite number"
    )
  }
)

test_that(
  "prepare_data validates bin_size is not negative",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "bins",
        bin_size = -500,
        number_of_shifts = NULL,
        rand = NULL
      ),
      "wrong sign in 'by' argument"
    )
  }
)

test_that(
  "prepare_data validates bin_size is not NA",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "bins",
        bin_size = NA,
        number_of_shifts = NULL,
        rand = NULL
      ),
      "'bin_size' must be one of the following: 'numeric'"
    )
  }
)

test_that(
  "prepare_data validates bin_size is integer",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )
    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "bins",
        bin_size = 500.5,
        number_of_shifts = NULL,
        rand = NULL
      ),
      "' bin_size ' must be a an integer"
    )
  }
)

test_that(
  "prepare_data validates bin_size is not 0",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )
    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "bins",
        bin_size = 0,
        number_of_shifts = NULL,
        rand = NULL
      ),
      "invalid"
    )
  }
)

test_that(
  "prepare_data works with minimum bin_size (1)",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )
    expect_no_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "bins",
        bin_size = 1,
        number_of_shifts = NULL,
        rand = NULL
      )
    )
  }
)

# 1.3 "MW"
test_that(
  "prepare_data validates bin_size is numeric (working_units=WM)",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "MW",
        bin_size = "500",
        number_of_shifts = 5,
        rand = NULL
      ),
      "'bin_size' must be one of the following: 'numeric'"
    )
  }
)

test_that(
  "prepare_data validates bin_size is not Inf (working_units=WM)",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "MW",
        bin_size = Inf,
        number_of_shifts = 5,
        rand = NULL
      ),
      "'to' must be a finite number"
    )
  }
)

test_that(
  "prepare_data validates bin_size is not negative (working_units=WM)",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "MW",
        bin_size = -500,
        number_of_shifts = 5,
        rand = NULL
      ),
      "wrong sign in 'by' argument"
    )
  }
)

test_that(
  "prepare_data validates bin_size is not NA (working_units=WM)",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "MW",
        bin_size = NA,
        number_of_shifts = 5,
        rand = NULL
      ),
      "'bin_size' must be one of the following: 'numeric'"
    )
  }
)

test_that(
  "prepare_data validates bin_size is integer (working_units=WM)",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )
    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "MW",
        bin_size = 500.5,
        number_of_shifts = 5,
        rand = NULL
      ),
      "' bin_size ' must be a an integer"
    )
  }
)

test_that(
  "prepare_data validates bin_size is not 0 (working_units=WM)",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )
    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "MW",
        bin_size = 0,
        number_of_shifts = 5,
        rand = NULL
      ),
      "invalid"
    )
  }
)

test_that(
  "prepare_data works with minimum bin_size (1) and MW",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )
    expect_no_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "MW",
        bin_size = 1,
        number_of_shifts = 5,
        rand = NULL
      )
    )
  }
)


test_that(
  "prepare_data handles extremely large bin_size",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = FALSE
      )

    age_range <-
      max(example_data$age$age) - min(example_data$age$age)
    huge_bin_size <-
      age_range * 10

    expect_no_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "bins",
        bin_size = huge_bin_size
      )
    )
  }
)

test_that(
  "prepare_data handles small bin_size efficiently",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = FALSE
      )

    # This will create many bins
    expect_no_error(
      result <-
        prepare_data(
          data_source_prep = example_data,
          working_units = "bins",
          bin_size = 1
        )
    )

    # Should still return valid structure
    expect_type(result, "list")
    expect_true(length(result) >= 1)
  }
)

# ---------------------------------------------------------- #
#                Number of shifts Input Tests                #
# ---------------------------------------------------------- #

test_that(
  "prepare_data validates number_of_shifts is present (working_units=WM)",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = NULL,
        rand = NULL
      ),
      "'number_of_shifts' must be one of the following: 'numeric'"
    )
  }
)

test_that(
  "prepare_data validates number_of_shifts is not Inf (working_units=WM)",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = Inf,
        rand = NULL
      ),
      "result would be too long a vector"
    )
  }
)

test_that(
  "prepare_data validates number_of_shifts is not negative (working_units=WM)",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = -5,
        rand = NULL
      ),
      "invalid 'times' argument"
    )
  }
)

test_that(
  "prepare_data validates number_of_shifts is not NA (working_units=WM)",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = NA,
        rand = NULL
      ),
      "'number_of_shifts' must be one of the following: 'numeric'"
    )
  }
)


test_that(
  "prepare_data validates number_of_shifts is integer (working_units=WM)",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )
    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5.5,
        rand = NULL
      ),
      "' number_of_shifts ' must be a an integer"
    )
  }
)

test_that(
  "prepare_data works with number_of_shifts = 0 (working_units=WM)",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )
    expect_no_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 0,
        rand = NULL
      ),
    )
  }
)

test_that(
  "prepare_data handles very large number_of_shifts",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = FALSE
      )

    expect_no_error(
      result <-
        prepare_data(
          data_source_prep = example_data,
          working_units = "MW",
          bin_size = 500,
          number_of_shifts = 100
        )
    )

    # Should create 100 different shifts
    expect_equal(length(result[[1]]), 100)
  }
)

# ---------------------------------------------------------- #
#                 Rand Input Tests                           #
# ---------------------------------------------------------- #

test_that(
  "prepare_data validates rand parameter when provided",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 5,
        rand = "invalid"
      ),
      "'rand' must be one of the following: 'NULL', 'numeric'"
    )
  }
)

test_that(
  "prepare_data validates rand is integer when numeric",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 5,
        rand = 5.5
      ),
      "' rand ' must be a an integer"
    )
  }
)

test_that(
  "prepare_data produces N = rand random samples",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_no_error(
      dat_5 <-
        prepare_data(
          data_source_prep = example_data,
          working_units = "levels",
          bin_size = 500,
          number_of_shifts = 5,
          rand = 5
        )
    )

    dat_3 <-
      prepare_data(
        data_source_prep = example_data,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 5,
        rand = 3
      )

    expect_equal(length(dat_3), 3)
    expect_equal(length(dat_5), 5)
  }
)

test_that(
  "prepare_data produces randomized samples with every trial",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    dat_5 <-
      prepare_data(
        data_source_prep = example_data,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 5,
        rand = 5
      )

    dat_3 <-
      prepare_data(
        data_source_prep = example_data,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 5,
        rand = 3
      )
    expect_false(identical(dat_3, dat_5[1:3]))
  }
)

test_that(
  "prepare_data throws error with rand = 0",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "levels",
        rand = 0
      ),
      "subscript out of bounds"
    )
  }
)


test_that(
  "prepare_data throws error with negative rand",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "levels",
        rand = -5
      ),
      "invalid 'size' argument"
    )
  }
)

test_that(
  "prepare_data handles age_uncertainty with all identical values (returns identical age across randomizations)",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    # Make all age uncertainty values identical
    example_data$age_un[] <-
      1000

    expect_no_error(
      result <-
        prepare_data(
          data_source_prep = example_data,
          working_units = "levels",
          rand = 10
        )
    )

    # All randomizations should be identical
    expect_true(
      all(
        sapply(
          result[-1], function(x) {
            identical(x[[1]]$data$age, result[[1]][[1]]$data$age)
          }
        )
      )
    )
  }
)


# ------------------------------------------ #
#         Reproducibility tests              #
# ------------------------------------------ #

test_that(
  "prepare_data samples reproducibly if rand > 1 and seed set manually",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    set.seed(123)
    res1 <-
      prepare_data(
        data_source_prep = example_data,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 5,
        rand = 5
      )

    set.seed(123)
    res2 <-
      prepare_data(
        data_source_prep = example_data,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 5,
        rand = 5
      )

    expect_identical(res1, res2)
  }
)

# fails currently
test_that(
  "prepare_data samples reproducibly if age_un is supplied, rand = 1 and no seed set manually",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    res1 <-
      prepare_data(
        data_source_prep = example_data,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 5,
        rand = 1
      )

    res2 <-
      prepare_data(
        data_source_prep = example_data,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 5,
        rand = 1
      )

    expect_identical(res1, res2)
  }
)

test_that(
  "prepare_data samples reproducibly if age_un = NULL and rand > 1 and no seed set manually",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = NULL,
        verbose = FALSE
      )

    res1 <-
      prepare_data(
        data_source_prep = example_data,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 5,
        rand = 1
      )

    res2 <-
      prepare_data(
        data_source_prep = example_data,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 5,
        rand = 1
      )

    expect_identical(res1, res2)
  }
)

test_that(
  "prepare_data does not change age if no age_uncertainty in data",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = NULL,
        verbose = FALSE
      )

    res <-
      prepare_data(
        data_source_prep = example_data,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 1,
        rand = 1
      )

    expect_identical(
      res[[1]][[1]]$data$age,
      example_data$age$age
    )
  }
)

test_that(
  "prepare_data changes age if age_uncertainty is in data",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    res <-
      prepare_data(
        data_source_prep = example_data,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 1,
        rand = 1
      )

    expect_false(
      identical(
        res[[1]][[1]]$data$age,
        example_data$age$age
      )
    )
  }
)

test_that(
  "prepare_data changes age if age_uncertainty is in data",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    res <-
      prepare_data(
        data_source_prep = example_data,
        working_units = "levels",
        rand = 1
      )

    expect_false(
      identical(
        res[[1]][[1]]$data$age,
        example_data$age$age
      )
    )
  }
)

test_that(
  "prepare_data changes age if age_uncertainty is in data",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    res <-
      prepare_data(
        data_source_prep = example_data,
        working_units = "bins",
        bin_size = 500,
        rand = 1
      )

    expect_false(
      identical(
        res[[1]][[1]]$data$age,
        example_data$age$age
      )
    )
  }
)


# ---------------------------------------------------------- #
#                  Output structure Tests                    #
# ---------------------------------------------------------- #

# with valid input:
test_that(
  "prepare_data returns consistent structure across different working_units",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = FALSE
      )

    result_levels <-
      prepare_data(
        example_data,
        working_units = "levels"
      )
    result_bins <-
      prepare_data(
        example_data,
        working_units = "bins",
        bin_size = 500
      )
    result_mw <-
      prepare_data(
        example_data,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 3
      )

    # All should have same top-level structure
    expect_identical(
      names(result_levels),
      names(result_bins)
    )
    expect_identical(
      names(result_levels),
      names(result_mw)
    )

    # All should have same internal structure
    expect_identical(
      names(result_levels[[1]][[1]]),
      names(result_bins[[1]][[1]])
    )
    expect_identical(
      names(result_levels[[1]][[1]]),
      names(result_mw[[1]][[1]])
    )

    # Bins should have same column structure
    expect_identical(
      names(result_levels[[1]][[1]]$bins),
      names(result_bins[[1]][[1]]$bins)
    )
    expect_identical(
      names(result_levels[[1]][[1]]$bins), names(result_mw[[1]][[1]]$bins)
    )
  }
)

# Fails because assertion for bin_size is in the wrong position in the function
test_that(
  "prepare_data returns the correct output structure with default parameters (levels)",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    res <-
      prepare_data(
        data_source_prep = example_data,
        working_units = "levels",
        bin_size = NULL,
        number_of_shifts = NULL,
        rand = NULL
      )

    # General output structure tests:
    expect_type(
      res, "list"
    )

    # Named list elements
    expect_identical(
      names(
        res[[1]][[1]]
      ),
      c("data", "bins")
    )

    # Named data object
    expect_identical(
      names(
        res[[1]][[1]]$data
      ),
      c("age", names(example_data$community))
    )

    # Named bins object
    expect_true(
      all(
        c("name", "shift", "age_diff", "start", "end", "res_age", "label")
        %in%
          colnames(res[[1]][[1]]$bins)
      )
    )

    # $data is not all 0 or empty
    expect_false(
      all(
        res[[1]][[1]]$data == 0
      )
    )
    expect_false(
      length(
        res[[1]][[1]]$data
      ) == 0
    )
    # $bins is not all 0 or empty
    expect_false(
      all(
        res[[1]][[1]]$bins == 0
      )
    )
    expect_false(
      nrow(res[[1]][[1]]$bins) == 0
    )
  }
)

test_that(
  "prepare_data returns the correct output structure with default parameters (MW)",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    res <-
      prepare_data(
        data_source_prep = example_data,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      )

    # General output structure tests:
    expect_type(
      res, "list"
    )

    # Named list elements
    expect_identical(
      names(
        res[[1]][[1]]
      ),
      c("data", "bins")
    )

    # Named data object
    expect_identical(
      names(
        res[[1]][[1]]$data
      ),
      c("age", names(example_data$community))
    )

    # Named bins object
    expect_true(
      all(
        c("name", "shift", "age_diff", "start", "end", "res_age", "label")
        %in%
          colnames(res[[1]][[1]]$bins)
      )
    )

    # Correct types
    expect_type(res[[1]][[1]]$bins$name, "double")
    expect_type(res[[1]][[1]]$bins$shift, "integer")
    expect_type(res[[1]][[1]]$bins$age_diff, "double")
    expect_type(res[[1]][[1]]$bins$start, "double")
    expect_type(res[[1]][[1]]$bins$end, "double")
    expect_type(res[[1]][[1]]$bins$res_age, "double")
    expect_type(res[[1]][[1]]$bins$label, "character")


    # $data is not all 0 or empty
    expect_false(all(res[[1]][[1]]$data == 0))
    expect_false(length(res[[1]][[1]]$data) == 0)
    # $bins is not all 0 or empty
    expect_false(all(res[[1]][[1]]$bins == 0))
    expect_false(nrow(res[[1]][[1]]$bins) == 0)
  }
)

test_that(
  "prepare_data returns the correct output structure with default parameters (bins)",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    res <-
      prepare_data(
        data_source_prep = example_data,
        working_units = "bins",
        bin_size = 500,
        rand = NULL
      )

    # General output structure tests:
    expect_type(
      res, "list"
    )

    # Named list elements
    expect_identical(
      names(
        res[[1]][[1]]
      ),
      c("data", "bins")
    )

    # Named data object
    expect_identical(
      names(
        res[[1]][[1]]$data
      ),
      c("age", names(example_data$community))
    )

    # Named bins object
    expect_true(
      all(
        c("name", "shift", "age_diff", "start", "end", "res_age", "label")
        %in%
          colnames(res[[1]][[1]]$bins)
      )
    )

    # Correct types
    expect_type(res[[1]][[1]]$bins$name, "double")
    expect_type(res[[1]][[1]]$bins$shift, "double")
    expect_type(res[[1]][[1]]$bins$age_diff, "double")
    expect_type(res[[1]][[1]]$bins$start, "double")
    expect_type(res[[1]][[1]]$bins$end, "double")
    expect_type(res[[1]][[1]]$bins$res_age, "double")
    expect_type(res[[1]][[1]]$bins$label, "character")

    # $data is not all 0 or empty
    expect_false(all(res[[1]][[1]]$data == 0))
    expect_false(length(res[[1]][[1]]$data) == 0)
    # $bins is not all 0 or empty
    expect_false(all(res[[1]][[1]]$bins == 0))
    expect_false(nrow(res[[1]][[1]]$bins) == 0)
  }
)

# ---------------------------------------------------------- #
#                Integration with make_bins Tests            #
# ---------------------------------------------------------- #

test_that(
  "prepare_data correctly calls make_bins for each working_units type",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = FALSE
      )

    # Test that bins are created correctly for each type
    result_levels <-
      prepare_data(example_data, working_units = "levels")
    result_bins <-
      prepare_data(example_data, working_units = "bins", bin_size = 500)
    result_mw <-
      prepare_data(example_data, working_units = "MW", bin_size = 500, number_of_shifts = 5)

    # Levels should have one bin per sample
    expect_equal(nrow(result_levels[[1]][[1]]$bins), nrow(example_data$community))

    # MW should have multiple shifts
    expect_true(max(result_mw[[1]][[1]]$bins$shift) <= 5)
    expect_true(length(result_mw[[1]]) == 5)

    # Bins should have consistent bin_size
    bins_result <-
      result_bins[[1]][[1]]$bins
    expect_true(all(bins_result$age_diff == 500))
  }
)







# ------------------------------------------------- #
# Tests using internal estimate_roc workflow:       #
#                   Errors                          #
# ------------------------------------------------- #

# Invalid parameter combinations:
test_that(
  "prepare_data throws message if invalid combinations of parameters are silently overridden",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]

    age <-
      RRatepol::example_data$sample_age[[1]]

    data_extract <-
      extract_data(
        community,
        age,
      )

    data_smoothed <-
      smooth_community_data(
        data_extract,
        smooth_method = "age.w"
      )

    data_work <-
      reduce_data(
        data_smoothed
      )

    expect_message(
      data_prepared <-
        prepare_data(
          data_work,
          working_units = "levels",
          number_of_shifts = 5,
          rand = 1
        ),
      # not programmed into the function yet.
      # e.g.,
      # "Invalid parameter 'number_of_shifts = 5' to working_units != 'MW'. number_of_shifts is automatically set to = 1"
    )
  }
)

test_that(
  "prepare_data sets number_of_shifts =1 with user-supplied: 'MW' and number_of_shifts = 0",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]

    age <-
      RRatepol::example_data$sample_age[[1]]

    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_extract <-
      extract_data(
        community,
        age,
        age_uncertainty
      )

    data_smoothed <-
      smooth_community_data(
        data_extract,
        smooth_method = "age.w"
      )

    data_work <-
      reduce_data(
        data_smoothed
      )

    data_prepared <-
      prepare_data(
        data_work,
        working_units = "MW",
        number_of_shifts = 0
      )

    expect_true(
      unique(
        data_prepared[[1]][[1]]$bins$shift
      ) == 1
    )
  }
)

# Test combinations of rand and age_un.
# Q: is it the same whether if the user sets rand = 1 or if the function does it if rand = NULL?
# expect identical results for rand = 1 and rand = NULL.
test_that(
  "Prepare data returns identical results for rand=NULL and rand=1",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]

    age <-
      RRatepol::example_data$sample_age[[1]]

    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_extract <-
      extract_data(
        community,
        age,
        age_uncertainty
      )

    data_smoothed <-
      smooth_community_data(
        data_extract,
        smooth_method = "age.w"
      )

    data_work <-
      reduce_data(
        data_smoothed
      )

    data_prepared_rand_1 <-
      prepare_data(
        data_work,
        working_units = "MW",
        number_of_shifts = 2,
        rand = 1
      )

    data_prepared_rand_null <-
      prepare_data(
        data_work,
        working_units = "MW",
        number_of_shifts = 2,
        rand = NULL
      )

    expect_identical(
      data_prepared_rand_1,
      data_prepared_rand_null
    )
  }
)

test_that(
  "Prepare data returns identical results for rand=NULL and rand=1 (without age_un)",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]

    age <-
      RRatepol::example_data$sample_age[[1]]

    data_extract <-
      extract_data(
        community,
        age,
      )

    data_smoothed <-
      smooth_community_data(
        data_extract,
        smooth_method = "age.w"
      )

    data_work <-
      reduce_data(
        data_smoothed
      )

    data_prepared_rand_1 <-
      prepare_data(
        data_work,
        working_units = "MW",
        number_of_shifts = 2,
        rand = 1
      )

    data_prepared_rand_null <-
      prepare_data(
        data_work,
        working_units = "MW",
        number_of_shifts = 2,
        rand = NULL
      )

    expect_identical(
      data_prepared_rand_1,
      data_prepared_rand_null
    )
  }
)

# Is rand ignored if no age_un? - no. but results are identical.
test_that(
  "prepare_data ignores rand parameter if no age_un is supplied",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]

    age <-
      RRatepol::example_data$sample_age[[1]]

    data_extract <-
      extract_data(
        community,
        age
      )

    data_smoothed <-
      smooth_community_data(
        data_extract,
        smooth_method = "age.w"
      )

    data_work <-
      reduce_data(
        data_smoothed
      )

    data_prepared_rand_1 <-
      prepare_data(
        data_work,
        working_units = "MW",
        number_of_shifts = 1,
        rand = 1
      )

    data_prepared_rand_100 <-
      prepare_data(
        data_work,
        working_units = "MW",
        number_of_shifts = 1,
        rand = 100
      )

    # Expect that rand is ignored if no age_un is supplied:
    expect_true(
      length(data_prepared_rand_100) == length(data_prepared_rand_1)
    )

    # # both results are identical
    # expect_identical(
    #   data_prepared_rand_1[[1]][[1]]$bins,
    #   data_prepared_rand_100[[2]][[1]]$bins
    # )

    # # all "randomized" samples are identical inside the same result
    # expect_identical(
    #   data_prepared_rand_100[[1]][[1]]$bins,
    #   data_prepared_rand_100[[2]][[1]]$bins
    # )

    # expect_identical(
    #   data_prepared_rand_100[[1]][[1]]$bins,
    #   data_prepared_rand_100[[100]][[1]]$bins
    # )
  }
)

test_that(
  "prepare_data returns identical results across randomizations if age_un is not supplied but rand is",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]

    age <-
      RRatepol::example_data$sample_age[[1]]

    data_extract <-
      extract_data(
        community,
        age
      )

    data_smoothed <-
      smooth_community_data(
        data_extract,
        smooth_method = "age.w"
      )

    data_work <-
      reduce_data(
        data_smoothed
      )

    data_prepared_rand_100 <-
      prepare_data(
        data_work,
        working_units = "MW",
        number_of_shifts = 1,
        rand = 100
      )

    # Test if all elements in data_prepared_rand_100 are identical
    expect_true(
      all(
        vapply(
          data_prepared_rand_100[-1],
          function(
              x) {
            identical(
              x, data_prepared_rand_100[[1]]
            )
          },
          logical(1)
        )
      )
    )
  }
)

test_that(
  "prepare_data handles only 1 age_uncertainty row with high rand values - expect identical results across randomizations",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]

    # Test 1: What happens with only 1 row of age_un but rand = 100?
    age_un_1row <-
      RRatepol::example_data$age_uncertainty[[1]][1, 1:63, drop = FALSE]

    data_extract_1row <-
      extract_data(
        community,
        age,
        age_un_1row
      )

    data_smoothed_1row <-
      smooth_community_data(
        data_extract_1row,
        smooth_method = "age.w"
      )

    data_work_1row <-
      reduce_data(
        data_smoothed_1row
      )

    data_prepared_1row <-
      prepare_data(
        data_work_1row,
        working_units = "MW",
        number_of_shifts = 1,
        rand = 100
      )

    # With only 1 row, ALL samples should be identical (sampling with replacement from 1 value)
    expect_equal(
      length(
        data_prepared_1row
      ),
      100
    )

    expect_true(
      all(
        vapply(
          data_prepared_1row[-1],
          function(
              x) {
            identical(
              x[[1]]$data$age, data_prepared_1row[[1]][[1]]$data$age
            )
          },
          logical(1)
        )
      )
    )
  }
)

test_that(
  "prepare_data handles only 2 age_uncertainty rows with high rand values",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    # Test 2: What happens with 2 rows of age_un but rand = 100?
    age_un_2rows <-
      RRatepol::example_data$age_uncertainty[[1]][1:2, 1:63]

    data_extract_2rows <-
      extract_data(
        community,
        age,
        age_un_2rows
      )
    data_smoothed_2rows <-
      smooth_community_data(data_extract_2rows, smooth_method = "age.w")
    data_work_2rows <-
      reduce_data(data_smoothed_2rows)

    data_prepared_2rows <-
      prepare_data(
        data_work_2rows,
        working_units = "MW",
        number_of_shifts = 1,
        rand = 100
      )

    # With 2 rows, we should see some variation (not all identical)
    expect_equal(length(data_prepared_2rows), 100)
    expect_false(
      all(
        vapply(
          data_prepared_2rows[-1], # compare all to the first
          function(x) {
            identical(
              x[[1]]$data$age, data_prepared_2rows[[1]][[1]]$data$age
            )
          },
          logical(1)
        )
      )
    )

    # Test 3: Verify sampling distribution with 2 rows
    # Extract all age vectors and count unique patterns
    age_vectors <-
      lapply(data_prepared_2rows, function(x) x[[1]]$data$age)
    unique_age_patterns <-
      length(unique(age_vectors))

    # Should have at most 2 unique patterns (since only 2 rows available)
    expect_true(unique_age_patterns <= 2)
    expect_true(unique_age_patterns >= 1)
  }
)

test_that(
  "prepare_data: with only 2 age_uncertainty rows and high rand, produces at most 2 unique randomizations, and full uncertainty produces more unique samples",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    # Test 2: What happens with 2 rows of age_un but rand = 100?
    age_un_2rows <-
      RRatepol::example_data$age_uncertainty[[1]][1:2, 1:63]

    data_extract_2rows <-
      extract_data(
        community,
        age,
        age_un_2rows
      )
    data_smoothed_2rows <-
      smooth_community_data(data_extract_2rows, smooth_method = "age.w")
    data_work_2rows <-
      reduce_data(data_smoothed_2rows)

    data_prepared_2rows <-
      prepare_data(
        data_work_2rows,
        working_units = "MW",
        number_of_shifts = 1,
        rand = 100
      )

    # Extract all age vectors and count unique patterns
    age_vectors <-
      lapply(data_prepared_2rows, function(x) x[[1]]$data$age)
    unique_age_patterns <-
      length(unique(age_vectors))

    # Test 4: Control - compare against full uncertainty data
    age_un_full <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_extract_full <-
      extract_data(
        community,
        age,
        age_un_full
      )
    data_smoothed_full <-
      smooth_community_data(
        data_extract_full,
        smooth_method = "age.w"
      )

    data_work_full <-
      reduce_data(
        data_smoothed_full
      )

    data_prepared_full <-
      prepare_data(
        data_work_full,
        working_units = "MW",
        number_of_shifts = 1,
        rand = 100
      )

    # Full uncertainty should have more variation than limited rows
    age_vectors_full <-
      lapply(
        data_prepared_full,
        function(
            x) {
          x[[1]]$data$age
        }
      )
    unique_patterns_full <-
      length(
        unique(
          age_vectors_full
        )
      )

    expect_true(
      unique_patterns_full > unique_age_patterns
    )
  }
)

test_that(
  "prepare_data with seed produces reproducible results with rand = 100",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]

    # Control - compare against full uncertainty data
    age_un_full <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_extract_full <-
      extract_data(
        community,
        age,
        age_un_full
      )
    data_smoothed_full <-
      smooth_community_data(
        data_extract_full,
        smooth_method = "age.w"
      )
    data_work_full <-
      reduce_data(
        data_smoothed_full
      )
    # Test 5: Verify first sample is identical regardless of uncertainty data size
    # (when using same random seed, first sample should use same row)
    set.seed(123)
    data_prep_full_seeded1 <-
      prepare_data(
        data_work_full,
        working_units = "MW",
        number_of_shifts = 1,
        rand = 100
      )
    set.seed(123)
    data_prep_full_seeded2 <-
      prepare_data(
        data_work_full,
        working_units = "MW",
        number_of_shifts = 1,
        rand = 100
      )

    # First random sample should be from same row index
    expect_identical(
      data_prep_full_seeded1[[1]][[1]]$data$age,
      data_prep_full_seeded2[[1]][[1]]$data$age
    )
  }
)

# -------------------------------------------------------------------- #
# 2. Test edge-cases with zeros and NAs in the raw data
# full estimate_roc() workflow:
# -------------------------------------------------------------------- #
# 2.1 With 0 or NA in samples
## 2.1.1a) one row in community 0
# ---- levels ---- #
test_that(
  "prepare_data, working_units = 'levels', smooth_method = 'm.avg' works with one row 0 in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      0
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_mavg <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "m.avg",
        smooth_n_points = 5
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_mavg,
        working_units = "levels",
        rand = NULL
      )

    # test if result$bins has only valid samples
    expect_identical(
      res[[1]][[1]]$bins$name,
      valid_samples
    )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data, working_units = 'levels', smooth_method = 'shep' works with one row 0 in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      0
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_shep <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "shep",
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_shep,
        working_units = "levels",
        rand = NULL
      )

    # test if result$bins has only valid samples
    expect_identical(
      res[[1]][[1]]$bins$name,
      valid_samples
    )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data, working_units = 'levels', smooth_method = 'age.w' works with one row 0 in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      0
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_agew <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "age.w",
        smooth_n_points = 5,
        smooth_age_range = 500
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_agew,
        working_units = "levels",
        rand = NULL
      )

    # test if result$bins has only valid samples
    expect_identical(
      res[[1]][[1]]$bins$name,
      valid_samples
    )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data, working_units = 'levels', smooth_method = 'grim' works with one row 0 in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      0
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_grim <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_grim,
        working_units = "levels",
        rand = NULL
      )

    # test if result$bins has only valid samples
    expect_identical(
      res[[1]][[1]]$bins$name,
      valid_samples
    )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

# --- bins --- #
test_that(
  "prepare_data, working_units = 'bins', smooth_method = 'm.avg' works with one row 0 in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      0
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_mavg <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "m.avg",
        smooth_n_points = 5
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_mavg,
        working_units = "bins",
        bin_size = 500,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data, working_units = 'bins', smooth_method = 'shep' works with one row 0 in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      0
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_shep <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "shep",
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_shep,
        working_units = "bins",
        bin_size = 500,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data, working_units = 'bins', smooth_method = 'age.w' works with one row 0 in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      0
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_agew <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "age.w",
        smooth_n_points = 5,
        smooth_age_range = 500
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_agew,
        working_units = "bins",
        bin_size = 500,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data, working_units = 'bins', smooth_method = 'grim' works with one row 0 in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      0
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_grim <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_grim,
        working_units = "bins",
        bin_size = 500,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

# --- MW --- #
test_that(
  "prepare_data, working_units = 'MW', smooth_method = 'm.avg' works with one row 0 in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      0
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_mavg <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "m.avg",
        smooth_n_points = 5
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_mavg,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data, working_units = 'MW', smooth_method = 'shep' works with one row 0 in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      0
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_shep <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "shep",
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_shep,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data, working_units = 'MW', smooth_method = 'age.w' works with one row 0 in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      0
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_agew <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "age.w",
        smooth_n_points = 5,
        smooth_age_range = 500
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_agew,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data, working_units = 'MW', smooth_method = 'grim' works with one row 0 in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      0
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_grim <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_grim,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

## 2.1.1b) one row in community NA
test_that(
  "prepare_data, working_units = 'levels', smooth_method = 'm.avg' works with one row NA in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      NA
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_mavg <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "m.avg",
        smooth_n_points = 5
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_mavg,
        working_units = "levels",
        rand = NULL
      )

    # test if result$bins has only valid samples
    expect_identical(
      res[[1]][[1]]$bins$name,
      valid_samples
    )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)


test_that(
  "prepare_data, working_units = 'levels', smooth_method = 'shep' works with one row NA in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      NA
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_shep <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "shep",
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_shep,
        working_units = "levels",
        rand = NULL
      )

    # test if result$bins has only valid samples
    expect_identical(
      res[[1]][[1]]$bins$name,
      valid_samples
    )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data and working_units = 'levels' works with one row NA in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      NA
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_agew <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "age.w",
        smooth_n_points = 5,
        smooth_age_range = 500
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_agew,
        working_units = "levels",
        rand = NULL
      )

    # test if result$bins has only valid samples
    expect_identical(
      res[[1]][[1]]$bins$name,
      valid_samples
    )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data and working_units = 'levels' works with one row NA in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      NA
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_grim <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_grim,
        working_units = "levels",
        rand = NULL
      )

    # test if result$bins has only valid samples
    expect_identical(
      res[[1]][[1]]$bins$name,
      valid_samples
    )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

# --- bins --- #
test_that(
  "prepare_data, working_units = 'bins', smooth_method = 'm.avg' works with one row NA in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      NA
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_mavg <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "m.avg",
        smooth_n_points = 5
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_mavg,
        working_units = "bins",
        bin_size = 500,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data, working_units = 'bins', smooth_method = 'shep' works with one row NA in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      NA
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_shep <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "shep",
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_shep,
        working_units = "bins",
        bin_size = 500,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data, working_units = 'bins', smooth_method = 'age.w' works with one row NA in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      NA
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_agew <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "age.w",
        smooth_n_points = 5,
        smooth_age_range = 500
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_agew,
        working_units = "bins",
        bin_size = 500,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data, working_units = 'bins', smooth_method = 'grim' works with one row NA in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      NA
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_grim <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_grim,
        working_units = "bins",
        bin_size = 500,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

# --- MW --- #
test_that(
  "prepare_data, working_units = 'MW' works with one row NA in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      NA
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_mavg <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "m.avg",
        smooth_n_points = 5
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_mavg,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data, working_units = 'MW' works with one row NA in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      NA
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_shep <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "shep",
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_shep,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data, working_units = 'MW' works with one row NA in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      NA
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_agew <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "age.w",
        smooth_n_points = 5,
        smooth_age_range = 500
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_agew,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data, working_units = 'MW' works with one row NA in community",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    community[1, -1] <-
      NA
    valid_samples <-
      community[-1, 1]$sample_id

    data_work_grim <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_grim,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

## 2.1.2a) one row in age 0
# -> no effect or reaction from function. Is just ignored.

### ------ Review required ------- ###

## Below are some tests where I manually fix some bugs
## from upstream functions within each unit test and compare the
## results to the expected results after fixing the bug.

# To Do: Double-check below (NAs in age.)

## 2.1.2b) one row in age NA
# ---- levels ---- #
test_that(
  "prepare_data, working_units = 'levels', smooth_method = 'm.avg' works with one row NA in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    age$age[1] <-
      NA

    data_work_mavg <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "m.avg",
        smooth_n_points = 5
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_mavg,
        working_units = "levels",
        rand = NULL
      )

    # Control for when NA row is exculded from all data sources:
    age_dropped <-
      age %>%
      na.omit()

    community_dropped <-
      community[-1, ]

    age_un_dropped <-
      age_uncertainty[, -1]

    data_dropped <-
      extract_data(
        data_community_extract = community_dropped,
        data_age_extract = age_dropped,
        age_uncertainty = age_un_dropped
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "age.w",
        smooth_n_points = 5,
        smooth_age_range = 500
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res_dropped <-
      prepare_data(
        data_source_prep = data_dropped,
        working_units = "levels",
        rand = NULL
      )

    expect_identical(
      res,
      res_dropped
    )
  }
)

test_that(
  "prepare_data, working_units = 'levels', smooth_method = 'shep' works with one row NA in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    age$age[1] <-
      NA

    data_work_shep <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "shep",
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_shep,
        working_units = "levels",
        rand = NULL
      )

    # Control for when NA row is exculded from all data sources:
    age_dropped <-
      age %>%
      na.omit()
    community_dropped <-
      community[-1, ]
    age_un_dropped <-
      age_uncertainty[, -1]

    data_dropped <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "age.w",
        smooth_n_points = 5,
        smooth_age_range = 500
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res_dropped <-
      prepare_data(
        data_source_prep = data_dropped,
        working_units = "levels",
        rand = NULL
      )

    expect_identical(
      res,
      res_dropped
    )
  }
)

test_that(
  "prepare_data and working_units = 'levels' works with one row NA in age (correct results)",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    age$age[1] <-
      NA

    data_work_agew <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "age.w",
        smooth_n_points = 5,
        smooth_age_range = 500
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res <-
      prepare_data(
        data_source_prep = data_work_agew,
        working_units = "levels",
        rand = NULL
      )

    # Control for when NA row is exculded from all data sources:
    age_dropped <-
      age %>%
      na.omit()
    community_dropped <-
      community[-1, ]
    age_un_dropped <-
      age_uncertainty[, -1]

    data_dropped <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "age.w",
        smooth_n_points = 5,
        smooth_age_range = 500
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    res_dropped <-
      prepare_data(
        data_source_prep = data_dropped,
        working_units = "levels",
        rand = NULL
      )

    expect_identical(
      res,
      res_dropped
    )
  }
)

# fails because smooth_community_data cannot handle NA
test_that(
  "prepare_data and working_units = 'levels' works with one row NA in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    age$age[1] <-
      NA

    data_work_grim <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    valid_samples <-
      data_work_grim$age %>%
      na.omit() %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      pull(sample_id)

    # Simulating to fix the bugs in extract_data() and reduce_data():
    # drop NAs from age manually
    data_work_grim$age <-
      data_work_grim$age %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      filter(sample_id %in% valid_samples) %>%
      select(-sample_id)

    # perform two-way matching with community and age_un manually
    data_work_grim$community <-
      data_work_grim$community %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      filter(sample_id %in% valid_samples) %>%
      select(-sample_id)

    data_work_grim$age_un <-
      data_work_grim$age_un %>%
      select(valid_samples)

    res <-
      prepare_data(
        data_source_prep = data_work_grim,
        working_units = "levels",
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )

    # expect no more NAs in age
    expect_false(
      any(
        is.na(
          res[[1]][[1]]$data$age
        )
      )
    )

    # expect no NAs in age-derivates
    expect_false(
      any(
        is.na(
          res[[1]][[1]]$bin$age_diff
        )
      )
    )

    # expect no NAs in age-derivates
    expect_false(
      any(
        is.na(
          res[[1]][[1]]$bin$res_age
        )
      )
    )
  }
)


# --- bins --- #
test_that(
  "prepare_data, working_units = 'bins', smooth_method = 'm.avg' works with one row NA in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    age$age[1] <-
      NA

    data_work_mavg <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "m.avg",
        smooth_n_points = 5
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    valid_samples <-
      data_work_mavg$age %>%
      na.omit() %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      pull(sample_id)

    # Simulating to fix the bugs in extract_data() and reduce_data():
    # drop NAs from age manually
    data_work_mavg$age <-
      data_work_mavg$age %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      filter(sample_id %in% valid_samples) %>%
      select(-sample_id)

    # perform two-way matching with community and age_un manually
    data_work_mavg$community <-
      data_work_mavg$community %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      filter(sample_id %in% valid_samples) %>%
      select(-sample_id)

    data_work_mavg$age_un <-
      data_work_mavg$age_un %>%
      select(valid_samples)

    # use modified result for prepare_data()
    res <-
      prepare_data(
        data_source_prep = data_work_mavg,
        working_units = "bins",
        bin_size = 500,
        rand = NULL
      )

    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data, working_units = 'bins', smooth_method = 'shep' works with one row NA in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    age$age[1] <-
      NA

    data_work_shep <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "shep",
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    valid_samples <-
      data_work_shep$age %>%
      na.omit() %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      pull(sample_id)

    # Simulating to fix the bugs in extract_data() and reduce_data():
    # drop NAs from age manually
    data_work_shep$age <-
      data_work_shep$age %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      filter(sample_id %in% valid_samples) %>%
      select(-sample_id)

    # perform two-way matching with community and age_un manually
    data_work_shep$community <-
      data_work_shep$community %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      filter(sample_id %in% valid_samples) %>%
      select(-sample_id)

    data_work_shep$age_un <-
      data_work_shep$age_un %>%
      select(valid_samples)

    res <-
      prepare_data(
        data_source_prep = data_work_shep,
        working_units = "bins",
        bin_size = 500,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data, working_units = 'bins', smooth_method = 'age.w' works with one row NA in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    age$age[1] <-
      NA

    data_work_agew <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "age.w",
        smooth_n_points = 5,
        smooth_age_range = 500
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    valid_samples <-
      data_work_agew$age %>%
      na.omit() %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      pull(sample_id)

    # Simulating to fix the bugs in extract_data() and reduce_data():
    # drop NAs from age manually
    data_work_agew$age <-
      data_work_agew$age %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      filter(sample_id %in% valid_samples) %>%
      select(-sample_id)

    # perform two-way matching with community and age_un manually
    data_work_agew$community <-
      data_work_agew$community %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      filter(sample_id %in% valid_samples) %>%
      select(-sample_id)

    data_work_agew$age_un <-
      data_work_agew$age_un %>%
      select(valid_samples)

    res <-
      prepare_data(
        data_source_prep = data_work_agew,
        working_units = "bins",
        bin_size = 500,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data, working_units = 'bins', smooth_method = 'grim' works with one row NA in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    age$age[1] <-
      NA

    data_work_grim <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    valid_samples <-
      data_work_grim$age %>%
      na.omit() %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      pull(sample_id)

    # Simulating to fix the bugs in extract_data() and reduce_data():
    # drop NAs from age manually
    data_work_grim$age <-
      data_work_grim$age %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      filter(sample_id %in% valid_samples) %>%
      select(-sample_id)

    # perform two-way matching with community and age_un manually
    data_work_grim$community <-
      data_work_grim$community %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      filter(sample_id %in% valid_samples) %>%
      select(-sample_id)

    data_work_grim$age_un <-
      data_work_grim$age_un %>%
      select(valid_samples)

    res <-
      prepare_data(
        data_source_prep = data_work_grim,
        working_units = "bins",
        bin_size = 500,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

# --- MW --- #

test_that(
  "prepare_data and working_units = 'MW' works with one row NA in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    age$age[1] <-
      NA
    valid_samples <-
      age %>%
      na.omit() %>%
      pull(sample_id)

    data_work_mavg <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "m.avg",
        smooth_n_points = 5
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    valid_samples <-
      data_work_mavg$age %>%
      na.omit() %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      pull(sample_id)

    # Simulating to fix the bugs in extract_data() and reduce_data():
    # drop NAs from age manually
    data_work_mavg$age <-
      data_work_mavg$age %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      filter(sample_id %in% valid_samples) %>%
      select(-sample_id)

    # perform two-way matching with community and age_un manually
    data_work_mavg$community <-
      data_work_mavg$community %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      filter(sample_id %in% valid_samples) %>%
      select(-sample_id)

    data_work_mavg$age_un <-
      data_work_mavg$age_un %>%
      select(valid_samples)

    res <-
      prepare_data(
        data_source_prep = data_work_mavg,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data, working_units = 'MW' works with one row NA in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    age$age[1] <-
      NA

    data_work_shep <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "shep",
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    valid_samples <-
      data_work_shep$age %>%
      na.omit() %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      pull(sample_id)

    # Simulating to fix the bugs in extract_data() and reduce_data():
    # drop NAs from age manually
    data_work_shep$age <-
      data_work_shep$age %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      filter(sample_id %in% valid_samples) %>%
      select(-sample_id)

    # perform two-way matching with community and age_un manually
    data_work_shep$community <-
      data_work_shep$community %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      filter(sample_id %in% valid_samples) %>%
      select(-sample_id)

    data_work_shep$age_un <-
      data_work_shep$age_un %>%
      select(valid_samples)

    res <-
      prepare_data(
        data_source_prep = data_work_shep,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data, working_units = 'MW' works with one row NA in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    age$age[1] <-
      NA
    valid_samples <-
      age %>%
      na.omit() %>%
      pull(sample_id)

    data_work_agew <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "age.w",
        smooth_n_points = 5,
        smooth_age_range = 500
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    valid_samples <-
      data_work_agew$age %>%
      na.omit() %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      pull(sample_id)

    # Simulating to fix the bugs in extract_data() and reduce_data():
    # drop NAs from age manually
    data_work_agew$age <-
      data_work_agew$age %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      filter(sample_id %in% valid_samples) %>%
      select(-sample_id)

    # perform two-way matching with community and age_un manually
    data_work_agew$community <-
      data_work_agew$community %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      filter(sample_id %in% valid_samples) %>%
      select(-sample_id)

    data_work_agew$age_un <-
      data_work_agew$age_un %>%
      select(valid_samples)

    res <-
      prepare_data(
        data_source_prep = data_work_agew,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

test_that(
  "prepare_data, working_units = 'MW' works with one row NA in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    age$age[1] <-
      NA

    data_work_grim <-
      extract_data(
        data_community_extract = community,
        data_age_extract = age,
        age_uncertainty = age_uncertainty
      ) %>%
      smooth_community_data(
        .,
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9
      ) %>%
      reduce_data(
        check_taxa = TRUE,
        check_levels = TRUE
      )

    valid_samples <-
      data_work_grim$age %>%
      na.omit() %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      pull(sample_id)

    # Simulating to fix the bugs in extract_data() and reduce_data():
    # drop NAs from age manually
    data_work_grim$age <-
      data_work_grim$age %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      filter(sample_id %in% valid_samples) %>%
      select(-sample_id)

    # perform two-way matching with community and age_un manually
    data_work_grim$community <-
      data_work_grim$community %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      filter(sample_id %in% valid_samples) %>%
      select(-sample_id)

    data_work_grim$age_un <-
      data_work_grim$age_un %>%
      select(valid_samples)

    res <-
      prepare_data(
        data_source_prep = data_work_grim,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5,
        rand = NULL
      )

    # test if result$data has only valid samples
    expect_identical(
      rownames(
        res[[1]][[1]]$data
      ),
      valid_samples
    )
  }
)

# ---------------- Check for NAs in output ------------------ #
# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]

    age$age[1:10] <-
      NA
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_extract <-
      extract_data(
        community,
        age,
        age_uncertainty
      )

    data_smoothed <-
      smooth_community_data(
        data_extract,
        smooth_method = "shep"
      )

    data_work <-
      reduce_data(
        data_smoothed
      )

    # empty working_units uses 'levels'
    data_prepared <-
      prepare_data(
        data_work
      )

    # expect no more NAs in bins
    expect_false(
      isTRUE(
        any(
          is.na(
            data_prepared[[1]][[1]]$bins
          )
        )
      )
    )
  }
)

# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]

    age$age[1:10] <-
      NA
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_extract <-
      extract_data(
        community,
        age,
        age_uncertainty
      )

    data_smoothed <-
      smooth_community_data(
        data_extract,
        smooth_method = "m.avg"
      )

    data_work <-
      reduce_data(
        data_smoothed
      )

    # empty working_units uses 'levels'
    data_prepared <-
      prepare_data(
        data_work
      )

    # expect no more NAs in bins
    expect_false(
      isTRUE(
        any(
          is.na(
            data_prepared[[1]][[1]]$bins
          )
        )
      )
    )
  }
)


# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]

    age$age[1:10] <-
      NA
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_extract <-
      extract_data(
        community,
        age,
        age_uncertainty
      )

    data_smoothed <-
      smooth_community_data(
        data_extract,
        smooth_method = "grim"
      )

    data_work <-
      reduce_data(
        data_smoothed
      )

    # empty working_units uses 'levels'
    data_prepared <-
      prepare_data(
        data_work
      )

    # expect no more NAs in bins
    expect_false(
      isTRUE(
        any(
          is.na(
            data_prepared[[1]][[1]]$bins
          )
        )
      )
    )
  }
)

# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]

    age$age[1:10] <-
      NA
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_extract <-
      extract_data(
        community,
        age,
        age_uncertainty
      )

    data_smoothed <-
      smooth_community_data(
        data_extract,
        smooth_method = "age.w"
      )

    data_work <-
      reduce_data(
        data_smoothed
      )

    # empty working_units uses 'levels'
    data_prepared <-
      prepare_data(
        data_work
      )

    # expect no more NAs in bins
    expect_false(
      isTRUE(
        any(
          is.na(
            data_prepared[[1]][[1]]$bins
          )
        )
      )
    )
  }
)


# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]

    age$age[1:10] <-
      NA
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_extract <-
      extract_data(
        community,
        age,
        age_uncertainty
      )

    data_smoothed <-
      smooth_community_data(
        data_extract,
        smooth_method = "shep"
      )

    data_work <-
      reduce_data(
        data_smoothed
      )

    data_prepared <-
      prepare_data(
        data_work,
        working_units = "bins"
      )

    # expect no more NAs in bins
    expect_false(
      isTRUE(
        any(
          is.na(
            data_prepared[[1]][[1]]$bins
          )
        )
      )
    )
  }
)

# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]

    age$age[1:10] <-
      NA
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_extract <-
      extract_data(
        community,
        age,
        age_uncertainty
      )

    data_smoothed <-
      smooth_community_data(
        data_extract,
        smooth_method = "m.avg"
      )

    data_work <-
      reduce_data(
        data_smoothed
      )

    data_prepared <-
      prepare_data(
        data_work,
        working_units = "bins"
      )

    # expect no more NAs in bins
    expect_false(
      isTRUE(
        any(
          is.na(
            data_prepared[[1]][[1]]$bins
          )
        )
      )
    )
  }
)


# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]

    age$age[1:10] <-
      NA
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_extract <-
      extract_data(
        community,
        age,
        age_uncertainty
      )

    data_smoothed <-
      smooth_community_data(
        data_extract,
        smooth_method = "grim"
      )

    data_work <-
      reduce_data(
        data_smoothed
      )

    data_prepared <-
      prepare_data(
        data_work,
        working_units = "bins"
      )

    # expect no more NAs in bins
    expect_false(
      isTRUE(
        any(
          is.na(
            data_prepared[[1]][[1]]$bins
          )
        )
      )
    )
  }
)

# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]

    age$age[1:10] <-
      NA
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_extract <-
      extract_data(
        community,
        age,
        age_uncertainty
      )

    data_smoothed <-
      smooth_community_data(
        data_extract,
        smooth_method = "age.w"
      )

    data_work <-
      reduce_data(
        data_smoothed
      )

    data_prepared <-
      prepare_data(
        data_work,
        working_units = "bins"
      )

    # expect no more NAs in bins
    expect_false(
      isTRUE(
        any(
          is.na(
            data_prepared[[1]][[1]]$bins
          )
        )
      )
    )
  }
)

# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]

    age$age[1:10] <-
      NA
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_extract <-
      extract_data(
        community,
        age,
        age_uncertainty
      )

    data_smoothed <-
      smooth_community_data(
        data_extract,
        smooth_method = "shep"
      )

    data_work <-
      reduce_data(
        data_smoothed
      )

    data_prepared <-
      prepare_data(
        data_work,
        working_units = "MW"
      )

    # expect no more NAs in bins
    expect_false(
      isTRUE(
        any(
          is.na(
            data_prepared[[1]][[1]]$bins
          )
        )
      )
    )
  }
)

# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]

    age$age[1:10] <-
      NA
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_extract <-
      extract_data(
        community,
        age,
        age_uncertainty
      )

    data_smoothed <-
      smooth_community_data(
        data_extract,
        smooth_method = "m.avg"
      )

    data_work <-
      reduce_data(
        data_smoothed
      )

    data_prepared <-
      prepare_data(
        data_work,
        working_units = "MW"
      )

    # expect no more NAs in bins
    expect_false(
      isTRUE(
        any(
          is.na(
            data_prepared[[1]][[1]]$bins
          )
        )
      )
    )
  }
)


# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]

    age$age[1:10] <-
      NA
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_extract <-
      extract_data(
        community,
        age,
        age_uncertainty
      )

    data_smoothed <-
      smooth_community_data(
        data_extract,
        smooth_method = "grim"
      )

    data_work <-
      reduce_data(
        data_smoothed
      )

    data_prepared <-
      prepare_data(
        data_work,
        working_units = "MW"
      )

    # expect no more NAs in bins
    expect_false(
      isTRUE(
        any(
          is.na(
            data_prepared[[1]][[1]]$bins
          )
        )
      )
    )
  }
)

# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age",
  {
    community <-
      RRatepol::example_data$pollen_data[[1]]
    age <-
      RRatepol::example_data$sample_age[[1]]

    age$age[1:10] <-
      NA
    age_uncertainty <-
      RRatepol::example_data$age_uncertainty[[1]]

    data_extract <-
      extract_data(
        community,
        age,
        age_uncertainty
      )

    data_smoothed <-
      smooth_community_data(
        data_extract,
        smooth_method = "age.w"
      )

    data_work <-
      reduce_data(
        data_smoothed
      )

    data_prepared <-
      prepare_data(
        data_work,
        working_units = "MW"
      )

    # expect no more NAs in bins
    expect_false(
      isTRUE(
        any(
          is.na(
            data_prepared[[1]][[1]]$bins
          )
        )
      )
    )
  }
)

#####################
## To Dos:
# - Double check if the expectation shouldn't be the reverse below for age_un
#####################

## 2.1.3a) one column in age_un 0
# ## Age_uncertainty data
# # Zeros
# test_that(
#   "prepare_data with working_units='levels' and rand=1 works if there are 0s in the age uncertainty",
#   {
#     example_data <-
#       extract_data(
#         data_community_extract = RRatepol::example_data$pollen_data[[1]],
#         data_age_extract = RRatepol::example_data$sample_age[[1]],
#         age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
#         verbose = FALSE
#       )

#     example_data$age_un[1:3] <-
#       0

#     expect_no_error(
#       res <-
#         prepare_data(
#           example_data,
#           working_units = "levels",
#           rand = 1
#         )
#     )

#     expect_false(
#       identical(
#         res[[1]][[1]]$data$age,
#         example_data$age$age
#       )
#     )
#   }
# )

# test_that(
#   "prepare_data with working_units='bins' and rand=1 works if there are 0s in the age uncertainty",
#   {
#     example_data <-
#       extract_data(
#         data_community_extract = RRatepol::example_data$pollen_data[[1]],
#         data_age_extract = RRatepol::example_data$sample_age[[1]],
#         age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
#         verbose = FALSE
#       )

#     example_data$age_un[1:3] <-
#       0

#     expect_no_error(
#       res <-
#         prepare_data(
#           example_data,
#           working_units = "bins",
#           bin_size = 500,
#           rand = 1
#         )
#     )

#     expect_false(
#       identical(
#         res[[1]][[1]]$data$age,
#         example_data$age$age
#       )
#     )
#   }
# )


# test_that(
#   "prepare_data with working_units='MW' and rand=1 works if there are 0s in the age uncertainty",
#   {
#     example_data <-
#       extract_data(
#         data_community_extract = RRatepol::example_data$pollen_data[[1]],
#         data_age_extract = RRatepol::example_data$sample_age[[1]],
#         age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
#         verbose = FALSE
#       )

#     example_data$age_un[1:3] <-
#       0

#     expect_no_error(
#       res <-
#         prepare_data(
#           example_data,
#           working_units = "MW",
#           bin_size = 500,
#           number_of_shifts = 5,
#           rand = 1
#         )
#     )

#     expect_false(
#       identical(
#         res[[1]][[1]]$data$age,
#         example_data$age$age
#       )
#     )

#     expect_identical(
#       res[[1]][[1]]$data$age[1:3],
#       c(0, 0, 0)
#     )
#   }
# )


## 2.1.3b) one column in age_un NA
# # NAs in Age_uncertainty
# test_that(
#   "prepare_data with working_units='levels' works if there are NAs in age_un data",
#   {
#     example_data <-
#       extract_data(
#         data_community_extract = RRatepol::example_data$pollen_data[[1]],
#         data_age_extract = RRatepol::example_data$sample_age[[1]],
#         age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
#         verbose = FALSE
#       )

#     example_data$age_un[1:3] <-
#       NA

#     expect_no_error(
#       res <-
#         prepare_data(
#           example_data,
#           working_units = "levels",
#           rand = 1
#         )
#     )

#     expect_identical(
#       res[[1]][[1]]$data$age[1:3],
#       as.numeric(c(NA, NA, NA))
#     )
#   }
# )

# test_that(
#   "prepare_data with working_units='bins' works if there are NAs in age_un data",
#   {
#     example_data <-
#       extract_data(
#         data_community_extract = RRatepol::example_data$pollen_data[[1]],
#         data_age_extract = RRatepol::example_data$sample_age[[1]],
#         age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
#         verbose = FALSE
#       )

#     example_data$age_un[1:3] <-
#       NA

#     expect_no_error(
#       prepare_data(
#         example_data,
#         working_units = "bins",
#         bin_size = 500,
#         rand = 1
#       )
#     )
#   }
# )

# test_that(
#   "prepare_data with working_units='MW' returns NA in data if there is NA in age_uncertainty",
#   {
#     example_data <-
#       extract_data(
#         data_community_extract = RRatepol::example_data$pollen_data[[1]],
#         data_age_extract = RRatepol::example_data$sample_age[[1]],
#         age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
#         verbose = FALSE
#       )

#     example_data$age_un[1:3] <-
#       NA

#     expect_no_error(
#       res <-
#         prepare_data(
#           example_data,
#           working_units = "MW",
#           bin_size = 500,
#           number_of_shifts = 5,
#           rand = 1
#         )
#     )

#     # returns NA in output$age
#     expect_true(any(is.na(res[[1]][[1]]$data$age)))
#   }
# )
