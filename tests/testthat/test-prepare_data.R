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
  "prepare_data rejects incomplete list(data.frame) in data_source",
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
  "prepare_data throws error with data_source_prep with wrong structure",
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

# no error without age_uncertainty
test_that(
  "prepare_data works if no age_uncertainty in data_source",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = NULL,
        verbose = FALSE
      )

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

    expect_equal(
      nrow(
        result[[1]][[1]]$bins), 
      1
      )
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
      ncol(result[[1]][[1]]$data) - 1, # -1 for age column
      1
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
  "prepare_data validates working_units parameter is not NA",
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
# ---------------------------------------------------------- #
#                   Bin_size Input Tests                     #
# ---------------------------------------------------------- #

# 1.1 "levels"
test_that(
  "prepare_data takes bin_size=NULL if working_units = levels",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )
    
    expect_silent(
      prepare_data(
        data_source_prep = example_data,
        working_units = "levels",
        bin_size = NULL,
        number_of_shifts = NULL,
        rand = NULL
      )
    )
  }
)

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
  "prepare_data throws warning if age_uncertainty has all identical values",
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

    expect_warning(
      result <-
        prepare_data(
          data_source_prep = example_data,
          working_units = "levels",
          rand = 10
        ),
      #none programmed into the function yet
      # e.g., 
      #"Warning: age_uncertainty has only 1 unique value. 
      # randomizations will result in identical results"
    )

    # All randomizations should be identical
    # (returns identical age across randomizations)
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
  "prepare_data changes age if age_uncertainty is in data and working_units = MW",
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
  "prepare_data changes age if age_uncertainty is in data and working_units = levels",
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
  "prepare_data changes age if age_uncertainty is in data and working_units = bins",
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
        #bin_size = NULL,
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
    
    set.seed(123)
    data_prepared_rand_1 <-
      prepare_data(
        data_work,
        working_units = "MW",
        number_of_shifts = 2,
        rand = 1
      )

    set.seed(123)
    data_prepared_rand_null <-
      prepare_data(
        data_work,
        working_units = "MW",
        number_of_shifts = 2,
        rand = NULL
      )

    expect_true(identical(
      data_prepared_rand_1,
      data_prepared_rand_null
    ))
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
  "prepare_data throws warning if rand parameter is supplied but not age_un",
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

    expect_warning(
    data_prepared_rand_100 <-
      prepare_data(
        data_work,
        working_units = "MW",
        number_of_shifts = 1,
        rand = 100
      ),
    #none programmed into function yet
    # e.g., 
    # "Warning: setting rand != NULL without age_un will result 
    # in identical randomization results."
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

