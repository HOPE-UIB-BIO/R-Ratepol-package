# ---------------------------------------------------------- #
#               Data_source_prep Input Tests                 #
# ---------------------------------------------------------- #

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




## Community data
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

# NAs in Age
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

# NAs in Age
test_that(
  "prepare_data with working_units='levels' works if there are NAs in age data",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    example_data$age$age[1:3] <-
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
  "prepare_data with working_units='bins' fails if there are NAs in age data",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    example_data$age$age[1:3] <-
      NA

    expect_error(
      prepare_data(
        example_data,
        working_units = "bins",
        bin_size = 500
      ),
      "'from' must be a finite number"
    )
  }
)


test_that(
  "prepare_data with working_units='MW' fails if there are NAs in age data",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    example_data$age$age[1:3] <-
      NA

    expect_error(
      prepare_data(
        example_data,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5
      ),
      "'from' must be a finite number"
    )
  }
)


## Age_uncertainty data
# Zeros
test_that(
  "prepare_data with working_units='levels' and rand=1 works if there are 0s in the age uncertainty",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    example_data$age_un[1:3] <-
      0

    expect_no_error(
      res <-
        prepare_data(
          example_data,
          working_units = "levels",
          rand = 1
        )
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
  "prepare_data with working_units='bins' and rand=1 works if there are 0s in the age uncertainty",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    example_data$age_un[1:3] <-
      0

    expect_no_error(
      res <-
        prepare_data(
          example_data,
          working_units = "bins",
          bin_size = 500,
          rand = 1
        )
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
  "prepare_data with working_units='MW' and rand=1 works if there are 0s in the age uncertainty",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    example_data$age_un[1:3] <-
      0

    expect_no_error(
      res <-
        prepare_data(
          example_data,
          working_units = "MW",
          bin_size = 500,
          number_of_shifts = 5,
          rand = 1
        )
    )

    expect_false(
      identical(
        res[[1]][[1]]$data$age,
        example_data$age$age
      )
    )
  }
)


# NAs in Age_uncertainty
test_that(
  "prepare_data with working_units='levels' works if there are NAs in age_un data",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    example_data$age_un[1:3] <-
      NA

    expect_no_error(
      prepare_data(
        example_data,
        working_units = "levels",
        rand = 1
      )
    )
  }
)

test_that(
  "prepare_data with working_units='bins' works if there are NAs in age_un data",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    example_data$age_un[1:3] <-
      NA

    expect_no_error(
      prepare_data(
        example_data,
        working_units = "bins",
        bin_size = 500,
        rand = 1
      )
    )
  }
)


test_that(
  "prepare_data with working_units='MW' returns NA in data if there is NA in age_uncertainty",
  {
    example_data <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    example_data$age_un[1:3] <-
      NA

    expect_no_error(
      res <-
        prepare_data(
          example_data,
          working_units = "MW",
          bin_size = 500,
          number_of_shifts = 5,
          rand = 1
        )
    )

    # returns NA in output$age
    expect_true(any(is.na(res[[1]][[1]]$data$age)))
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


# test fails because function silently uses the default parameters
test_that(
  "prepare_data rejects no working_unit as input",
  {
    data_source_bins <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expect_error(
      prepare_data(
        data_source_bins
      )
      # no error programmed into function yet
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



### Bin Size ###

# ---------------------------------------------------------- #
#                   Bin_size Input Tests                     #
# ---------------------------------------------------------- #

# 1.1 "levels"

test_that(
  "prepare_data validates bin_size is not NULL with working_unit='levels'",
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
        bin_size = NULL,
        number_of_shifts = NULL,
        rand = NULL
      ),
      "'bin_size' must be one of the following: 'numeric'"
    )
  }
)

test_that(
  "prepare_data validates bin_size is not Inf with working_unit='levels'",
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
        bin_size = Inf,
        number_of_shifts = NULL,
        rand = NULL
      ),
      "'to' must be a finite number"
    )
  }
)


test_that(
  "prepare_data validates bin_size is not negative with working_unit='levels'",
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
        bin_size = -500,
        number_of_shifts = NULL,
        rand = NULL
      ),
      "wrong sign in 'by' argument"
    )
  }
)

test_that(
  "prepare_data validates bin_size is not NA with working_unit='levels'",
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
        bin_size = NA,
        number_of_shifts = NULL,
        rand = NULL
      ),
      "'bin_size' must be one of the following: 'numeric'"
    )
  }
)

test_that(
  "prepare_data validates bin_size is integer with working_unit='levels'",
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
        bin_size = 500.5,
        number_of_shifts = NULL,
        rand = NULL
      ),
      "' bin_size ' must be a an integer"
    )
  }
)

test_that(
  "prepare_data validates bin_size is not 0 with working_unit='levels'",
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
        bin_size = 0,
        number_of_shifts = NULL,
        rand = NULL
      ),
      "invalid"
    )
  }
)


test_that(
  "prepare_data works with minimum bin_size (1) with working_unit='levels'",
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
        working_units = "levels",
        bin_size = 1,
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
  "prepare_data produces _randomized_ samples with every trial",
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
  "prepare_data samples reproducibly if age_un is supplied, rand > 1 and no seed set manually",
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

# ---------------------------------------------------------- #
#                Additional Edge Case Tests                  #
# ---------------------------------------------------------- #

test_that(
  "prepare_data  throws error with empty data_source_prep list",
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
    x = 1:5),
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

#fails
test_that(
  "prepare_data  throws error with mismatched rownames between community and age",
  {
    example_data <-
      extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      verbose = FALSE
    )
    
    # Change rownames to create mismatch
    rownames(example_data$community)[1:5] <-
      paste0("wrong_", 1:5)
    
    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "levels"
      ),
      # none programmed into the function yet
    )
  }
)

test_that(
  "prepare_data  throws error with age_uncertainty with wrong dimensions",
  {
    example_data <-
      extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      verbose = FALSE
    )
    
    # Add age_uncertainty with wrong number of columns
    example_data$age_un <-
      matrix(1:100, nrow = 10, ncol = 10)
    
    expect_error(
      prepare_data(
        data_source_prep = example_data,
        working_units = "levels",
        rand = 5
      ),
      "`age` must be size 63 or 1, not 10."
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
  "prepare_data handles extremely small bin_size efficiently",
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

# ---------------------------------------------------------- #
#                Age Uncertainty Edge Cases                  #
# ---------------------------------------------------------- #

test_that(
  "prepare_data handles age_uncertainty with all identical values",
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
    result[-1], function(
    x) 
          identical(x[[1]]$data$age, result[[1]][[1]]$data$age)
        )
      )
    )
  }
)




# WIP:
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
      prepare_data(example_data, working_units = "levels")
    result_bins <-
      prepare_data(example_data, working_units = "bins", bin_size = 500)
    result_mw <-
      prepare_data(example_data, working_units = "MW", bin_size = 500, number_of_shifts = 3)
    
    # All should have same top-level structure
    expect_identical(names(result_levels), names(result_bins))
    expect_identical(names(result_levels), names(result_mw))
    
    # All should have same internal structure
    expect_identical(names(result_levels[[1]][[1]]), names(result_bins[[1]][[1]]))
    expect_identical(names(result_levels[[1]][[1]]), names(result_mw[[1]][[1]]))
    
    # Bins should have same column structure
    expect_identical(names(result_levels[[1]][[1]]$bins), names(result_bins[[1]][[1]]$bins))
    expect_identical(names(result_levels[[1]][[1]]$bins), names(result_mw[[1]][[1]]$bins))
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
    expect_false(all(res[[1]][[1]]$data == 0))
    expect_false(length(res[[1]][[1]]$data) == 0)
    # $bins is not all 0 or empty
    expect_false(all(res[[1]][[1]]$bins == 0))
    expect_false(nrow(res[[1]][[1]]$bins) == 0)
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


# ---------------------------------------------------------- #
#                Data Structure Validation Tests             #
# ---------------------------------------------------------- #


test_that(
  "prepare_data handles single-row community data",
  {
    example_data <-
      extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]][1, , drop = FALSE],
      data_age_extract = RRatepol::example_data$sample_age[[1]][1, , drop = FALSE],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]][,1 , drop = FALSE],
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
    result[[1]][[1]]$data) - 1, 1) # -1 for age column
  }
)






















# WIP



# --------------------------------------------------------- #
# Tests for the full functionality within estimate_roc()
# --------------------------------------------------------- #

# working_units = levels

# # 1. default with correct data
# test_that(
#   "prepare_data works correctly within estimate_roc() workflow", {
#   community <-
#     RRatepol::example_data$pollen_data[[1]]
#   age <-
#     RRatepol::example_data$sample_age[[1]]
#   age_uncertainty <-
#     RRatepol::example_data$age_uncertainty[[1]]
# 
#   data_extract <-
#     extract_data(
#       community,
#       age,
#       age_uncertainty
#     )
# 
#   data_smoothed <-
#     smooth_community_data(
#       data_extract,
#       smooth_method = "shep"
#     )
# 
#   data_work <-
#     reduce_data(
#       data_smoothed
#     )
# 
#   # empty working_units uses 'levels'
#   data_prepared <-
#     prepare_data(
#       data_work
#     )
# })

# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age", {
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
})

# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age", {
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
})


# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age", {
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
})

# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age", {
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
})


# NAs in community
test_that(
  "prepare_data works correctly within the workflow of estimate_roc() with NAs in community", {
  community <-
    RRatepol::example_data$pollen_data[[1]]

  community[1:10, -1] <-
    NA

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

  expect_identical(
    data_prepared[[1]][[1]]$bins$name,
    rownames(data_work$community)
  )
})

# NAs in community
test_that(
  "prepare_data works correctly within the workflow of estimate_roc() with NAs in community", {
  community <-
    RRatepol::example_data$pollen_data[[1]]

  community[1:10, -1] <-
    NA

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

  expect_identical(
    data_prepared[[1]][[1]]$bins$name,
    rownames(data_work$community)
  )
})

# NAs in community
test_that(
  "prepare_data works correctly within the workflow of estimate_roc() with NAs in community", {
  community <-
    RRatepol::example_data$pollen_data[[1]]

  community[1:10, -1] <-
    NA

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

  expect_identical(
    data_prepared[[1]][[1]]$bins$name,
    rownames(data_work$community)
  )
})

# NAs in community
test_that(
  "prepare_data works correctly within the workflow of estimate_roc() with NAs in community", {
  community <-
    RRatepol::example_data$pollen_data[[1]]

  community[1:10, -1] <-
    NA

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

  # empty working_units uses 'levels'
  data_prepared <-
    prepare_data(
      data_work
    )

  expect_identical(
    data_prepared[[1]][[1]]$bins$name,
    rownames(data_work$community)
  )
})



# working_units = bins

# 1. default with correct data
# test_that(
#   "prepare_data works correctly within estimate_roc() workflow", {
#   community <-
#     RRatepol::example_data$pollen_data[[1]]
#   age <-
#     RRatepol::example_data$sample_age[[1]]
#   age_uncertainty <-
#     RRatepol::example_data$age_uncertainty[[1]]
# 
#   data_extract <-
#     extract_data(
#       community,
#       age,
#       age_uncertainty
#     )
# 
#   data_smoothed <-
#     smooth_community_data(
#       data_extract,
#       smooth_method = "shep"
#     )
# 
#   data_work <-
#     reduce_data(
#       data_smoothed
#     )
# 
#   data_prepared <-
#     prepare_data(
#       data_work,
#       working_units = "bins"
#     )
# })

# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age", {
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
})

# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age", {
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
})


# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age", {
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
})

# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age", {
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
})


# NAs in community
# test_that(
#   "prepare_data works correctly within the workflow of estimate_roc() with NAs in community", {
#   community <-
#     RRatepol::example_data$pollen_data[[1]]
# 
#   community[1:10, -1] <-
#     NA
# 
#   age <-
#     RRatepol::example_data$sample_age[[1]]
# 
#   age_uncertainty <-
#     RRatepol::example_data$age_uncertainty[[1]]
# 
#   data_extract <-
#     extract_data(
#       community,
#       age,
#       age_uncertainty
#     )
# 
#   data_smoothed <-
#     smooth_community_data(
#       data_extract,
#       smooth_method = "shep"
#     )
# 
#   data_work <-
#     reduce_data(
#       data_smoothed
#     )
# 
#   data_prepared <-
#     prepare_data(
#       data_work,
#       working_units = "bins"
#     )
# 
#   # expect_identical(
#   #   data_prepared[[1]][[1]]$bins$name,
#   #   rownames(data_work$community)
#   # )
# })
# 
# NAs in community
# test_that(
#   "prepare_data works correctly within the workflow of estimate_roc() with NAs in community", {
#   community <-
#     RRatepol::example_data$pollen_data[[1]]
# 
#   community[1:10, -1] <-
#     NA
# 
#   age <-
#     RRatepol::example_data$sample_age[[1]]
# 
#   age_uncertainty <-
#     RRatepol::example_data$age_uncertainty[[1]]
# 
#   data_extract <-
#     extract_data(
#       community,
#       age,
#       age_uncertainty
#     )
# 
#   data_smoothed <-
#     smooth_community_data(
#       data_extract,
#       smooth_method = "m.avg"
#     )
# 
#   data_work <-
#     reduce_data(
#       data_smoothed
#     )
# 
#   data_prepared <-
#     prepare_data(
#       data_work,
#       working_units = "bins"
#     )
# 
#   # expect_identical(
#   #   data_prepared[[1]][[1]]$bins$name,
#   #   rownames(data_work$community)
#   # )
# })
# 
# NAs in community
# test_that(
#   "prepare_data works correctly within the workflow of estimate_roc() with NAs in community", {
#   community <-
#     RRatepol::example_data$pollen_data[[1]]
# 
#   community[1:10, -1] <-
#     NA
# 
#   age <-
#     RRatepol::example_data$sample_age[[1]]
# 
#   age_uncertainty <-
#     RRatepol::example_data$age_uncertainty[[1]]
# 
#   data_extract <-
#     extract_data(
#       community,
#       age,
#       age_uncertainty
#     )
# 
#   data_smoothed <-
#     smooth_community_data(
#       data_extract,
#       smooth_method = "grim"
#     )
# 
#   data_work <-
#     reduce_data(
#       data_smoothed
#     )
# 
#   data_prepared <-
#     prepare_data(
#       data_work,
#       working_units = "bins"
#     )
# 
#   # expect_identical(
#   #   data_prepared[[1]][[1]]$bins$name,
#   #   rownames(data_work$community)
#   # )
# })
# 
# NAs in community
# test_that(
#   "prepare_data works correctly within the workflow of estimate_roc() with NAs in community", {
#   community <-
#     RRatepol::example_data$pollen_data[[1]]
# 
#   community[1:10, -1] <-
#     NA
# 
#   age <-
#     RRatepol::example_data$sample_age[[1]]
# 
#   age_uncertainty <-
#     RRatepol::example_data$age_uncertainty[[1]]
# 
#   data_extract <-
#     extract_data(
#       community,
#       age,
#       age_uncertainty
#     )
# 
#   data_smoothed <-
#     smooth_community_data(
#       data_extract,
#       smooth_method = "age.w"
#     )
# 
#   data_work <-
#     reduce_data(
#       data_smoothed
#     )
# 
#   data_prepared <-
#     prepare_data(
#       data_work,
#       working_units = "bins"
#     )
# 
#   # expect_identical(
#   #   data_prepared[[1]][[1]]$bins$name,
#   #   rownames(data_work$community)
#   # )
# })
# 
# 
# 
# working_units = MW

# 1. default with correct data
# test_that(
#   "prepare_data works correctly within estimate_roc() workflow", {
#   community <-
#     RRatepol::example_data$pollen_data[[1]]
#   age <-
#     RRatepol::example_data$sample_age[[1]]
#   age_uncertainty <-
#     RRatepol::example_data$age_uncertainty[[1]]
# 
#   data_extract <-
#     extract_data(
#       community,
#       age,
#       age_uncertainty
#     )
# 
#   data_smoothed <-
#     smooth_community_data(
#       data_extract,
#       smooth_method = "shep"
#     )
# 
#   data_work <-
#     reduce_data(
#       data_smoothed
#     )
# 
#   data_prepared <-
#     prepare_data(
#       data_work,
#       working_units = "MW"
#     )
# })
# 
# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age", {
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
})

# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age", {
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
})


# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age", {
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
})

# NAs in age (known issue)
test_that(
  "prepare_data works correctly within estimate_roc() workflow with NAs in age", {
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
})


# NAs in community
# test_that(
#   "prepare_data works correctly within the workflow of estimate_roc() with NAs in community", {
#   community <-
#     RRatepol::example_data$pollen_data[[1]]
# 
#   community[1:10, -1] <-
#     NA
# 
#   age <-
#     RRatepol::example_data$sample_age[[1]]
# 
#   age_uncertainty <-
#     RRatepol::example_data$age_uncertainty[[1]]
# 
#   data_extract <-
#     extract_data(
#       community,
#       age,
#       age_uncertainty
#     )
# 
#   data_smoothed <-
#     smooth_community_data(
#       data_extract,
#       smooth_method = "shep"
#     )
# 
#   data_work <-
#     reduce_data(
#       data_smoothed
#     )
# 
#   data_prepared <-
#     prepare_data(
#       data_work,
#       working_units = "MW"
#     )
# 
#   # expect_identical(
#   #   data_prepared[[1]][[1]]$bins$name,
#   #   rownames(data_work$community)
#   # )
# })
# 
# NAs in community
# test_that(
#   "prepare_data works correctly within the workflow of estimate_roc() with NAs in community", {
#   community <-
#     RRatepol::example_data$pollen_data[[1]]
# 
#   community[1:10, -1] <-
#     NA
# 
#   age <-
#     RRatepol::example_data$sample_age[[1]]
# 
#   age_uncertainty <-
#     RRatepol::example_data$age_uncertainty[[1]]
# 
#   data_extract <-
#     extract_data(
#       community,
#       age,
#       age_uncertainty
#     )
# 
#   data_smoothed <-
#     smooth_community_data(
#       data_extract,
#       smooth_method = "m.avg"
#     )
# 
#   data_work <-
#     reduce_data(
#       data_smoothed
#     )
# 
#   data_prepared <-
#     prepare_data(
#       data_work,
#       working_units = "MW"
#     )
# 
#   # expect_identical(
#   #   data_prepared[[1]][[1]]$bins$name,
#   #   rownames(data_work$community)
#   # )
# })
# 
# NAs in community
# test_that(
#   "prepare_data works correctly within the workflow of estimate_roc() with NAs in community", {
#   community <-
#     RRatepol::example_data$pollen_data[[1]]
# 
#   community[1:10, -1] <-
#     NA
# 
#   age <-
#     RRatepol::example_data$sample_age[[1]]
# 
#   age_uncertainty <-
#     RRatepol::example_data$age_uncertainty[[1]]
# 
#   data_extract <-
#     extract_data(
#       community,
#       age,
#       age_uncertainty
#     )
# 
#   data_smoothed <-
#     smooth_community_data(
#       data_extract,
#       smooth_method = "grim"
#     )
# 
#   data_work <-
#     reduce_data(
#       data_smoothed
#     )
# 
#   data_prepared <-
#     prepare_data(
#       data_work,
#       working_units = "MW"
#     )
# 
#   # expect_identical(
#   #   data_prepared[[1]][[1]]$bins$name,
#   #   rownames(data_work$community)
#   # )
# })
# 
# NAs in community
# test_that(
#   "prepare_data works correctly within the workflow of estimate_roc() with NAs in community", {
#   community <-
#     RRatepol::example_data$pollen_data[[1]]
# 
#   community[1:10, -1] <-
#     NA
# 
#   age <-
#     RRatepol::example_data$sample_age[[1]]
# 
#   age_uncertainty <-
#     RRatepol::example_data$age_uncertainty[[1]]
# 
#   data_extract <-
#     extract_data(
#       community,
#       age,
#       age_uncertainty
#     )
# 
#   data_smoothed <-
#     smooth_community_data(
#       data_extract,
#       smooth_method = "age.w"
#     )
# 
#   data_work <-
#     reduce_data(
#       data_smoothed
#     )
# 
#   data_prepared <-
#     prepare_data(
#       data_work,
#       working_units = "MW"
#     )
# 
#   # expect_identical(
#   #   data_prepared[[1]][[1]]$bins$name,
#   #   rownames(data_work$community)
#   # )
# })
# 
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



# tests for specific parameter combinations:

# Invalid parameter combinations:

test_that(
  "prepare_data throws message if invalid combinations of parameters are silently overridden", {
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
    # "Invalid parameter 'number_of_shifts = 5' to working_units = 'levels'. number_of_shifts is automatically set to = 1"
    )
})


test_that(
  "prepare_data sets number_of_shifts =1 with user-supplied: 'MW' and number_of_shifts = 0", {
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
})

# Test combinations of rand and age_un.
# Q: is it the same whether if the user sets rand = 1 or if the function does it if rand = NULL?
# expect identical results for rand = 1 and rand = NULL.

test_that(
  "Prepare data returns identical results for rand=NULL and rand=1", {
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
})


test_that(
  "Prepare data returns identical results for rand=NULL and rand=1 (without age_un)", {
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
})


# Is rand ignored if no age_un?
test_that(
  "prepare_data ignores rand parameter if no age_un is supplied", {
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
})


test_that(
  "prepare_data returns identical results across randomizations if age_un is not supplied but rand is", {
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
    x) identical(
    x, data_prepared_rand_100[[1]]),
        logical(1)
      )
      )
    )
})


# Q: What happens if age_un has only 1 row but rand = 100?

test_that(
  "prepare_data returns randomized results if age_un has only 2 rows but rand = 100", {
  community <-
    RRatepol::example_data$pollen_data[[1]]

  age <-
    RRatepol::example_data$sample_age[[1]]

  age_un <-
    RRatepol::example_data$age_uncertainty[[1]][1:2, 1:63]


  data_extract <-
    extract_data(
      community,
      age,
      age_un
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

  data_prepared_rand_100_2rows <-
    prepare_data(
      data_work,
      working_units = "MW",
      number_of_shifts = 1,
      rand = 100
    )

  expect_false(
    identical(
      data_prepared_rand_100_2rows[[100]][[1]]$bins,
      data_prepared_rand_100_2rows[[7]][[1]]$bins
    )
  )


  # Control:
  age_un <-
    RRatepol::example_data$age_uncertainty[[1]]


  data_extract <-
    extract_data(
      community,
      age,
      age_un
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

  data_prepared_rand_100_full <-
    prepare_data(
      data_work,
      working_units = "MW",
      number_of_shifts = 1,
      rand = 100
    )
  

  expect_identical(
    data_prepared_rand_100_full[[1]],
    data_prepared_rand_100_2rows[[1]]
  )

})







### WIP:
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
    data_prepared_1row),
      100
    )

    expect_true(
      all(
        vapply(
          data_prepared_1row[-1],
          function(
    x) identical(
    x[[1]]$data$age, data_prepared_1row[[1]][[1]]$data$age),
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
          data_prepared_2rows[-1],
          function(
    x) identical(
    x[[1]]$data$age, data_prepared_2rows[[1]][[1]]$data$age),
          logical(1)
        )
      )
    )

    # Test 3: Compare specific samples to ensure they can differ
    expect_false(
      identical(
        data_prepared_2rows[[100]][[1]]$data$age,
        data_prepared_2rows[[7]][[1]]$data$age
      )
    )

    # Test 4: Verify sampling distribution with 2 rows
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

      # Test 5: Control - compare against full uncertainty data
  age_un_full <-
    RRatepol::example_data$age_uncertainty[[1]]
  
  data_extract_full <-
    extract_data(community, age, age_un_full)
  data_smoothed_full <-
    smooth_community_data(data_extract_full, smooth_method = "age.w")
  data_work_full <-
    reduce_data(data_smoothed_full)
  
  data_prepared_full <-
    prepare_data(
    data_work_full,
    working_units = "MW",
    number_of_shifts = 1,
    rand = 100
  )
  
  # Full uncertainty should have more variation than limited rows
  age_vectors_full <-
    lapply(data_prepared_full, function(x) x[[1]]$data$age)
  unique_patterns_full <-
    length(unique(age_vectors_full))
  
  expect_true(unique_patterns_full > unique_age_patterns)

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
      extract_data(community, age, age_un_full)
    data_smoothed_full <-
      smooth_community_data(data_extract_full, smooth_method = "age.w")
    data_work_full <-
      reduce_data(data_smoothed_full)


    # Test 6: Verify first sample is identical regardless of uncertainty data size
    # (when using same random seed, first sample should use same row)
    set.seed(123)
    data_prep_full_seeded1 <-
      prepare_data(data_work_full, working_units = "MW", number_of_shifts = 1, rand = 1)
    set.seed(123)
    data_prep_full_seeded2 <-
      prepare_data(data_work_full, working_units = "MW", number_of_shifts = 1, rand = 1)

    # First random sample should be from same row index
    expect_identical(
      data_prep_full_seeded1[[1]][[1]]$data$age,
      data_prep_full_seeded2[[1]][[1]]$data$age
    )
  }
)



# Q: What if age_un is empty or has NAs?

