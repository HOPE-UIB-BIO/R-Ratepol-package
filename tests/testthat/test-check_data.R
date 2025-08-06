test_that(
  "function throws error if wrong input is supplied",
  {
    data_source_check <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    class(data_source_check) <-
      "data.frame"

    expect_error(
      check_data(
        data_source_check,
        "'data_source_check' must be one of the following: 'list'"
      )
    )
  }
)


test_that(
  "check_data() tells the user when extract_data() returns empty data",
  {
    # Create deliberately empty but valid-structured inputs
    empty_community <-
      data.frame(sample_id = character(0))
    empty_age <-
      data.frame(sample_id = character(0), age = numeric(0))
    empty_uncertainty <-
      matrix(nrow = 0, ncol = 0)

    # Get result from extract_data() (assuming it doesn’t crash)
    result_empty <-
      extract_data(
        data_community_extract = empty_community,
        data_age_extract = empty_age,
        age_uncertainty = empty_uncertainty,
        verbose = FALSE
      )

    msg <-
      suppressWarnings(
        capture.output(
          check_data(
            result_empty
          ),
          type = "message"
        )
      )


    # Now check that check_data() reproduces the empty input data
    expect_true(
      any(
        grepl(
          "Community data have 0 taxa and 0 samples.  Age data have 0 samples",
          msg
        )
      )
    )

    expect_true(
      any(
        grepl(
          "Age data has values of min Inf, max -Inf, mean NaN, and median NA",
          msg
        )
      )
    )
  }
)



test_that(
  "does not return NA in min, max, mean, median age if there is NA in the data",
  {
    # introduce NAs:
    example_community_NA <-
      RRatepol::example_data$pollen_data[[1]]

    example_community_NA[5:10, 3] <-
      NA

    example_age_NA <-
      RRatepol::example_data$sample_age[[1]]

    example_age_NA[5:10, 3] <-
      NA

    example_uncertainty_NA <-
      RRatepol::example_data$age_uncertainty[[1]]

    example_uncertainty_NA[5:10, 5] <-
      NA

    result_NA <-
      extract_data(
        example_community_NA,
        example_age_NA,
        example_uncertainty_NA
      )

    msg <-
      capture.output(
        check_data(
          result_NA
        ),
        type = "message"
      )

    expect_false(
      any(
        grepl(
          "Age data has values of min NA, max NA, mean NA, and median NA",
          msg
        )
      )
    )
  }
)
