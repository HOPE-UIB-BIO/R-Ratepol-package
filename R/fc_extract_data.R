fc_extract_data <-
  function(data_community_extract,
           data_age_extract,
           age_uncertainty = NULL,
           verbose = FALSE) {
    util_check_class("verbose", "logical")

    if (verbose == TRUE) {
      util_output_comment(
        paste(
          "Data extraction started",
          Sys.time()
        )
      )
    }

    # 1. Initial tests -----

    # 1.1 Data types -----
    util_check_class("data_community_extract", "data.frame")

    util_check_class("data_age_extract", "data.frame")

    # 1.2. Sample id -----
    # community

    if (
      "sample.id" %in% names(data_community_extract)
    ) {
      usethis::ui_oops(
        paste(
          "'sample.id' was detected but 'sample_id' is prefered.",
          "Recommend renaming your data"
        )
      )

      data_community_extract <-
        data_community_extract %>%
        dplyr::rename(sample_id = .data$sample.id)
    }

    util_check_col_names("data_community_extract", "sample_id")

    assertthat::assert_that(
      "character" %in% class(data_community_extract$sample_id),
      msg = "Variable 'sample_id' in 'data_source_community' must
    be a 'character'"
    )

    # age
    if ("sample.id" %in% names(data_age_extract)) {
      usethis::ui_oops(
        paste(
          "'sample.id' was detected but 'sample_id' is prefered.",
          "Recomend renaming your data"
        )
      )

      data_age_extract <-
        data_age_extract %>%
        dplyr::rename(sample_id = .data$sample.id)
    }

    util_check_col_names("data_age_extract", "sample_id")

    assertthat::assert_that(
      all(data_community_extract$sample_id %in% data_age_extract$sample_id) &&
        all(data_age_extract$sample_id %in% data_community_extract$sample_id),
      msg = "Variable 'sample_id' must have same values in
    'data_source_age' and 'data_source_community'"
    )

    # 1.3. Age test -----

    util_check_col_names("data_age_extract", "sample_id")

    assertthat::assert_that(
      "numeric" %in% class(data_age_extract$age),
      msg = "Variable 'age' in 'data_source_age' must be a 'numeric'"
    )

    if (
      is.unsorted(data_age_extract$age) == FALSE
    ) {
      # order of the age
      data_age_extract <-
        data_age_extract %>%
        dplyr::arrange(.data$age)
    }

    # 1.4. Size test -----

    n_samples_com <- nrow(data_community_extract)
    n_samples_age <- nrow(data_age_extract)

    assertthat::assert_that(
      n_samples_com == n_samples_age,
      msg = "Object 'data_source_community' and 'data_source_age'
    must have the same number of levels"
    )

    # 2. Processes data -----

    # 2.1 Extract data  -----

    # extract both important tables if stored as tibbles
    dat_age <-
      as.data.frame(data_age_extract)
    dat_community <-
      as.data.frame(data_community_extract)

    # select only variable 'age' and 'sample id'
    dat_age <-
      dat_age %>%
      dplyr::select(.data$sample_id, .data$age)

    # 2.2 Row.names  -----

    # add row.names to commity, age, and uncertainty data
    dat_community <-
      dat_community %>%
      tibble::column_to_rownames("sample_id")

    dat_age <-
      dat_age %>%
      tibble::column_to_rownames("sample_id")

    # 2.3 Age uncertainty  -----

    # if age_uncertainty is used
    if (
      is.null(age_uncertainty) == FALSE
    ) {
      util_check_class("age_uncertainty", "matrix")

      n_samples_un <- ncol(age_uncertainty)

      assertthat::assert_that(
        n_samples_age == n_samples_un,
        msg = "Object 'data_source_age' and 'age_uncertainty' must have
      the same number of levels. 'age_uncertainty' must have samples stored as
      columns"
      )

      # save as dataframe
      age_un <- data.frame(age_uncertainty)

      names(age_un) <- dat_age$sample_id
    } else {
      age_un <- NULL
    }

    # 2.4 Missing values ---
    if (
      any(is.na(dat_community))
    ) {
      util_output_comment(
        paste(
          "Missing data has been replaces with '0'",
          "in community data"
        )
      )

      dat_community <-
        dat_community %>%
        dplyr::mutate(
          dplyr::across(
            dplyr::everything(),
            ~ tidyr::replace_na(.x, 0)
          )
        )
    }

    dat_age <-
      dat_age %>%
      tidyr::drop_na(.data$age)

    # 2.5 Summary  ----

    # create list of 3 variables pollen, age, age_un
    dat_merge <-
      list(
        community = dat_community,
        age = dat_age,
        age_un = age_un
      )

    #  exclude redundnat rows and columns
    dat_merge <-
      fc_reduce(
        data_source = dat_merge,
        check_taxa = TRUE,
        check_levels = TRUE
      )

    if (verbose == TRUE) {
      fc_check_data(
        data_source_check = dat_merge
      )

      util_output_comment(
        paste(
          "Data extraction completed",
          Sys.time()
        )
      )
    }

    return(dat_merge)
  }
