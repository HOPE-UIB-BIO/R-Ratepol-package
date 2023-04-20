#' @title Extract and check the data
#'
#' @param data_community_extract
#' Data.frame. Community data with species as columns and
#' levels (samples) as rows. First column should be `sample_id` (character).
#' @param data_age_extract
#' Data.frame with two columns:
#' \itemize{
#' \item `sample_id` - unique ID of each level (character)
#' \item `age` - age of level (numeric)
#' }
#' @param age_uncertainty
#' Usage of age uncertainty form Age-depth models. Either:
#' \itemize{
#' \item matrix with number of columns as number of samples. Each column is one sample,
#'  each row is one age sequence from age-depth model. Age sequence is randomly
#'  sampled from age-depth model uncertainties at the beginning of each run.
#'  \item `NULL` - Age uncertainties are not available and, therefore, will not be used.
#' }
#' @param verbose DESCRIPTION.
#' Logical. If `TRUE`, function will output messages about internal processes
#' @description
#' Function for general preparation of input data
extract_data <-
  function(data_community_extract,
           data_age_extract,
           age_uncertainty = NULL,
           verbose = FALSE) {

    # 1. Initial tests -----

    # 1.1 Data types -----

    RUtilpol::check_class("verbose", "logical")

    if (
      isTRUE(verbose)
    ) {
      RUtilpol::output_heading(
        paste(
          "Data extraction started",
          Sys.time()
        ),
        size = "h2"
      )
    }

    RUtilpol::check_class("data_community_extract", "data.frame")

    RUtilpol::check_class("data_age_extract", "data.frame")

    RUtilpol::check_class("age_uncertainty", c("NULL", "matrix"))

    # 1.2. Sample id -----
    # community

    if (
      "sample.id" %in% names(data_community_extract)
    ) {
      usethis::ui_oops(
        paste(
          "'sample.id' was detected in 'data_community'",
          "but 'sample_id' is prefered.",
          "Recommend renaming your data"
        )
      )

      data_community_extract <-
        data_community_extract %>%
        dplyr::rename(sample_id = .data$sample.id)
    }

    RUtilpol::check_col_names("data_community_extract", "sample_id")

    assertthat::assert_that(
      "character" %in% class(data_community_extract$sample_id),
      msg = "Variable 'sample_id' in 'data_community' must
    be a 'character'"
    )

    # age
    if (
      "sample.id" %in% names(data_age_extract)
    ) {
      usethis::ui_oops(
        paste(
          "'sample.id' was detected in 'data_age' but 'sample_id' is prefered.",
          "Recomend renaming your data"
        )
      )

      data_age_extract <-
        data_age_extract %>%
        dplyr::rename(sample_id = .data$sample.id)
    }

    RUtilpol::check_col_names("data_age_extract", "sample_id")

    assertthat::assert_that(
      all(data_community_extract$sample_id %in% data_age_extract$sample_id) &&
        all(data_age_extract$sample_id %in% data_community_extract$sample_id),
      msg = "Variable 'sample_id' must have same values in
    'data_age' and 'data_community'"
    )

    # 1.3. Age test -----

    RUtilpol::check_col_names("data_age_extract", "sample_id")

    assertthat::assert_that(
      "numeric" %in% class(data_age_extract$age),
      msg = "Variable 'age' in 'data_source_age' must be a 'numeric'"
    )

    if (
      isFALSE(is.unsorted(data_age_extract$age))
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
      msg = "Object 'data_community' and 'data_age'
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
      dplyr::select("sample_id", "age")

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
      isFALSE(is.null(age_uncertainty))
    ) {
      RUtilpol::check_class("age_uncertainty", "matrix")

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
      RUtilpol::output_warning(
        paste(
          "Missing data has been detected in community data",
          "and automatically replaces with '0'"
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

    if (
      any(is.na(dat_community$age))
    ) {
      RUtilpol::output_warning(
        paste(
          "Missing 'age' values has detected in age data",
          "Such levels has been filtered out"
        )
      )

      dat_age <-
        dat_age %>%
        tidyr::drop_na(.data$age)
    }

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
      reduce_data(
        data_source_reduce = dat_merge,
        check_taxa = TRUE,
        check_levels = TRUE
      )

    if (
      isTRUE(verbose)
    ) {
      check_data(
        data_source_check = dat_merge
      )

      RUtilpol::output_heading(
        paste(
          "Data extraction completed",
          Sys.time()
        ),
        size = "h2"
      )
    }

    return(dat_merge)
  }
