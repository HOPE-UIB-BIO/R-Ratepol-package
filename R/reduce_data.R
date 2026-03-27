#' @title Reduce datasets
#' @param data_source_reduce `list` with `community`, `age`, and `age_un`
#' @param check_taxa
#' `logical`. Should columns be checked for all-zero taxa?
#' @param check_levels
#' `logical`. Should rows be checked for empty levels?
#' @description
#' Check the community dataset for all-zero taxa and empty levels and
#' filter them out.
#' @return
#' The input `list` with redundant columns and rows removed.
#' @keywords internal
reduce_data <-
  function(data_source_reduce,
           check_taxa = TRUE,
           check_levels = TRUE) {
    util_check_class(data_source_reduce, "list")

    util_check_class(check_taxa, "logical")

    util_check_class(check_levels, "logical")

    assertthat::assert_that(
      !base::is.null(data_source_reduce$community),
      msg = "'data_source_reduce$community' must not be NULL"
    )

    assertthat::assert_that(
      !base::is.null(data_source_reduce$age),
      msg = "'data_source_reduce$age' must not be NULL"
    )

    if (
      isTRUE(check_taxa)
    ) {
      valid_taxa <-
        (colSums(data_source_reduce$community, na.rm = TRUE) > 0)

      data_source_reduce$community <-
        data_source_reduce$community %>%
        dplyr::select(
          dplyr::any_of(
            c(names(valid_taxa[valid_taxa]))
          )
        )
    }

    if (
      isTRUE(check_levels) # if filter out samples without individuals
    ) {
      valid_samples_community <-
        data_source_reduce$community %>%
        tibble::rownames_to_column("sample_id") %>%
        dplyr::mutate(
          row_sum = base::rowSums(
            dplyr::pick(-"sample_id"),
            na.rm = TRUE
          )
        ) %>%
        dplyr::filter(.data$row_sum > 0) %>%
        dplyr::pull("sample_id")

      valid_levels_age_comm <-
        intersect(
          valid_samples_community, rownames(data_source_reduce$age)
        )

      data_source_reduce$community <-
        data_source_reduce$community %>%
        tibble::rownames_to_column("sample_id") %>%
        dplyr::filter(.data$sample_id %in% valid_levels_age_comm) %>%
        tibble::column_to_rownames("sample_id")

      data_source_reduce$age <-
        data_source_reduce$age %>%
        tibble::rownames_to_column("sample_id") %>%
        dplyr::filter(.data$sample_id %in% valid_levels_age_comm) %>%
        tibble::column_to_rownames("sample_id")

      if (
        isFALSE(is.null(data_source_reduce$age_un))
      ) {
        valid_samples_all <-
          intersect(
            valid_levels_age_comm,
            colnames(data_source_reduce$age_un)
          )

        data_source_reduce$community <-
          data_source_reduce$community %>%
          tibble::rownames_to_column("sample_id") %>%
          dplyr::filter(.data$sample_id %in% valid_samples_all) %>%
          tibble::column_to_rownames("sample_id")

        data_source_reduce$age <-
          data_source_reduce$age %>%
          tibble::rownames_to_column("sample_id") %>%
          dplyr::filter(.data$sample_id %in% valid_samples_all) %>%
          tibble::column_to_rownames("sample_id")

        data_source_reduce$age_un <-
          data_source_reduce$age_un[, valid_samples_all, drop = FALSE]
      }
    }

    return(data_source_reduce)
  }
