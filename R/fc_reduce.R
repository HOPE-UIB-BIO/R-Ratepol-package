#' @title Reduce datasets
#' @param data_source List with `community`, `age`, and `age_un`
#' @param check_taxa Logical. Should columns be check for redundnat data?
#' @param check_rows Logical. Should rows be check for redundnat data?
#' @description
#' Check the community dataset for redundnat taxa and levels
#' and filter them out.
fc_reduce <-
  function(data_source,
           check_taxa = TRUE,
           check_levels = TRUE) {
    util_check_class("data_source", "list")

    util_check_class("check_taxa", "logical")

    util_check_class("check_levels", "logical")

    if (
      check_taxa == TRUE
    ) {
      valid_taxa <-
        (colSums(data_source$community, na.rm = TRUE) > 0)

      data_source$community <-
        data_source$community %>%
        dplyr::select(
          dplyr::any_of(
            c(names(valid_taxa[valid_taxa]))
          )
        )
    }

    if (
      check_levels == TRUE # if filter out samples without individuals
    ) {
      valid_levels <-
        (rowSums(data_source$community, na.rm = TRUE) > 0)

      data_source$age <-
        data_source$age %>%
        tibble::rownames_to_column("sample_id") %>%
        dplyr::filter(
          names(valid_levels) %in% .data$sample_id
        ) %>%
        tibble::column_to_rownames("sample_id")

      data_source$community <-
        data_source$community %>%
        tibble::rownames_to_column("sample_id") %>%
        dplyr::filter(
          names(valid_levels) %in% .data$sample_id
        ) %>%
        tibble::column_to_rownames("sample_id")

      if (
        is.null(data_source$age_un) == FALSE
      ) {
        data_source$age_un <-
          data_source$age_un[, valid_levels]
      }
    }

    return(data_source)
  }
