#' @title Reduce datasets
#' @param data_source_reduce List with `community`, `age`, and `age_un`
#' @param check_taxa Logical. Should columns be check for redundnat data?
#' @param check_levels Logical. Should rows be check for redundnat data?
#' @description
#' Check the community dataset for redundnat taxa and levels
#' and filter them out.
#' @keywords internal
fc_reduce <-
  function(data_source_reduce,
           check_taxa = TRUE,
           check_levels = TRUE) {
    RUtilpol::check_class("data_source_reduce", "list")

    RUtilpol::check_class("check_taxa", "logical")

    RUtilpol::check_class("check_levels", "logical")

    if (
      check_taxa == TRUE
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
      check_levels == TRUE # if filter out samples without individuals
    ) {
      valid_levels <-
        (rowSums(data_source_reduce$community, na.rm = TRUE) > 0)

      data_source_reduce$age <-
        data_source_reduce$age %>%
        tibble::rownames_to_column("sample_id") %>%
        dplyr::filter(
          names(valid_levels) %in% .data$sample_id
        ) %>%
        tibble::column_to_rownames("sample_id")

      data_source_reduce$community <-
        data_source_reduce$community %>%
        tibble::rownames_to_column("sample_id") %>%
        dplyr::filter(
          names(valid_levels) %in% .data$sample_id
        ) %>%
        tibble::column_to_rownames("sample_id")

      if (
        is.null(data_source_reduce$age_un) == FALSE
      ) {
        data_source_reduce$age_un <-
          data_source_reduce$age_un[, valid_levels]
      }
    }

    return(data_source_reduce)
  }
