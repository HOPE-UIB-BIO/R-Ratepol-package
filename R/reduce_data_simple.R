#' @title Reduce datasets in merged format
#' @param data_source_reduce List with `community`, `age`, and `age_un`
#' @param ommit_vars
#' Character vector with names of columns to omit in community data.
#' @inheritParams reduce_data
#' @description
#' Check the community dataset for redundnat taxa and levels
#' and filter them out. This function is simplified due to performance.
#' @keywords internal
reduce_data_simple <-
    function(data_source_reduce,
             ommit_vars = c("label", "res_age", "age_diff"),
             check_taxa = TRUE,
             check_levels = TRUE) {
        util_check_class(data_source_reduce, "data.frame")

        assertthat::assert_that(
            base::nrow(data_source_reduce) > 0,
            msg = "'data_source_reduce' must not be empty"
        )

        util_check_class(check_taxa, "logical")

        assertthat::assert_that(
            base::length(check_taxa) == 1,
            msg = "'check_taxa' must be a single TRUE or FALSE"
        )

        util_check_class(check_levels, "logical")

        assertthat::assert_that(
            base::length(check_levels) == 1,
            msg = "'check_levels' must be a single TRUE or FALSE"
        )

        data_com <-
            subset_community(
                data_source_reduce,
                ommit_vars = ommit_vars
            )

        if (
            isTRUE(check_taxa)
        ) {
            valid_taxa <-
                (colSums(data_com, na.rm = TRUE) > 0)

            data_source_reduce <-
                data_source_reduce %>%
                dplyr::select(
                    dplyr::any_of(
                        c(
                            ommit_vars,
                            names(valid_taxa[valid_taxa])
                        )
                    )
                )
        }

        if (
            isTRUE(check_levels)
        ) {
            valid_levels <-
                (rowSums(data_com, na.rm = TRUE) > 0)

            data_source_reduce <-
                data_source_reduce[valid_levels, ]
        }

        return(data_source_reduce)
    }
