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
