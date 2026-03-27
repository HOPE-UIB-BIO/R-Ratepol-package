
#' @title Subset data to only contain community (taxa)
#'
#' @param data_source
#' `data.frame` to clean.
#' @param ommit_vars
#' `character` vector of column names to exclude.
#' @return
#' A `data.frame` containing only taxon columns.
#' @keywords internal
subset_community <-
    function(data_source,
             ommit_vars = c("label", "res_age", "age_diff", "age")) {
        util_check_class(data_source, "data.frame")

        assertthat::assert_that(
            base::nrow(data_source) > 0,
            msg = "'data_source' must not be empty"
        )

        data_source %>%
            dplyr::select(
                !dplyr::any_of(
                    c(ommit_vars)
                )
            ) %>%
            return()
    }
