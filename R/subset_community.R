
#' @title Subset data to only contain community (taxa)
#'
#' @param data_source Data.frame wich shuld be cleaned.
#' @param ommit_vars Vector with names of columns to ommit.
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
