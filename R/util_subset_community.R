util_subset_community <-
    function(data_source,
             ommit_vars = c("label", "res_age", "age_diff", "age")) {
        data_source %>%
            dplyr::select(
                !dplyr::any_of(
                    c(ommit_vars)
                )
            ) %>%
            return()
    }
