fc_prepare_data_bins <-
    function(data_source_prep,
             bin_size = 500,
             Number_of_shifts = 5,
             rand = NULL) {
        util_check_class("data_source_prep", "list")

        util_check_class("bin_size", "numeric")

        util_check_if_integer("bin_size")

        util_check_class("rand", c("NULL", "numeric"))

        # create age bins
       if (
        Number_of_shifts == 0
       ) {
         bins <-
            fc_make_bins(
                data_source_bin = data_source_prep,
                bin_size = bin_size
            )
       } else {
        bins <-
            fc_make_mw(
                data_source_mw = data_source_prep,
                bin_size = bin_size,
                Number_of_shifts = Number_of_shifts
            )
       }

        if (
            is.null(rand) == FALSE
        ) {
            util_check_if_integer("rand")

               if (
                is.null(data_source_prep$age_un) == TRUE
            ) {
                res <-
                    rep(
                        list(
                            data = data_source_prep,
                            bins = bins
                        ),
                        rand
                    ) %>%
                    rlang::set_names(
                        nm = 1:rand
                    )
            } else {
                res <-
                    c(1:rand) %>%
                    rlang::set_names(
                        nm = 1:rand
                    ) %>%
                    purrr::map(
                        .f = ~ {
                            data_temp <-
                                data_source_prep

                            random_value <-
                                sample(
                                    c(
                                        1:max(1, nrow(data_temp$age_un))
                                    ),
                                    1
                                )

                            data_temp$age$age <-
                                as.numeric(data_temp$age_un[random_value, ])

                            list(
                                data = data_temp,
                                bins = bins
                            ) %>%
                                return()
                        }
                    )
            }
        } else {
            res <-
                list(
                    list(
                        data = data_source_prep,
                        bins = bins
                    )
                ) %>%
                rlang::set_names(
                    nm = 1
                )
        }
        return(res)
    }
