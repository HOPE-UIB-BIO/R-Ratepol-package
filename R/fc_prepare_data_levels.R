fc_prepare_data_levels <-
    function(data_source_prep,
             rand = NULL) {
        util_check_class("data_source_prep", "list")

        util_check_class("rand", c("NULL", "numeric"))

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
                            data = data_source_prep
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
                                data = data_temp
                            ) %>%
                                return()
                        }
                    )
            }
        } else {
            res <-
                list(
                    data = data_source_prep
                ) %>%
                rlang::set_names(
                    nm = 1
                )
        }
        return(res)
    }
