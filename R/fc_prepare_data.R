fc_prepare_data <-
    function(data_source_prep,
             Working_Units = c("levels", "bins", "MW"),
             bin_size = 500,
             Number_of_shifts = 5,
             rand = NULL) {
        util_check_class("data_source_prep", "list")

        util_check_class("Working_Units", "character")

        util_check_vector_values("Working_Units", c("levels", "bins", "MW"))

        Working_Units <- match.arg(Working_Units)

        util_check_class("bin_size", "numeric")

        util_check_if_integer("bin_size")

        util_check_class("rand", c("NULL", "numeric"))

        # check the condition
        is_rand_present <-
            (is.null(rand) == FALSE)

        if (
            is_rand_present == TRUE
        ) {
            util_check_if_integer("rand")
        }

        is_shift_present <-
            Working_Units == "MW" && (Number_of_shifts != 0)

        if (
            Working_Units == "levels"
        ) {
            bin_dummy <-
                fc_make_bins(
                    data_source_bins = data_source_prep,
                    Working_Units = "levels"
                )
        } else if (
            is_shift_present == TRUE
        ) {
            bin_dummy <-
                fc_make_bins(
                    data_source_bins = data_source_prep,
                    Working_Units = "MW",
                    bin_size = bin_size,
                    Number_of_shifts = Number_of_shifts
                )
        } else {
            bin_dummy <-
                fc_make_bins(
                    data_source_bins = data_source_prep,
                    Working_Units = "bins",
                    bin_size = bin_size
                )
        }

        is_uncertit_present <-
            (is.null(data_source_prep$age_un) == FALSE)

        if (
            is_uncertit_present == TRUE && is_rand_present == TRUE
        ) {
            random_value <-
                sample(
                    c(
                        1:max(1, nrow(data_source_prep$age_un))
                    ),
                    rand
                )
        }

        data_merge <-
            data_source_prep$community %>%
            tibble::rownames_to_column("row_name") %>%
            dplyr::left_join(
                data_source_prep$age %>%
                    tibble::rownames_to_column("row_name"),
                by = "row_name"
            ) %>%
            dplyr::relocate(age) %>%
            tibble::column_to_rownames("row_name")

        # test individual possibilities
        if (
            is_rand_present == FALSE
        ) {
            list(
                list(
                    data = data_merge,
                    bins = bin_dummy
                )
            ) %>%
                rlang::set_names(
                    nm = 1
                ) %>%
                return()
        } else if (
            # is_rand_present is TRUE
            is_uncertit_present == FALSE
        ) {
            rep(
                list(
                    list(
                        data = data_merge,
                        bins = bin_dummy
                    )
                ),
                rand
            ) %>%
                rlang::set_names(
                    nm = 1:rand
                ) %>%
                return()
        } else {
            # is_uncertit_present is TRUE
            c(1:rand) %>%
                rlang::set_names(
                    nm = 1:rand
                ) %>%
                purrr::map(
                    .f = ~ {
                        random_age <-
                            as.numeric(
                                data_source_prep$age_un[random_value[.x], ]
                            )

                        data_with_rand_age <-
                            data_merge %>%
                                dplyr::mutate(
                                    age = random_age
                                )
                        list(
                            data = data_with_rand_age,
                            bins = bin_dummy
                        ) %>%
                            return()
                    }
                ) %>%
                return()
        }
        # stop("Not a valid setting for creating bins")
    }
