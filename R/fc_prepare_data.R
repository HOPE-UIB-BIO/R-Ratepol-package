#' @title Prepare data for single RoC estimation
#'
#' @inheritParams fc_estimate_RoC
#' @param data_source_prep
#' List with `community` and `age`
#' @description
#' Create a single list with all information needed to estimate RoC.
#' This is done because such list can be then evaluated in parallel.
#' @keywords internal
fc_prepare_data <-
    function(data_source_prep,
             Working_Units = c("levels", "bins", "MW"),
             bin_size = 500,
             Number_of_shifts = 5,
             rand = NULL) {
        RUtilpol::check_class("data_source_prep", "list")

        RUtilpol::check_class("Working_Units", "character")

        RUtilpol::check_vector_values("Working_Units", c("levels", "bins", "MW"))

        Working_Units <- match.arg(Working_Units)

        RUtilpol::check_class("bin_size", "numeric")

        RUtilpol::check_if_integer("bin_size")

        RUtilpol::check_class("rand", c("NULL", "numeric"))

        # check the condition
        is_shift_present <-
            Working_Units == "MW" && (Number_of_shifts != 0)

        if (
            is_shift_present == FALSE
        ) {
            Number_of_shifts <- 1
        } else {
            RUtilpol::check_class("Number_of_shifts", "numeric")
            RUtilpol::check_if_integer("Number_of_shifts")
        }

        is_rand_present <-
            (is.null(rand) == FALSE)

        if (
            is_rand_present == TRUE
        ) {
            RUtilpol::check_if_integer("rand")
        } else {
            rand <- 1
        }

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
            is_uncertit_present == TRUE
        ) {
            if (
                is_rand_present == TRUE
            ) {
                random_value <-
                    sample(
                        c(
                            1:max(1, nrow(data_source_prep$age_un))
                        ),
                        rand
                    )

                random_age <-
                    purrr::map(
                        .x = seq_along(random_value),
                        .f = ~ as.numeric(
                            data_source_prep$age_un[random_value[.x], ]
                        )
                    )
            } else {
                random_age <-
                    list(data_source_prep$age$age)
            }
        }

        data_merge <-
            data_source_prep$community %>%
            tibble::rownames_to_column("row_name") %>%
            dplyr::left_join(
                data_source_prep$age %>%
                    tibble::rownames_to_column("row_name"),
                by = "row_name"
            ) %>%
            dplyr::relocate(.data$age) %>%
            tibble::column_to_rownames("row_name")

        shift_vec <-
            c(1:Number_of_shifts) %>%
            rlang::set_names(
                nm = as.character(1:Number_of_shifts)
            )

        rand_vec <-
            c(1:rand) %>%
            rlang::set_names(
                nm = as.character(1:rand)
            )

        if (
            # is_rand_present is TRUE
            is_uncertit_present == FALSE
        ) {
            rand_vec %>%
                purrr::map(
                    .f = ~ purrr::map(
                        .x = shift_vec,
                        .f = ~
                            list(
                                data = data_merge,
                                bins = bin_dummy %>%
                                    dplyr::filter(shift == .x)
                            )
                    )
                ) %>%
                return()
        } else {
            # is_uncertit_present is TRUE
            rand_vec %>%
                purrr::map(
                    .f = ~ {
                        data_with_random_age <-
                            data_merge %>%
                            dplyr::mutate(
                                age = random_age[[.x]]
                            )

                        purrr::map(
                            .x = shift_vec,
                            .f = ~
                                list(
                                    data = data_with_random_age,
                                    bins = bin_dummy %>%
                                        dplyr::filter(shift == .x)
                                ) %>%
                                    return()
                        ) %>%
                            return(res)
                    }
                ) %>%
                return()
        }
        # stop("Not a valid setting for creating bins")
    }
