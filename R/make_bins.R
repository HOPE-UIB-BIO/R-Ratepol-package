#' @title Create template for RoC estimations
#'
#' @inheritParams estimate_roc
#' @param data_source_bins
#' List with `community` and `age`
#' @keywords internal
make_bins <-
    function(data_source_bins,
             working_units = c("levels", "bins", "MW"),
             bin_size = 500,
             number_of_shifts = 5) {
        RUtilpol::check_class("data_source_bins", "list")

        RUtilpol::check_class("working_units", "character")

        RUtilpol::check_vector_values(
            "working_units",
            c("levels", "bins", "MW")
        )

        working_units <- match.arg(working_units)

        age_dat <-
            data_source_bins$age

        if (
            working_units == "levels"
        ) {
            n_levels <-
                nrow(age_dat)

            age_dat_longer <-
                dplyr::bind_rows(
                    age_dat,
                    data.frame(
                        age = Inf,
                        row.names = c("Inf")
                    )
                )

            age_vec <-
                age_dat_longer[, "age"]

            res <-
                data.frame(
                    name = rownames(age_dat_longer)[1:n_levels],
                    shift = 1
                ) %>%
                dplyr::mutate(
                    age_diff = abs(
                        age_vec[1 + (1:n_levels)] -
                            age_vec[1:n_levels]
                    ),
                    age_diff = ifelse(.data$age_diff == 0, 0.1, .data$age_diff),
                    start = as.character(.data$name),
                    end = as.character(
                        rownames(age_dat_longer)[1 + 1:n_levels]
                    ),
                    res_age = age_vec[1:n_levels],
                    label = paste(.data$start, .data$end, sep = "-")
                )
            return(res)
        }

        RUtilpol::check_class("bin_size", "numeric")

        RUtilpol::check_if_integer("bin_size")

        bin_oldest <-
            ceiling(
                max(age_dat$age)
            ) + bin_size

        bin_youngest <-
            floor(
                min(age_dat$age)
            )

        bin_breaks <- seq(
            from = bin_youngest,
            to = bin_oldest,
            by = bin_size
        )

        if (
            working_units == "bins"
        ) {
            res_df <-
                data.frame(
                    name = bin_breaks,
                    shift = 1
                )
        } else {
            RUtilpol::check_class("number_of_shifts", "numeric")

            RUtilpol::check_if_integer("number_of_shifts")
        }

        if (
            working_units == "MW"
        ) {
            shift_value <- bin_size / number_of_shifts

            shift_increment <-
                shift_value * sort(
                    rep(0:(number_of_shifts - 1), length(bin_breaks))
                )

            res_df <-
                data.frame(
                    name = (rep(bin_breaks, number_of_shifts) +
                        shift_increment),
                    shift = sort(rep(c(1:number_of_shifts), length(bin_breaks)))
                )
        }

        res <-
            res_df %>%
            dplyr::mutate(
                age_diff = bin_size,
                start = .data$name,
                end = .data$start + bin_size,
                res_age = (.data$end + .data$start) / 2,
                label = paste(.data$start, .data$end, sep = "-")
            )

        return(res)
    }
