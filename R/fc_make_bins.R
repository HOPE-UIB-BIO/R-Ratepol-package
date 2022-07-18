fc_make_bins <-
    function(data_source_bins,
             Working_Units = c("levels", "bins", "MW"),
             bin_size = 500,
             Number_of_shifts = 5) {
        util_check_class("data_source_bins", "list")

        util_check_class("Working_Units", "character")

        util_check_vector_values(
            "Working_Units",
            c("levels", "bins", "MW")
        )

        Working_Units <- match.arg(Working_Units)

        age_dat <-
            data_source_bins$age

        if (
            Working_Units == "levels"
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
                    name = rownames(age_dat_longer)[1:n_levels]
                ) %>%
                dplyr::mutate(
                    age_diff = abs(
                        age_vec[1 + (1:n_levels)] -
                            age_vec[1:n_levels]
                    ),
                    start = as.character(name),
                    end = as.character(
                        rownames(age_dat_longer)[1 + 1:n_levels]),
                    res_age = age_vec[1:n_levels],
                    label = paste(start, end, sep = "-")
                )
            return(res)
        }

        util_check_class("bin_size", "numeric")

        util_check_if_integer("bin_size")

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
            Working_Units == "bins"
        ) {
            res_df <-
                data.frame(
                    name = bin_breaks
                )
        } else {
            util_check_class("Number_of_shifts", "numeric")

            util_check_if_integer("Number_of_shifts")
        }

        if (
            Working_Units == "MW"
        ) {
            shift_value <- bin_size / Number_of_shifts

            shift_increment <-
                shift_value * sort(
                    rep(0:(Number_of_shifts - 1), length(bin_breaks))
                )

            res_df <-
                data.frame(
                    name = (rep(bin_breaks, Number_of_shifts) +
                        shift_increment),
                    shift = sort(rep(c(1:Number_of_shifts), length(bin_breaks)))
                )
        }

        res <-
            res_df %>%
            dplyr::mutate(
                age_diff = bin_size,
                start = name,
                end = start + bin_size,
                res_age = (end + start) / 2,
                label = paste(start, end, sep = "-")
            )

        return(res)
    }
