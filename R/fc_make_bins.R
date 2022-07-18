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

        if (
            Working_Units == "levels"
        ) {
            n_levels <-
                nrow(data_source_prep$age)

            n_res <- n_levels - 1

            res <-
                data.frame(
                    name = rownames(data_source_prep$age)[1:n_res]
                ) %>%
                dplyr::mutate(
                    start = as.character(name),
                    end = as.character(
                        rownames(data_source_prep$age)[2:n_levels]
                    ),
                    age = data_source_prep$age$age[1:n_res],
                    label = paste(start, end, sep = "-")
                )
            return(res)
        }

        util_check_class("bin_size", "numeric")

        util_check_if_integer("bin_size")

        bin_oldest <-
            ceiling(
                max(data_source_bins$age$age)
            )

        bin_youngest <-
            floor(
                min(data_source_bins$age$age)
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
                start = name,
                end = start + bin_size,
                age = (end + start) / 2,
                label = paste(start, end, sep = "-")
            )

        return(res)
    }
