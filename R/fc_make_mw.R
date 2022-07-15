fc_make_mw <-
    function(data_source_mw,
             bin_size = 500,
             Number_of_shifts = 5) {
        util_check_class("data_source_mw", "list")

        util_check_class("bin_size", "numeric")

        util_check_if_integer("bin_size")

        util_check_class("Number_of_shifts", "numeric")

        util_check_if_integer("Number_of_shifts")

        bin_oldest <-
            ceiling(
                max(data_source_mw$age$age)
            )

        bin_youngest <-
            floor(
                min(data_source_mw$age$age)
            )

        shift_value <- bin_size / Number_of_shifts

        bin_breaks <- seq(
            from = bin_youngest,
            to = bin_oldest,
            by = bin_size
        )
        shift_increment <-
            sort(
                rep(0:(Number_of_shifts - 1), length(bin_breaks))
            ) * shift_value

        data.frame(
            name = (rep(bin_breaks, Number_of_shifts) +
                shift_increment),
            shift = sort(rep(c(1:Number_of_shifts), length(bin_breaks)))
        ) %>%
            return()
    }
