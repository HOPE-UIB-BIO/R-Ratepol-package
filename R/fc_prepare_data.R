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

        if (
            Working_Units != "levels"
        ) {
            util_check_class("bin_size", "numeric")

            util_check_if_integer("bin_size")

            if (
                Working_Units == "MW"
            ) {
                util_check_class("Number_of_shifts", "numeric")

                util_check_if_integer("Number_of_shifts")
            }
        }

        util_check_class("rand", c("NULL", "numeric"))

        if (
            is.null(rand) == FALSE
        ) {
            util_check_if_integer("rand")
        }

        res <-
            switch(Working_Units,
                "levels" = {
                    fc_prepare_data_levels(
                        data_source_prep = data_source_prep,
                        rand = rand
                    )
                },
                "bins" = {
                    fc_prepare_data_bins(
                        data_source_prep = data_source_prep,
                        bin_size = bin_size,
                        Number_of_shifts = 0,
                        rand = rand
                    )
                },
                "MW" = {
                    fc_prepare_data_bins(
                        data_source_prep = data_source_prep,
                        bin_size = bin_size,
                        Number_of_shifts = Number_of_shifts,
                        rand = rand
                    )
                }
            )

        return(res)
    }
