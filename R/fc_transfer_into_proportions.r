
#' @title Transform community data into proportions
#' @param data_source_trans
#' Data.frame with `label`, `res_age`, and all community data
#' @param sel_method
#' variable to select result as either proportions (`percentages`) or
#' percentage (`percentages`).
#' @description Tranform pollen data into proportions (or percentages)
fc_transfer_into_proportions <-
    function(data_source_trans,
             sel_method = c("proportions", "percentages"),
             verbose = FALSE) {
        util_check_class("data_source_trans", "data.frame")

        util_check_class("sel_method", "character")

        util_check_vector_values("sel_method", c("percentages", "proportions"))

        sel_method <- match.arg(sel_method)

        util_check_class("verbose", "logical")

        if (verbose == TRUE) {
            util_output_comment(
                "Community data values are being converted to proportions"
            )
        }

        data_com <-
            util_subset_community(data_source_trans)

        # convert the values community data to proportion of sum of each sample
        data_rowsums <-
            rowSums(data_com, na.rm = TRUE)

        data_com <-
            data_com / data_rowsums *
                switch(sel_method,
                    "percentages" = 100,
                    "proportions" = 1
                )
        data_source_trans[, names(data_com)] <-
            data_com

        return(data_source_trans)
    }
