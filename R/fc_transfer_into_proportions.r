
#' @title Transform community data into proportions
#' @param data_source List with `community` data.frame with pollen data.
#'  Each row represent one level (sample) and each column represent one taxon. 
#' able must contain sample_id` as rownames.
#' @param sel_method
#' variable to select result as either proportions (`percentages`) or
#' percentage (`percentages`).
#' @description Tranform pollen data into proportions (or percentages)
fc_transfer_into_proportions  <-
    function(data_source,
             sel_method = c("proportions", "percentages"),
             verbose = FALSE) {
        util_check_class("data_source", "list")

        util_check_class("sel_method", "character")

        util_check_vector_values("sel_method", c("percentages", "proportions"))

        sel_method <- match.arg(sel_method)

        util_check_class("verbose", "logical")

        if (verbose == TRUE) {
            util_output_comment(
                "community data values are being converted to proportions"
            )
        }

        # convert the values community data to proportion of sum of each sample
        data_rownames <-
            rowSums(data_source@community, na.rm = TRUE)

        data_source@community <-
            data_rownames / data_rownames * switch(sel_method,
                "percentages" = 100,
                "proportions" = 1
            )
        return(data_source)
    }
