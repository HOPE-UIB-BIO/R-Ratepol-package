#' @title Transform community data into proportions
#' @param data_source_trans
#' Data.frame with `label`, `res_age`, and all community data
#' @param sel_method
#' variable to select result as either proportions (`percentages`) or
#' percentage (`percentages`).
#' @param verbose Logical. Should additional information be output?
#' @description Tranform pollen data into proportions (or percentages)
#' @keywords internal
transform_into_proportions <- function(
  data_source_trans,
  sel_method = c("proportions", "percentages"),
  verbose = FALSE,
  silent =FALSE
) {
  RUtilpol::check_class("data_source_trans", "data.frame")

  assertthat::assert_that(
    base::nrow(data_source_trans) > 0,
    msg = "'data_source_trans' must not be empty"
  )

  RUtilpol::check_class("sel_method", "character")

  RUtilpol::check_vector_values(
    "sel_method",
    c("percentages", "proportions")
  )

  assertthat::assert_that(
    base::length(sel_method) == 1,
    msg = "'sel_method' must be a single value"
  )

  sel_method <- match.arg(sel_method)

  RUtilpol::check_class("verbose", "logical")

  RUtilpol::check_class("silent", "logical")

  if (isFALSE(silent) && isTRUE(verbose)) {
    RUtilpol::output_comment(
      "Community data values are being converted to proportions"
    )
  }

  data_com <-
    subset_community(data_source_trans)

  # convert the values community data to proportion of sum of each sample
  data_rowsums <-
    rowSums(data_com, na.rm = TRUE)

  data_com <-
    data_com /
    data_rowsums *
    switch(sel_method, "percentages" = 100, "proportions" = 1)
  data_source_trans[, names(data_com)] <-
    data_com

  return(data_source_trans)
}
