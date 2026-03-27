#' @title Transform community data into proportions
#' @param data_source_trans
#' `data.frame` with `label`, `res_age`, and community count columns.
#' @param sel_method
#' `character`. Scale of the output:
#' \itemize{
#' \item `"proportions"` - values sum to 1 per sample.
#' \item `"percentages"` - values sum to 100 per sample.
#' }
#' @param verbose `logical`. If `TRUE`, print progress messages.
#' @description
#' Transform community count data into proportions or percentages.
#' @return
#' The input `data.frame` with community columns replaced by the
#' transformed values.
#' @keywords internal
transform_into_proportions <- function(
  data_source_trans,
  sel_method = c("proportions", "percentages"),
  verbose = FALSE,
  silent =FALSE
) {
  util_check_class(data_source_trans, "data.frame")

  assertthat::assert_that(
    base::nrow(data_source_trans) > 0,
    msg = "'data_source_trans' must not be empty"
  )

  util_check_class(sel_method, "character")

  util_check_vector_values(
    sel_method,
    c("percentages", "proportions")
  )

  assertthat::assert_that(
    base::length(sel_method) == 1,
    msg = "'sel_method' must be a single value"
  )

  sel_method <- match.arg(sel_method)

  util_check_class(verbose, "logical")

  util_check_class(silent, "logical")

  if (isFALSE(silent) && isTRUE(verbose)) {
    util_output_comment(
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
