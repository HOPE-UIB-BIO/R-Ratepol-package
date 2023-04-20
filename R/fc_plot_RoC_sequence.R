#' @title Function to plot the Rate-of-Change sequence
#'
#' @inheritParams plot_roc
#' @param Roc_threshold
#' `r lifecycle::badge("deprecated")` `Roc_threshold` is no
#'   longer supported; please use `roc_threshold`
#' @param Peaks
#' `r lifecycle::badge("deprecated")` `Peaks` is no
#'   longer supported; please use `peaks`
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated, please use [plot_roc()].
#' @export
#' @seealso [plot_roc()]
#' @keywords internal
fc_plot_RoC_sequence <-
  function(data_source,
           age_threshold = NULL,
           roc_threshold = NULL,
           Roc_threshold = lifecycle::deprecated(),
           peaks = FALSE,
           Peaks = lifecycle::deprecated(),
           trend = NULL) {
    lifecycle::deprecate_warn("1.2.0", "fc_plot_RoC_sequence()", "plot_roc()")

    if (
      lifecycle::is_present(Roc_threshold)
    ) {
      lifecycle::deprecate_warn(
        "1.2.0",
        "estimate_roc(Roc_threshold)",
        "estimate_roc(roc_threshold)"
      )
      roc_threshold <- Roc_threshold
    }

    if (
      lifecycle::is_present(Peaks)
    ) {
      lifecycle::deprecate_warn(
        "1.2.0",
        "estimate_roc(Peaks)",
        "estimate_roc(peaks)"
      )
      peaks <- Peaks
    }

    plot_roc(
      data_source = data_source,
      age_threshold = age_threshold,
      roc_threshold = roc_threshold,
      peaks = peaks,
      trend = trend
    )
  }
