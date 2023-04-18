#' @title Detect significant peak points
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated, please use [detect_peak_points()].
#' @export
#' @seealso [detect_peak_points()]
#' @keywords internal
fc_detect_peak_points <-
  function(data_source,
           sel_method = c(
             "trend_linear", "trend_non_linear",
             "threshold", "GAM_deriv", "SNI"
           ),
           sd_threshold = 2) {
    detect_peak_points(
      data_source = data_source,
      sel_method = sel_method,
      sd_threshold = sd_threshold
    )
  }
