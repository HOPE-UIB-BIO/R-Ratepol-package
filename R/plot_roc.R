#' @title Plot the Rate-of-Change sequence
#' @description
#' Plot RoC scores through time with an upper uncertainty envelope and,
#' optionally, a superimposed trend curve and peak-point markers.
#' @param data_source
#' `tibble`. Output of [estimate_roc()] or [detect_peak_points()].
#' @param age_threshold
#' `numeric`. Optional. Upper (oldest) age cut-off; samples older than
#' this value are excluded from the plot.
#' @param roc_threshold
#' `numeric`. Optional. Upper RoC cut-off; values above this are clipped.
#' @param peaks
#' `logical`. If `TRUE` and a `Peak` column is present in `data_source`,
#' peak points are highlighted on the plot (default = `FALSE`).
#' @param trend
#' `character` or `NULL`. When `peaks = TRUE`, optionally overlay the
#' trend curve used during peak detection. One of `"threshold"`,
#' `"trend_linear"`, or `"trend_non_linear"`. `NULL` (default) shows no
#' trend line.
#' @param silent
#' `logical`. If `TRUE`, suppress all console output (default = `FALSE`).
#' @return
#' A `ggplot2` object.
#' @seealso [estimate_roc()], [detect_peak_points()]
#' @export
#' @examples
#' \dontrun{
#' data("example_data", package = "RRatepol")
#'
#' sequence_01 <-
#'   estimate_roc(
#'     data_source_community = example_data$pollen_data[[1]],
#'     data_source_age = example_data$sample_age[[1]],
#'     smooth_method = "shep",
#'     working_units = "MW",
#'     rand = 1e3,
#'     use_parallel = TRUE,
#'     dissimilarity_coefficient = "chisq"
#'   )
#'
#' plot_roc(
#'   sequence_01,
#'   age_threshold = 8e3,
#'   roc_threshold = 1
#' )
#' }
plot_roc <- function(
  data_source,
  age_threshold = NULL,
  roc_threshold = NULL,
  peaks = FALSE,
  trend = NULL,
  silent = FALSE
) {
  # age_threshold
  util_check_class(data_source, "data.frame")

  util_check_col_names(
    data_source,
    c("Age", "ROC", "ROC_up", "ROC_dw")
  )

  util_check_class(age_threshold, c("NULL", "numeric"))

  if (isTRUE(is.null(age_threshold))) {
    age_threshold <- max(data_source$Age)
  }

  data_source_filter <-
    data_source %>%
    dplyr::filter(.data$Age <= age_threshold)

  if (base::any(base::is.na(data_source_filter$ROC))) {
    warning(
      "NA values detected in 'ROC' column of 'data_source'.",
      call. = FALSE
    )
  }

  # roc_threshold
  util_check_class(roc_threshold, c("NULL", "numeric"))

  if (isTRUE(is.null(roc_threshold))) {
    roc_threshold <- max(data_source$ROC_up)
  }

  util_check_class(peaks, "logical")

  assertthat::assert_that(
    base::length(peaks) == 1,
    msg = "'peaks' must be a single value"
  )

  assertthat::assert_that(
    !base::is.na(peaks),
    msg = "'peaks' must not be NA"
  )

  util_check_class(trend, c("NULL", "character"))

  p_res <-
    ggplot2::ggplot(
      data_source_filter,
      mapping = ggplot2::aes(
        y = .data$ROC,
        x = .data$Age
      )
    ) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_continuous(trans = "reverse") +
    ggplot2::geom_vline(
      xintercept = seq(0, age_threshold, 2e3),
      colour = "gray90",
      linewidth = 0.1
    ) +
    ggplot2::coord_flip(
      xlim = c(age_threshold, 0),
      ylim = c(0, roc_threshold)
    ) +
    ggplot2::geom_ribbon(
      mapping = ggplot2::aes(
        ymin = .data$ROC_up,
        ymax = .data$ROC_dw
      ),
      fill = "gray90"
    ) +
    ggplot2::geom_line(
      alpha = 1,
      linewidth = 1,
      color = "gray30"
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "gray30",
      lty = 3
    ) +
    ggplot2::labs(
      x = "Age (cal yr BP)",
      y = "Rate of change score"
    )

  if (isFALSE(is.null(trend))) {
    util_check_vector_values(
      trend,
      c("threshold", "trend_linear", "trend_non_linear")
    )

    if (isFALSE(peaks)) {
      if (isFALSE(silent)) {
        util_output_comment(
          msg = paste(
            "'trend' has been set to NOT 'NULL',",
            "'peaks' will be plotted"
          )
        )
      }
      # set peaks
      peaks <- TRUE
    }

    if (trend == "threshold") {
      p_res <-
        p_res +
        ggplot2::geom_hline(
          yintercept = stats::median(data_source_filter$ROC),
          color = "blue",
          linewidth = 1
        )
    }

    if (trend == "trend_linear") {
      p_res <-
        p_res +
        ggplot2::geom_line(
          data = data.frame(
            ROC = make_trend(
              data_source = data_source,
              sel_method = "linear"
            ),
            Age = data_source$Age
          ),
          color = "blue",
          linewidth = 1
        )
    }

    if (trend == "trend_non_linear") {
      p_res <-
        p_res +
        ggplot2::geom_line(
          data = data.frame(
            ROC = make_trend(
              data_source = data_source,
              sel_method = "non_linear"
            ),
            Age = data_source$Age
          ),
          color = "blue",
          linewidth = 1
        )
    }
  }

  if (isTRUE(peaks)) {
    util_check_col_names(data_source, "Peak")

    p_res <-
      p_res +
      ggplot2::geom_point(
        data = data_source_filter %>%
          dplyr::filter(.data$Peak == TRUE),
        color = "green",
        alpha = 1,
        size = 3
      )
  }

  return(p_res)
}
