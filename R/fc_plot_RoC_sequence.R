#' @title Function to plot the Rate-of-Change sequence
#'
#' @param data_source
#' Data.frame. Output of `fc_estimate_RoC` function
#' @param age_threshold
#' Numeric. Cut-off value used as maximum age.
#' @param Roc_threshold
#' Numeric Cut-off value used as maximum RoC value.
#' @param Peaks
#' Logical. If peak-points are presented in the dataset and `Peaks` == `TRUE`,
#' then peak points will be displayed
#' @param trend
#' If peak-points are presented in the dataset and `Peaks` == `TRUE`,
#'  then one of the three method can be used to visualise the process of peak detection:
#'  \itemize{
#'  \item `"threshold"` - Each point in the RoC sequence is compared to a median
#'  of all RoC scores from the whole sequence (i.e. threshold value).
#'  The ROC value for a point is considered significant if the 95th quantile of
#'  the RoC scores from all calculations is higher than the threshold value.
#'  \item `"trend_linear"` -  A linear model is fitted between the RoC values
#'  and their ages. Differences between the model and each point are calculated
#'  (residuals). The standard deviation (SD) is calculated from all the residuals.
#'  A peak is considered significant if it is 2 SD higher than the model.
#'  \item `"trend_non_linear"` - A conservative generalised additive model (GAM)
#'  is fitted through the RoC scores and their ages (GAM = `RoC ~ s(age, k = 3)` using
#'  the `mgcv` package (Wood, 2011). The distance between each point and the
#'  fitted value is calculated (residuals). The standard deviation (SD) is
#'  calculated from all the residuals. A peak is considered significant if it
#'  is 2 SD higher than the model.
#'  }
#' @description Plot Rate-of-Change sequence with a error estimate and trend
#' and/or peak-points in present.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' example_data <- RRatepol::example_data
#'
#' # example 1
#' sequence_01 <-
#'   fc_estimate_RoC(
#'     data_source_community = example_data$pollen_data[[1]],
#'     data_source_age = example_data$sample_age[[1]],
#'     age_uncertainty = FALSE,
#'     smooth_method = "shep",
#'     Working_Units = "MW",
#'     rand = 1e3,
#'     treads = TRUE,
#'     DC = "chisq"
#'   )
#'
#' fc_plot_RoC_sequence(
#'   sequence_01,
#'   age_threshold = 8e3,
#'   Roc_threshold = 1
#' )
#'
#' # example 2
#' sequence_02 <-
#'   fc_estimate_RoC(
#'     data_source_community = example_data$pollen_data[[2]],
#'     data_source_age = example_data$sample_age[[2]],
#'     age_uncertainty = FALSE,
#'     smooth_method = "shep",
#'     Working_Units = "MW",
#'     rand = 1e3,
#'     treads = TRUE,
#'     DC = "chisq"
#'   )
#'
#' sequence_02_peak <-
#'   fc_detect_peak_points(sequence_01, sel_method = "trend_non_linear")
#'
#' fc_plot_RoC_sequence(
#'   sequence_02_peak,
#'   age_threshold = 8e3,
#'   Roc_threshold = 2,
#'   Peaks = TRUE,
#'   trend = "trend_non_linear"
#' )
#' }
fc_plot_RoC_sequence <-
  function(data_source,
           age_threshold = NULL,
           Roc_threshold = NULLL,
           Peaks = FALSE,
           trend = NULL) {

    # age_threshold
    util_check_class("data_source", "data.frame")

    util_check_col_names(
      "data_source",
      c("Age", "ROC", "ROC_up", "ROC_dw")
    )

    util_check_class("age_threshold", c("NULL", "numeric"))


    if (
      is.null(age_threshold) == TRUE
    ) {
      age_threshold <- max(data_source$Age)
    }

    data_source_filter <-
      data_source %>%
      dplyr::filter(.data$Age <= age_threshold)

    # Roc_threshold
    util_check_class("Roc_threshold", c("NULL", "numeric"))

    if (
      is.null(Roc_threshold) == TRUE
    ) {
      Roc_threshold <- max(data_source$ROC_up)
    }

    util_check_class("Peaks", "logical")

    util_check_class("trend", c("NULL", "character"))

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
        size = 0.1
      ) +
      ggplot2::coord_flip(
        xlim = c(age_threshold, 0),
        ylim = c(0, Roc_threshold)
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
        size = 1,
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

    if (
      Peaks == TRUE
    ) {
      util_check_col_names("data_source", "Peak")

      p_res <-
        p_res +
        ggplot2::geom_point(
          data = .data %>%
            dplyr::filter(Peak == TRUE),
          color = "green",
          alpha = 1,
          size = 1
        )
    }

    if (
      is.null(trend) == FALSE
    ) {
      util_check_vector_values(
        "trend",
        c("threshold", "trend_linear", "trend_non_linear")
      )

      trend <- match.arg(trend)

      if (
        trend == "threshold"
      ) {
        p_res <-
          p_res +
          ggplot2::geom_hline(
            yintercept = stats::median(data_source_filter$ROC),
            color = "blue",
            size = 1
          )
      }

      if (
        trend == "trend_linear"
      ) {
        p_res <-
          p_res +
          ggplot2::geom_line(
            data = data.frame(
              ROC = util_make_trend(
                data_source = data_source,
                sel_method = "linear"
              ),
              Age = data_source$Age
            ),
            color = "blue", size = 1
          )
      }

      if (
        trend == "trend_non_linear"
      ) {
        p_res <-
          p_res +
          ggplot2::geom_line(
            data = data.frame(
              ROC = util_make_trend(
                data_source = data_source,
                sel_method = "non_linear"
              ),
              Age = data_source$Age
            ),
            color = "blue",
            size = 1
          )
      }
    }

    return(p_res)
  }
