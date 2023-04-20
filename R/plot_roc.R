#' @title Function to plot the Rate-of-Change sequence
#'
#' @param data_source
#' Data.frame. Output of `estimate_roc` function
#' @param age_threshold
#' Numeric. Cut-off value used as maximum age.
#' @param roc_threshold
#' Numeric Cut-off value used as maximum RoC value.
#' @param peaks
#' Logical. If peak-points are presented in the dataset and `peaks` == `TRUE`,
#' then peak points will be displayed
#' @param trend
#' If peak-points are presented in the dataset and `peaks` == `TRUE`,
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
#'   estimate_roc(
#'     data_source_community = example_data$pollen_data[[1]],
#'     data_source_age = example_data$sample_age[[1]],
#'     age_uncertainty = FALSE,
#'     smooth_method = "shep",
#'     working_units = "MW",
#'     rand = 1e3,
#'     treads = TRUE,
#'     dissimilarity_coefficient = "chisq"
#'   )
#'
#' plot_roc(
#'   sequence_01,
#'   age_threshold = 8e3,
#'   roc_threshold = 1
#' )
#'
#' # example 2
#' sequence_02 <-
#'   estimate_roc(
#'     data_source_community = example_data$pollen_data[[2]],
#'     data_source_age = example_data$sample_age[[2]],
#'     age_uncertainty = FALSE,
#'     smooth_method = "shep",
#'     working_units = "MW",
#'     rand = 1e3,
#'     treads = TRUE,
#'     dissimilarity_coefficient = "chisq"
#'   )
#'
#' sequence_02_peak <-
#'   detect_peak_points(sequence_01, sel_method = "trend_non_linear")
#'
#' plot_roc(
#'   sequence_02_peak,
#'   age_threshold = 8e3,
#'   roc_threshold = 2,
#'   peaks = TRUE,
#'   trend = "trend_non_linear"
#' )
#' }
plot_roc <-
  function(data_source,
           age_threshold = NULL,
           roc_threshold = NULL,
           peaks = FALSE,
           trend = NULL) {

    # age_threshold
    RUtilpol::check_class("data_source", "data.frame")

    RUtilpol::check_col_names(
      "data_source",
      c("Age", "ROC", "ROC_up", "ROC_dw")
    )

    RUtilpol::check_class("age_threshold", c("NULL", "numeric"))


    if (
      isTRUE(is.null(age_threshold))
    ) {
      age_threshold <- max(data_source$Age)
    }

    data_source_filter <-
      data_source %>%
      dplyr::filter(.data$Age <= age_threshold)

    # roc_threshold
    RUtilpol::check_class("roc_threshold", c("NULL", "numeric"))

    if (
      isTRUE(is.null(roc_threshold))
    ) {
      roc_threshold <- max(data_source$ROC_up)
    }

    RUtilpol::check_class("peaks", "logical")

    RUtilpol::check_class("trend", c("NULL", "character"))

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
      isFALSE(is.null(trend))
    ) {
      RUtilpol::check_vector_values(
        "trend",
        c("threshold", "trend_linear", "trend_non_linear")
      )

      if (
        isFALSE(peaks)
      ) {
        RUtilpol::output_comment(
          msg = paste(
            "'trend' has been set to NOT 'NULL',",
            "'peaks' will be plotted"
          )
        )
        # set peaks
        peaks <- TRUE
      }

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
              ROC = make_trend(
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
              ROC = make_trend(
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

    if (
      isTRUE(peaks)
    ) {
      RUtilpol::check_col_names("data_source", "Peak")

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
