#' @title Detect significant peak points
#' @description
#' Detect points of sudden increase in Rate-of-Change values from the
#' output of [estimate_roc()].
#' @param data_source
#' `tibble`. Output of [estimate_roc()].
#' @param sel_method
#' `character`. Method to use for peak-point detection.
#' \itemize{
#' \item `"threshold"` - Each point is compared to the median of all RoC
#' scores. A point is significant if its 95th-quantile RoC exceeds the
#' median threshold.
#' \item `"trend_linear"` - A linear model is fitted between RoC values
#' and their ages. A peak is significant if it is `sd_threshold` SD
#' above the fitted value.
#' \item `"trend_non_linear"` - A conservative GAM
#' (`RoC ~ s(age, k = 3)`) is fitted. A peak is significant if it is
#' `sd_threshold` SD above the fitted value.
#' \item `"GAM_deriv"` - A smooth GAM (`RoC ~ s(age)`) is fitted and the
#' first derivative evaluated using the `gratia` package (Simpson, 2018).
#' A peak is significant if the confidence interval of the first
#' derivative excludes zero.
#' \item `"SNI"` - Signal-to-noise index adapted from Kelly et al.
#' (2011). A peak is significant if SNI > 3.
#' }
#' @param sd_threshold
#' `numeric`. Number of standard deviations above the trend required for
#' a point to be classified as a peak (default = 2). Used by
#' `"trend_linear"` and `"trend_non_linear"`.
#' @return
#' The input `tibble` with an additional `logical` column `Peak` that
#' is `TRUE` for samples identified as significant peak points.
#' @seealso [estimate_roc()], [plot_roc()]
#' @references
#' Kelly, R.F., Higuera, P.E., Barrett, C.M., Feng Sheng, H., 2011.
#' A signal-to-noise index to quantify the potential for peak detection
#' in sediment-charcoal records. Quat. Res. 75, 11-17.
#'
#' Simpson, G.L., 2018. Modelling palaeoecological time series using
#' generalised additive models. Front. Ecol. Evol. 6, 1-21.
#'
#' Wood, S.N., 2011. Fast stable restricted maximum likelihood and
#' marginal likelihood estimation of semiparametric generalized linear
#' models. J. R. Stat. Soc. Ser. B Stat. Methodol. 73, 3-36.
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
#' sequence_01_peak <-
#'   detect_peak_points(
#'     sequence_01,
#'     sel_method = "trend_non_linear",
#'     sd_threshold = 2
#'   )
#'
#' plot_roc(
#'   sequence_01_peak,
#'   age_threshold = 8e3,
#'   roc_threshold = 2,
#'   peaks = TRUE,
#'   trend = "trend_non_linear"
#' )
#' }
detect_peak_points <-
  function(data_source,
           sel_method = c(
             "trend_linear", "trend_non_linear",
             "threshold", "GAM_deriv", "SNI"
           ),
           sd_threshold = 2) {
    util_check_class(data_source, "data.frame")

    util_check_col_names(data_source, "ROC")

    util_check_class(sel_method, "character")

    util_check_vector_values(
      sel_method,
      c(
        "trend_linear", "trend_non_linear",
        "threshold", "GAM_deriv", "SNI"
      )
    )

    sel_method <- match.arg(sel_method)

    util_check_class(sd_threshold, "numeric")

    assertthat::assert_that(
      sd_threshold > 0,
      msg = "'sd_threshold' must be bigger than 0"
    )

    #----------------------------------------------------------#
    # 1. Median peak threshold -----
    #----------------------------------------------------------#
    if (
      sel_method == "threshold"
    ) {
      util_check_col_names(data_source, "ROC_dw")

      # threshold for RoC peaks is set as median of all RoC in dataset
      r_threshold <-
        stats::median(data_source$ROC)

      # mark peaks which have 95% quantile above the threshold as Peak
      data_source$Peak <-
        data_source$ROC_dw > r_threshold
    }

    #----------------------------------------------------------#
    # 2. Linear trend  -----
    #----------------------------------------------------------#
    if (
      sel_method == "trend_linear"
    ) {
      util_check_col_names(data_source, "Age")

      # mark points that are abowe the linear model
      #   (exactly sd_threshold SD higher than prediction)
      data_source$pred_linear <-
        make_trend(
          data_source = data_source,
          sel_method = "linear"
        )

      data_source$residuals <-
        data_source$ROC - data_source$pred_linear

      data_source$Peak <-
        (data_source$residuals) >
          (sd_threshold * stats::sd(data_source$residuals))
    }

    #----------------------------------------------------------#
    # 3. Non-linear trend  -----
    #----------------------------------------------------------#
    if (
      sel_method == "trend_non_linear"
    ) {
      util_check_col_names(data_source, "Age")
      # mark points that are abowe the GAM model
      #   (exactly sd_threshold SD higher than GAM prediction)
      data_source$pred_gam <-
        make_trend(
          data_source = data_source,
          sel_method = "non_linear"
        )

      data_source$residuals <-
        data_source$ROC - data_source$pred_gam

      data_source$Peak <-
        (data_source$residuals) >
          (sd_threshold * stats::sd(data_source$residuals))
    }

    #----------------------------------------------------------#
    # 4. Firts derivative of GAM model  -----
    #----------------------------------------------------------#
    if (
      sel_method == "GAM_deriv"
    ) {
      util_check_col_names(data_source, "Age")
      # fit gam well smother gam model and use first derivative of the function
      #   to detect signifiant increases in the function
      gam_model <-
        mgcv::gam(
          ROC ~ s(Age),
          data = data_source,
          family = mgcv::Tweedie(p = 2),
          method = "REML"
        )

      new_data <-
        tibble::tibble(Age = data_source$Age)

      gam_deriv <-
        gratia::derivatives(gam_model,
          data = new_data,
          n = 1000
        )

      data_source$Peak <-
        (gam_deriv$.lower_ci > 0)
    }

    #----------------------------------------------------------#
    # 5. Signal-to-Noise-ratio Index  -----
    #----------------------------------------------------------#
    if (
      sel_method == "SNI"
    ) {
      util_check_col_names(data_source, "Age")
      # set moving window of 5 times higher than average distance between samples
      mean_age_window <- 5 * mean(diff(data_source$Age))

      # create GAM
      pred_gam <-
        make_trend(
          data_source = data_source,
          sel_method = "non_linear"
        )

      # calculate SNI (singal to noise ratio)
      SNI_calc <-
        detect_sni(
          data.frame(
            data_source$Age,
            data_source$ROC,
            pred_gam
          ),
          mean_age_window
        )

      # mark points with SNI higher than 3
      data_source$Peak <-
        (SNI_calc$SNI > 3) & (data_source$ROC > pred_gam)
    }

    #----------------------------------------------------------#
    # 6. save result  -----
    #----------------------------------------------------------#

    data_result <-
      data_source %>%
      dplyr::select(
        "Working_Unit",
        "Age",
        "ROC",
        "ROC_up",
        "ROC_dw",
        "Peak"
      )

    return(data_result)
  }
