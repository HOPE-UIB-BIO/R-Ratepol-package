#' @title RRatepol: Estimate rate of change
#'
#' @inheritParams estimate_roc
#' @param smooth_N_points
#' `r lifecycle::badge("deprecated")` `smooth_N_points` is no
#'   longer supported; please use `smooth_n_points`
#' @param smooth_N_max
#' `r lifecycle::badge("deprecated")` `smooth_N_max` is no
#'   longer supported; please use `smooth_n_max`
#' @param Working_Units
#' `r lifecycle::badge("deprecated")` `Working_Units` is no
#'   longer supported; please use `working_units`
#' @param Number_of_shifts
#' `r lifecycle::badge("deprecated")` `Number_of_shifts` is no
#'   longer supported; please use `number_of_shifts`
#' @param N_individuals
#' `r lifecycle::badge("deprecated")` `N_individuals` is no
#'   longer supported; please use `n_individuals`
#' @param DC
#' `r lifecycle::badge("deprecated")` `DC` is no
#'   longer supported; please use `dissimilarity_coefficient`
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated, please use [estimate_roc()].
#' @seealso [estimate_roc()]
#' @keywords internal
#' @export 
fc_estimate_RoC <-
  function(data_source_community,
           data_source_age,
           age_uncertainty = NULL,
           smooth_method = c("none", "m.avg", "grim", "age.w", "shep"),
           smooth_n_points = 5,
           smooth_N_points = lifecycle::deprecated(),
           smooth_age_range = 500,
           smooth_n_max = 9,
           smooth_N_max = lifecycle::deprecated(),
           working_units = c("levels", "bins", "MW"),
           Working_Units = lifecycle::deprecated(),
           bin_size = 500,
           number_of_shifts = 5,
           Number_of_shifts = lifecycle::deprecated(),
           bin_selection = c("random", "first"),
           standardise = FALSE,
           n_individuals = 150,
           N_individuals = lifecycle::deprecated(),
           dissimilarity_coefficient = c("euc", "euc.sd", "chord", "chisq", "gower", "bray"),
           DC = lifecycle::deprecated(),
           tranform_to_proportions = TRUE,
           rand = NULL,
           use_parallel = FALSE,
           interest_threshold = NULL,
           time_standardisation = NULL,
           verbose = FALSE) {
    lifecycle::deprecate_warn("1.2.0", "fc_estimate_RoC()", "estimate_roc()")

    if (
      lifecycle::is_present(smooth_N_points)
    ) {
      lifecycle::deprecate_warn(
        "1.2.0",
        "estimate_roc(smooth_N_points)",
        "estimate_roc(smooth_n_points)"
      )
      smooth_n_points <- smooth_N_points
    }

    if (
      lifecycle::is_present(smooth_N_max)
    ) {
      lifecycle::deprecate_warn(
        "1.2.0",
        "estimate_roc(smooth_N_max)",
        "estimate_roc(smooth_n_max)"
      )
      smooth_n_max <- smooth_N_max
    }

    if (
      lifecycle::is_present(Working_Units)
    ) {
      lifecycle::deprecate_warn(
        "1.2.0",
        "estimate_roc(Working_Units)",
        "estimate_roc(working_units)"
      )
      working_units <- Working_Units
    }

    if (
      lifecycle::is_present(Number_of_shifts)
    ) {
      lifecycle::deprecate_warn(
        "1.2.0",
        "estimate_roc(Number_of_shifts)",
        "estimate_roc(number_of_shifts)"
      )
      Number_of_shifts <- Number_of_shifts
    }

    if (
      lifecycle::is_present(N_individuals)
    ) {
      lifecycle::deprecate_warn(
        "1.2.0",
        "estimate_roc(N_individuals)",
        "estimate_roc(n_individuals)"
      )
      n_individuals <- N_individuals
    }

    if (
      lifecycle::is_present(DC)
    ) {
      lifecycle::deprecate_warn(
        "1.2.0",
        "estimate_roc(DC)",
        "estimate_roc(disarity_coefficient)"
      )
      dissimilarity_coefficient <- DC
    }

    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = smooth_method,
      smooth_n_points = smooth_n_points,
      smooth_age_range = smooth_age_range,
      smooth_n_max = smooth_n_max,
      working_units = working_units,
      bin_size = bin_size,
      number_of_shifts = number_of_shifts,
      bin_selection = bin_selection,
      standardise = standardise,
      n_individuals = n_individuals,
      dissimilarity_coefficient = dissimilarity_coefficient,
      tranform_to_proportions = tranform_to_proportions,
      rand = rand,
      use_parallel = use_parallel,
      interest_threshold = interest_threshold,
      time_standardisation = time_standardisation,
      verbose = verbose
    )
  }
