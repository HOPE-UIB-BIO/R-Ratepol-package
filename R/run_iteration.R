#' @title Run a single iteration of RoC estimation
#'
#' @param data_source_run
#' `list` with `data` and `bins` prepared by [prepare_data()].
#' @inheritParams estimate_roc
#' @description
#' A single run is computed following these steps:
#' \itemize{
#' \item Subset levels into each bin to define Working Units.
#' \item Optionally standardise assemblage counts in each WU.
#' \item Calculate dissimilarity between consecutive WUs.
#' \item Standardise dissimilarity by age difference between WUs.
#' }
#' @return
#' A `tibble` with one row per WU pair containing columns `Age`, `ROC`,
#' `ROC_up`, `ROC_dw`, and `Working_Unit`.
#' @seealso [estimate_roc()]
#' @keywords internal
run_iteration <- function(
  data_source_run,
  bin_selection = "first",
  standardise = FALSE,
  n_individuals = 150,
  tranform_to_proportions = TRUE,
  dissimilarity_coefficient = "euc",
  time_standardisation = 500,
  verbose = FALSE,
  silent = FALSE
) {
  util_check_class(data_source_run, "list")

  util_check_class(standardise, "logical")

  assertthat::assert_that(
    base::length(standardise) == 1,
    msg = "'standardise' must be a single TRUE or FALSE"
  )

  util_check_class(tranform_to_proportions, "logical")

  assertthat::assert_that(
    base::length(tranform_to_proportions) == 1,
    msg = "'tranform_to_proportions' must be a single TRUE or FALSE"
  )

  util_check_class(verbose, "logical")

  assertthat::assert_that(
    base::length(verbose) == 1,
    msg = "'verbose' must be a single TRUE or FALSE"
  )

  assertthat::assert_that(
    !base::is.na(time_standardisation) && time_standardisation != 0,
    msg = "'time_standardisation' must not be 0 or NA"
  )

  #----------------------------------------------------------#
  # 4.1 Data subsetting -----
  #----------------------------------------------------------#

  # select one sample for each bin based on the age of the samples.
  # subset data
  data_subset <-
    subset_samples(
      data_source_subset = data_source_run$data,
      data_source_bins = data_source_run$bins,
      bin_selection = bin_selection
    )

  # reduce
  data_subset <-
    reduce_data_simple(
      data_source_reduce = data_subset
    )

  #----------------------------------------------------------#
  # 4.2 Data Standardisation -----
  #----------------------------------------------------------#

  # standardisation of community data to n_individuals
  if (isTRUE(standardise)) {
    # select only community data
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )

    # adjust the value to a minimal of presented values
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )

    # check if all samples has n_individuals of individuals
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]

    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )

    # standardisation
    data_sd <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      )

    assertthat::assert_that(
      all(
        n_individuals ==
          rowSums(
            subset_community(data_sd),
            na.rm = TRUE
          )
      ),
      msg = paste(
        "Data standardisation was unsuccessful,",
        "try 'standardise' = FALSE"
      )
    )
  } else {
    data_sd <- data_subset
  }

  # data reduce
  data_sd <-
    reduce_data_simple(
      data_source_reduce = data_sd
    )

  if (isTRUE(tranform_to_proportions)) {
    # tunr into proportion
    data_sd_prop <-
      transform_into_proportions(
        data_source_trans = data_sd,
        sel_method = "proportions",
        verbose = verbose,
        silent = silent
      )
  } else {
    data_sd_prop <-
      data_sd
  }

  #----------------------------------------------------------#
  # 4.3 dissimilarity_coefficient Calculation -----
  #----------------------------------------------------------#

  # calculate dissimilarity_coefficient between each subsequent samples/bins
  dc_res <-
    estimate_dissimilarity_coefficient(
      data_source_dc = data_sd_prop,
      dissimilarity_coefficient = dissimilarity_coefficient,
      verbose = verbose,
      silent = silent
    )

  #----------------------------------------------------------#
  # 4.4 Rate of Change -----
  #----------------------------------------------------------#

  #  calculate dissimilarity_coefficient standardise by time
  roc_res <-
    data_sd_prop[seq_along(dc_res), ] %>%
    dplyr::mutate(
      dc = dc_res,
      age_diff_st = .data$age_diff / time_standardisation,
      roc = .data$dc / .data$age_diff_st
    ) %>%
    dplyr::select("label", "res_age", "roc")

  #----------------------------------------------------------#
  # 4.6 Result of a single randomisation run -----
  #----------------------------------------------------------#

  if (nrow(roc_res) < 1) {
    stop("Estimation not succesfull")
  }

  return(roc_res)
}
