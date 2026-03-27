#' @title Estimate rate of change
#' @description
#' Estimate the Rate of Change (RoC) in community composition along a
#' temporal sequence. RoC is defined as the dissimilarity between consecutive
#' Working Units (WUs), standardised by the age difference between them.
#' @param data_source_community
#' `data.frame`. Community data with taxa as columns and samples as rows.
#' The first column must be named `sample_id` (`character`).
#' @param data_source_age
#' `data.frame` with two columns:
#' \itemize{
#' \item `sample_id` - unique identifier of each level (`character`)
#' \item `age` - age of the level (`numeric`)
#' }
#' @param age_uncertainty
#' Optional age-uncertainty matrix from an age-depth model. Either:
#' \itemize{
#' \item A `matrix` with one column per sample and one row per age sequence
#' drawn from a posterior age-depth model. One row is randomly selected at
#' the start of each randomisation run.
#' \item `NULL` (default) - age uncertainties are not used.
#' }
#' @param smooth_method
#' `character`. Smoothing method applied to each taxon before RoC is
#' computed.
#' \itemize{
#' \item `"none"` - no smoothing (default)
#' \item `"shep"` - Shepard's 5-term filter (Davis, 1986; Wilkinson, 2005)
#' \item `"m.avg"` - moving average
#' \item `"age.w"` - age-weighted average
#' \item `"grim"` - Grimm's smoothing (Grimm & Jacobson, 1992)
#' }
#' @param smooth_n_points
#' `numeric`. Number of points used for moving average, Grimm, and
#' age-weighted smoothing. Must be an odd number.
#' @param smooth_age_range
#' `numeric`. Maximum age range (in years) for Grimm and age-weighted
#' smoothing windows.
#' @param smooth_n_max
#' `numeric`. Maximum number of samples included in a Grimm smoothing
#' window.
#' @param working_units
#' `character`. Strategy used to define Working Units between which
#' dissimilarity is calculated.
#' \itemize{
#' \item `"levels"` - each stratigraphical level is its own WU.
#' \item `"bins"` - one representative level is selected from each time
#' bin of width `bin_size`.
#' \item `"MW"` - moving-window binning: selective binning is repeated
#' `number_of_shifts` times, shifting the window by
#' `bin_size / number_of_shifts` years each time. All results are
#' retained and summarised together.
#' }
#' @param bin_size
#' `numeric`. Width of each time bin in years. Used when `working_units`
#' is `"bins"` or `"MW"`.
#' @param number_of_shifts
#' `numeric`. Number of window shifts in moving-window binning
#' (`working_units = "MW"`).
#' @param bin_selection
#' `character`. Rule for selecting one level from each bin.
#' \itemize{
#' \item `"random"` (default) - a level is selected at random.
#' \item `"first"` - the level closest to the start of the bin is
#' selected.
#' }
#' @param standardise
#' `logical`. If `TRUE`, assemblage counts in each WU are rarefied to
#' `n_individuals` before dissimilarity is computed.
#' @param n_individuals
#' `numeric`. Number of individuals to rarefy to when
#' `standardise = TRUE`. Automatically reduced to the smallest count in
#' the sequence if any WU has fewer individuals.
#' @param dissimilarity_coefficient
#' `character`. Dissimilarity coefficient used to compare consecutive WUs.
#' See `vegan::vegdist()` for details.
#' \itemize{
#' \item `"euc"` - Euclidean distance
#' \item `"euc.sd"` - standardised Euclidean distance
#' \item `"chord"` - Chord distance
#' \item `"chisq"` - Chi-squared coefficient
#' \item `"gower"` - Gower's distance
#' \item `"bray"` - Bray-Curtis dissimilarity
#' }
#' @param tranform_to_proportions
#' `logical`. If `TRUE` (default), community data are converted to
#' proportions before dissimilarity is computed.
#' @param rand
#' `numeric`. Number of randomisation runs. Set to `NULL` (default) to
#' skip randomisation and use a single deterministic run.
#' @param use_parallel
#' Controls parallel computation of randomisation runs.
#' \itemize{
#' \item `FALSE` (default) - single core only.
#' \item `TRUE` - number of cores detected automatically.
#' \item A positive `numeric` - use that many cores.
#' }
#' @param interest_threshold
#' `numeric`. Optional. Samples older than this age are excluded from
#' the output.
#' @param time_standardisation
#' `numeric`. Time unit for RoC values. RoC is reported as dissimilarity
#' per `time_standardisation` years. Defaults to `bin_size` when `NULL`.
#' @param verbose
#' `logical`. If `TRUE`, print progress messages during computation.
#' @param silent
#' `logical`. If `TRUE`, suppress all console output (overrides
#' `verbose`).
#' @return
#' A `tibble` with one row per Working Unit pair, containing columns:
#' \itemize{
#' \item `Age` - mean age of the WU pair
#' \item `ROC` - median RoC score across all randomisation runs
#' \item `ROC_up` - 95th-quantile RoC score (upper uncertainty bound)
#' }
#' @details
#' RoC between two consecutive WUs is computed as dissimilarity divided
#' by the age difference between the WUs, scaled by
#' `time_standardisation`. When `rand > 1`, the full computation is
#' repeated `rand` times; in each run, one age sequence is drawn at
#' random from `age_uncertainty` (if supplied) and, when
#' `standardise = TRUE`, assemblage counts are independently rarefied.
#' The final RoC value for each WU pair is the median across all runs;
#' the 95th quantile is returned as an upper uncertainty bound.
#' @seealso
#' [detect_peak_points()], [plot_roc()]
#' @references
#' Birks, H.J.B., Gordon, A.D., 1985. Numerical Methods in Quaternary
#' Pollen Analysis. Academic Press, London.
#'
#' Davis, J.C., 1986. Statistics and Data Analysis in Geology, 2nd edn.
#' J. Wiley & Sons, New York.
#'
#' Grimm, E.C., Jacobson, G.L., 1992. Fossil-pollen evidence for abrupt
#' climate changes during the past 18000 years in eastern North America.
#' Clim. Dyn. 6, 179-184.
#'
#' Wilkinson, L., 2005. The Grammar of Graphics. Springer-Verlag,
#' New York.
#' @export
#' @examples
#' \dontrun{
#' data("example_data", package = "RRatepol")
#'
#' # `rand = NULL` uses a single deterministic run. For robust results,
#' # increase to e.g. `rand = 1e3` with `use_parallel = TRUE`.
#' sequence_01 <-
#'   estimate_roc(
#'     data_source_community = example_data$pollen_data[[1]],
#'     data_source_age = example_data$sample_age[[1]],
#'     dissimilarity_coefficient = "chisq",
#'     rand = NULL # increase to e.g. `rand = 1e3` with `use_parallel = TRUE` for robust results
#'   )
#'
#' plot_roc(
#'   data_source = sequence_01,
#'   age_threshold = 8e3,
#'   roc_threshold = 1
#' )
#' }
estimate_roc <- function(
  data_source_community,
  data_source_age,
  age_uncertainty = NULL,
  smooth_method = c("none", "m.avg", "grim", "age.w", "shep"),
  smooth_n_points = 5,
  smooth_age_range = 500,
  smooth_n_max = 9,
  working_units = c("levels", "bins", "MW"),
  bin_size = 500,
  number_of_shifts = 5,
  bin_selection = c("random", "first"),
  standardise = FALSE,
  n_individuals = 150,
  dissimilarity_coefficient = c(
    "euc",
    "euc.sd",
    "chord",
    "chisq",
    "gower",
    "bray"
  ),
  tranform_to_proportions = TRUE,
  rand = NULL,
  use_parallel = FALSE,
  interest_threshold = NULL,
  time_standardisation = NULL,
  verbose = FALSE,
  silent = FALSE
) {
  # Start of the code

  #----------------------------------------------------------#
  # 0. Arguments check -----
  #----------------------------------------------------------#

  assertthat::assert_that(
    !missing(data_source_community),
    msg = "Object 'data_source_community' must be included as a 'data.frame'"
  )

  assertthat::assert_that(
    !missing(data_source_age),
    msg = "Object 'data_source_age' must be included as a 'data.frame'"
  )

  util_check_class(data_source_community, "data.frame")

  util_check_class(data_source_age, "data.frame")

  util_check_class(age_uncertainty, c("NULL", "matrix"))

  util_check_class(working_units, "character")

  util_check_vector_values(working_units, c("levels", "bins", "MW"))

  working_units <- match.arg(working_units)

  if (is.null(time_standardisation)) {
    time_standardisation <- bin_size
  }
  util_check_class(time_standardisation, "numeric")

  util_check_if_integer(time_standardisation)

  if (working_units != "levels") {
    util_check_class(bin_size, "numeric")

    util_check_if_integer(bin_size)

    util_check_class(bin_selection, "character")

    util_check_vector_values(bin_selection, c("first", "random"))

    bin_selection <- match.arg(bin_selection)

    if (working_units == "MW") {
      util_check_class(number_of_shifts, "numeric")

      util_check_if_integer(number_of_shifts)
    }
  }

  util_check_class(standardise, "logical")

  assertthat::assert_that(
    base::length(standardise) == 1,
    msg = "'standardise' must be a single TRUE or FALSE"
  )

  if (isTRUE(standardise)) {
    util_check_class(n_individuals, "numeric")

    util_check_if_integer(n_individuals)
  }

  util_check_class(tranform_to_proportions, "logical")

  assertthat::assert_that(
    base::length(tranform_to_proportions) == 1,
    msg = "'tranform_to_proportions' must be a single TRUE or FALSE"
  )

  util_check_class(interest_threshold, c("NULL", "numeric"))

  if (!base::is.null(interest_threshold)) {
    assertthat::assert_that(
      base::length(interest_threshold) == 1,
      msg = "'interest_threshold' must be a single value"
    )
  }

  util_check_class(smooth_method, "character")

  util_check_vector_values(
    smooth_method,
    c("none", "m.avg", "grim", "age.w", "shep")
  )

  smooth_method <- match.arg(smooth_method)

  if (!smooth_method %in% c("none", "shep")) {
    assertthat::assert_that(
      smooth_n_points %% 2 != 0,
      msg = "'smooth_n_points' must be an odd number"
    )

    if (smooth_method != "m.avg") {
      util_check_class(smooth_age_range, "numeric")

      if (smooth_method == "grim") {
        assertthat::assert_that(
          smooth_n_max %% 2 != 0,
          msg = "'smooth_n_max' must be an odd number"
        )

        assertthat::assert_that(
          smooth_n_points < smooth_n_max,
          msg = "'smooth_n_max' must be bigger than 'smooth_n_points"
        )
      }
    }
  }

  util_check_class(dissimilarity_coefficient, "character")

  util_check_vector_values(
    dissimilarity_coefficient,
    c("euc", "euc.sd", "chord", "chisq", "gower", "bray")
  )

  dissimilarity_coefficient <- match.arg(dissimilarity_coefficient)

  util_check_class(rand, c("NULL", "numeric"))

  if (isFALSE(is.null(rand))) {
    util_check_if_integer(rand)
  }

  util_check_class(use_parallel, c("logical", "numeric"))

  if (is.numeric(use_parallel)) {
    util_check_if_integer(use_parallel)

    assertthat::assert_that(
      !base::is.na(use_parallel) && use_parallel != 0,
      msg = "'use_parallel' must not be 0 or NA when numeric"
    )
  }

  util_check_class(verbose, "logical")

  assertthat::assert_that(
    base::length(verbose) == 1,
    msg = "'verbose' must be a single TRUE or FALSE"
  )

  util_check_class(silent, "logical")

  #--------------------------------------------------#
  # 0.1. Report to user -----
  #--------------------------------------------------#

  start_time <- Sys.time()

  if (isFALSE(silent)) {
    util_output_heading(
      paste("RRatepol started", start_time),
      size = "h1"
    )
  }

  if (isFALSE(is.null(age_uncertainty))) {
    if (isFALSE(silent)) {
      util_output_comment(
        "'age_uncertainty' will be used for in the RoC estimation"
      )
    }

    if (rand < 100) {
      if (isFALSE(silent)) {
        util_output_warning(
          paste(
            "'age_uncertainty' was selected to be used with low number",
            "of replication. Recommend to increase 'rand'"
          )
        )
      }
    }
  }

  if (isFALSE(silent)) {
    switch(
      working_units,
      "levels" = {
        util_output_comment(
          "RoC will be estimated between individual subsequent levels"
        )
      },
      "bins" = {
        util_output_comment(
          paste(
            "RoC will be estimated using selective binning with",
            bin_size,
            "yr time bin"
          )
        )
      },
      "MW" = {
        util_output_comment(
          paste(
            "RoC will be estimated using 'binning with the mowing window' of",
            bin_size,
            "yr time bin over",
            number_of_shifts,
            "number of window shifts"
          )
        )
      }
    )
  }

  if (working_units != "levels") {
    if (bin_selection == "random") {
      if (isFALSE(silent)) {
        util_output_comment(
          "Sample will randomly selected for each bin"
        )
      }

      if (rand < 100) {
        if (isFALSE(silent)) {
          util_output_warning(
            paste(
              "'bin_selection' was selected as 'random' with low number",
              "of replication. Recommend to increase 'rand'"
            )
          )
        }
      }
    } else {
      if (isFALSE(silent)) {
        util_output_comment(
          "First sample of each time bin will selected"
        )
      }
    }
  }

  if (isFALSE(silent)) {
    util_output_comment(
      paste(
        "'time_standardisation' =",
        time_standardisation,
        ":",
        "RoC values will be reported as disimilarity per",
        time_standardisation,
        "years."
      )
    )
  }

  if (working_units != "levels" && time_standardisation != bin_size) {
    if (isFALSE(silent)) {
      util_output_comment(
        paste(
          "RoC values will be reported in different units than size of bin.",
          "Recommend to keep 'time_standardisation'",
          "and 'bin_size' as same values"
        )
      )
    }
  }

  if (isTRUE(standardise)) {
    if (isFALSE(silent)) {
      util_output_comment(
        paste(
          "Data will be standardise in each Working unit to",
          n_individuals,
          "or the lowest number detected in dataset"
        )
      )
    }

    if (rand < 100) {
      if (isFALSE(silent)) {
        util_output_warning(
          paste(
            "'standardise' was selected as 'TRUE' with low number of replication.",
            "Recommend to increase 'rand'"
          )
        )
      }
    }
  }

  #----------------------------------------------------------#
  # 1. Data extraction -----
  #----------------------------------------------------------#

  # extract data into working format
  # already include data check
  data_extract <-
    extract_data(
      data_community_extract = data_source_community,
      data_age_extract = data_source_age,
      age_uncertainty = age_uncertainty,
      verbose = verbose,
      silent = silent
    )

  if (ncol(data_extract$community) == 1 && isTRUE(tranform_to_proportions)) {
    if (isFALSE(silent)) {
      util_output_warning(
        msg = paste(
          "Community data has only 1 variable and `tranform_to_proportions`",
          "is set to `TRUE`.",
          "This will result in 0 RoC.",
          "Therefore, `tranform_to_proportions` will be set to `FALSE`"
        )
      )
    }

    tranform_to_proportions <- FALSE
  }

  #----------------------------------------------------------#
  # 2. Data smoothing -----
  #----------------------------------------------------------#

  if (smooth_method != "none") {
    # smooth data by selected smoothing type
    data_smooth <-
      smooth_community_data(
        data_source_smooth = data_extract,
        smooth_method = smooth_method,
        smooth_n_points = smooth_n_points,
        smooth_n_max = smooth_n_max,
        smooth_age_range = smooth_age_range,
        round_results = standardise,
        verbose = verbose,
        silent = silent
      )
  } else {
    data_smooth <- data_extract
  }

  # reduce data dimentions
  data_work <-
    reduce_data(
      data_source_reduce = data_smooth
    )

  if (isFALSE(silent) && isTRUE(verbose)) {
    check_data(data_work, silent = silent)
  }

  #----------------------------------------------------------#
  # 3. Crete datasets to use -----
  #----------------------------------------------------------#

  if (
    is.null(age_uncertainty) &&
      isFALSE(standardise) &&
      isFALSE(is.null(rand)) &&
      (working_units == "levels" || bin_selection == "first")
  ) {
    if (isFALSE(silent) && isTRUE(verbose)) {
      util_output_comment(
        msg = paste(
          "There is no need for randomisation.",
          "Changing `rand` to NULL"
        )
      )
    }

    rand <- NULL
  }

  data_prepared <-
    prepare_data(
      data_source_prep = data_work,
      working_units = working_units,
      bin_size = bin_size,
      number_of_shifts = number_of_shifts,
      rand = rand
    )

  data_to_run <-
    util_flatten_list_by_one(data_prepared)

  #----------------------------------------------------------#
  # 4. Estimation -----
  #----------------------------------------------------------#

  if (isFALSE(silent) && isTRUE(verbose)) {
    util_output_heading(
      msg = "Start of estimation",
      size = "h2"
    )

    util_output_comment(
      msg = paste(
        "Number of estimation set to",
        length(data_to_run)
      )
    )
  }

  # select the preferred number of cores for of cores for parallel computation
  if (isTRUE(use_parallel)) {
    if (methods::is(use_parallel, "numeric")) {
      n_cores <-
        as.numeric(use_parallel) # set value
    } else {
      n_cores <-
        parallel::detectCores() # detect number
    }

    # create cluster
    cl <-
      parallel::makeCluster(n_cores)

    # eval packages
    parallel::clusterEvalQ(cl, {
      library("tidyverse")
      library("RRatepol")
    })
  } else {
    cl <- NULL
  }

  pbapply::pboptions(type = "timer")

  if (isTRUE(silent)) {
    pbapply_setting <-
      pbapply::pboptions(type = "none")
  }

  # run the estimation with progress bar
  result_list <-
    pbapply::pblapply(
      X = data_to_run,
      FUN = run_iteration,
      cl = cl,
      bin_selection = bin_selection,
      standardise = standardise,
      n_individuals = n_individuals,
      tranform_to_proportions = tranform_to_proportions,
      dissimilarity_coefficient = dissimilarity_coefficient,
      time_standardisation = time_standardisation,
      verbose = verbose,
      silent = silent
    )

  if (isTRUE(silent)) {
    pbapply::pboptions(pbapply_setting)
  }

  # close progress bar and cluster
  if (isFALSE(is.null(cl))) {
    parallel::stopCluster(cl)
    cl <- NULL
  }
  gc(verbose = FALSE)

  #----------------------------------------------------------#
  # 5. Results Summary -----
  #----------------------------------------------------------#

  # create new dataframe with summary of randomisation results
  results_full <-
    purrr::map_dfr(
      .x = result_list,
      .f = as.data.frame,
      .id = "it"
    ) %>%
    dplyr::group_by(.data$label) %>%
    dplyr::summarise(
      .groups = "drop",
      Age = stats::median(.data$res_age, na.rm = TRUE),
      ROC = stats::median(.data$roc, na.rm = TRUE),
      ROC_up = stats::quantile(.data$roc, 0.975, na.rm = TRUE),
      ROC_dw = stats::quantile(.data$roc, 0.025, na.rm = TRUE)
    ) %>%
    dplyr::select(
      Working_Unit = "label",
      "Age",
      "ROC",
      "ROC_up",
      "ROC_dw"
    )

  # reduce results by the focus age time
  if (methods::is(interest_threshold, "numeric")) {
    results_full <-
      results_full %>%
      dplyr::filter(
        .data$Age <= interest_threshold
      )
  }

  # final tibble (sort samples by age)
  results_full_fin <-
    results_full %>%
    dplyr::arrange(.data$Age)

  # time duration output
  end_time <- Sys.time()
  time_duration <- end_time - start_time

  if (isFALSE(silent)) {
    util_output_heading(
      paste(
        "RRatepol finished",
        end_time,
        "taking",
        round(time_duration, 2),
        units(time_duration)
      ),
      size = "h1"
    )
  }

  return(results_full_fin)
}
# end of code
