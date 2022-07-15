#' @title RRatepol: Estimate rate of change
#'
#' @param data_source_community
#' Data.frame. Community data with species as columns and
#' levels (samples) as rows. First column should be `sample_id` (character).
#' @param data_source_age
#' Data.frame with two columns:
#' \itemize{
#' \item `sample_id` - unique ID of each level (character)
#' \item `age` - age of level (numeric)
#' }
#' @param age_uncertainty
#' Usage of age uncertainty form Age-depth models. Either:
#' \itemize{
#' \item matrix with number of columns as number of samples. Each column is one sample,
#'  each row is one age sequence from age-depth model. Age sequence is randomly
#'  sampled from age-depth model uncertainties at the beginning of each run.
#'  \item `NULL` - Age uncertainties are not available and, therefore, will not be used.
#' }
#' @param smooth_method
#' Character. type of smoothing applied for the each of the pollen type
#' \itemize{
#' \item `"none"` - Pollen data is not smoothed
#' \item `"m.avg"` - Moving average
#' \item `"grim"` - Grimm's smoothing
#' \item `"age.w""` - Age-weighted average
#' \item `"shep"` - Shepard's 5-term filter
#' }
#' @param smooth_N_points
#' Numeric. Number of points for used for moving average,
#'  Grimm and Age-Weighted smoothing (odd number)
#' @param smooth_age_range
#' Numeric. Maximal age range for both Grimm and Age-weight smoothing
#' @param smooth_N_max
#' Numeric. Maximal number of samples to look in Grimm smoothing
#' @param Working_Units
#' Character. Selection of units that the DC will be calculated between.
#' \itemize{
#' \item `"levels"` - individual levels are going to be used
#' \item `"bins"` - samples in predefined bins will be pooled together and one sample
#' will be selected from each time bin as a representation.
#' \item `"MW"` - Bins of selected size are created, starting from the beginning of the core.
#'  This is repeated many times, with each time bin (window) shifting by Z years forward.
#'   This is repeated X times, where X = bin size / Z.
#' }
#' @param bin_size
#' Numeric. Size of the time bin (in years)
#' @param Number_of_shifts
#' Numeric. Value determining the number of shifts of window used
#' in Moving window method
#' @param bin_selection
#' Character. Setting determining the the process of selection of
#' samples from bins.
#' \itemize{
#' \item `"first"` - sample closest to the beginning of the bin is selected
#' as a representation.
#' \item `"random"` - a random sample is selected as a representation.
#' }
#' @param standardise
#' Logical. If `standardise` == `TRUE`, then standardise
#' each Working Unit to certain number of individuals (using random resampling
#' without repetition)
#' @param N_individuals
#' Numeric. Number of grain to perform standardisation to.
#' The `N_individual` is automatically adjusted to the smallest number
#' of pollen grains in sequence.
#' @param DC
#' Character. Dissimilarity coefficient. Type of calculation of differences
#' between Working Units
#' \itemize{
#' \item `"euc"` - Euclidean distance
#' \item `"euc.sd"` - Standardised Euclidean distance
#' \item `"chord"` - Chord distance
#' \item `"chisq"` - Chi-squared coefficient
#' \item `"gower"` - Gower's distance
#' \item `"bray"` - Bray-Curtis distance
#' }
#' @param tranform_to_proportions
#' Logical. Should the community data be transformed to a
#' proportion during calculations?
#' @param rand
#' Numeric. Number of runs used in randomisation.
#' @param treads
#' `r lifecycle::badge("deprecated")`
#' @param use_parallel
#' Preference of usage of parallel computation of randomisation
#' \itemize{
#' \item `[value]` - selected number of cores
#' \item `TRUE` - automatically selected number of cores
#' \item `FALSE` - does not use parallel computation (only single core)
#' }
#' @param interest_threshold
#' Numeric. Optional. Age, after which all results of RoC are excluded.
#' @param only_subsequent
#' `r lifecycle::badge("deprecated")`
#' Logical.
#' \itemize{
#' \item `FALSE` - RoC between WUs can be calculated using every consecutive WU
#' \item `TRUE` -  calculation of RoC can be restricted to only directly adjacent WUs
#' }
#' Using the former increases the number of samples for which RoC can be calculated
#' within a sequence, which varies in terms of sample resolution, but may still
#' introduce biases related to the RoC estimation as a result of the varying inter-sample distances.
#' Recommended setting is `only_subsequent` = `TRUE`.
#' Only `only_subsequent` = `TRUE` will be kept in the next version.
#' @param time_standardisation
#' Numeric. Units scaling for result RoC values. For example,
#' if `time_standardisation` = 100, the RoC will be reported as
#' dissimilarity per 100 yr.
#' @param verbose
#' Logical. If `TRUE`, function will output messages about internal processes
#' @param Debug
#' `r lifecycle::badge("deprecated")`
#' Use `verbose` instead.
#'
#' @description A function to estimate Rate of change in community data in time series
#' @details R-Ratepol is written as an R package and includes a range of
#' possible settings including a novel method to evaluate RoC in a single
#' stratigraphical sequence using assemblage data and age uncertainties for
#' each level. There are multiple built-in dissimilarity coefficients (DC) for
#' different types of assemblage data, and various levels of data smoothing
#' that can be applied depending on the type and variance of the data.
#' In addition, R-Ratepol can use randomisation, accompanied by use of age
#' uncertainties of each level and taxon standardisation to detect RoC patterns
#' in datasets with high data noise or variability (i.e. numerous rapid changes
#' in composition or sedimentation rates).
#'
#' The computation of RoC in R-Ratepol is performed using the following steps:
#' \enumerate{
#' \item Assemblage and age-model data are extracted from the original source and
#' should be compiled together, i.e. depth, age, variable (taxon) 1, variable (taxon) 2, etc.
#' \item (optional) Smoothing of assemblage data: Each variable within the
#' assemblage data is smoothed using one of five in-built smoothing methods:
#' \itemize{
#' \item none (`smooth_method` = `"none"`)
#' \item Shepard's 5-term filter (`smooth_method` = `"shep"`; Davis, 1986; Wilkinson, 2005)
#' \item moving average ´(`smooth_method` = `"m.avg"}`)
#' \item age-weighted average (`smooth_method` = `"age.w"`)
#' \item Grimm's smoothing (`smooth_method` = `"grim"`; Grimm and Jacobson, 1992)
#' }
#' \item Creation of time bins: A template for all time bins in all window movements is created.
#' \item A single run (an individual loop) is computed:
#' \itemize{
#' \item (optional) Selection of one time series from age uncertainties (see section on randomisation)
#' \item Subsetting levels in each bin: Here the working units (WU) are defined
#' \item (optional) Standardisation of assemblage data in each WU
#' \item The summary of a single run is produced based on all moving windows
#' \item Calculation of RoC between WUs: RoC is calculated as the dissimilarity
#' coefficient (DC) standardised by age differences between WUs. Five in-built
#' dissimilarity coefficients are available:
#' \itemize{
#' \item Euclidean distance (`DC` = `"euc"`)
#' \item standardised Euclidean distance (`DC` = `"euc.sd"`)
#' \item Chord distance (`DC` = `"chord"`)
#' \item Chi-squared coefficient (`DC` = `"chisq"`; Prentice, 1980)
#' \item Gower's distance (`DC` = `"gower"`;Gower, 1971)
#' \item Bray-Curtis distance (`DC` = `"bray"`)
#' }
#' The choice of DC depends on the type of assemblage data. In addition, RoC
#' between WUs be calculated using every consecutive WU (`only_subsequent` = `FALSE`),
#' or alternatively, calculation of RoC can be restricted to only directly
#' adjacent WUs (`only_subsequent` = `TRUE`). Using the former increases
#' the number of samples for which RoC can be calculated within a sequence,
#' which varies in terms of sample resolution, but may still introduce
#' biases related to the RoC estimation as a result of the varying
#' inter-sample distances.
#' }
#' \item Step 4 is repeated multiple times (e.g. 10,000 times).
#' \item Validation and summary of results from all runs of RoC calculation are produced.
#' \item (Optional) Data beyond a certain age can be excluded.
#' }
#' ## Selection of working units (WU; Step 3)
#' RoC is calculated between consecutive Working Units (WU). Traditionally,
#' these WUs represent individual stratigraphical levels. However, changes in
#' sedimentation rates and sampling strategies can result in an uneven temporal
#' distribution of levels within a time sequence, which in turn makes
#' the comparison of RoC between sequences problematic. There are various methods
#' that attempt to minimise such problems. The first is interpolation of levels
#' to evenly spaced time intervals, and the use of the interpolated data as WUs.
#' This can lead to a loss of information when the density of levels is high.
#' Second is binning of levels: assemblage data are pooled into age brackets
#' of various size (i.e. time bins) and these serve as WUs. Here, the issue
#' is a lower resolution of WUs and their uneven size in terms of total
#' assemblage count (bins with more levels have higher assemblage counts).
#' Third is selective binning: like classical binning, bins of selected size
#' are created, but instead of pooling assemblage data together, only one
#' level per time bin is selected as representative of each bin. This results
#' in an even number of WUs in bins with a similar count size in the assemblage.
#' However, the issue of low resolution remains.
#' Therefore, we propose a new method of binning with a moving window,
#' which is a compromise between using individual levels and selective binning.
#' This method follows a simple sequence: time bins are created,
#' levels are selected as in selective binning, and RoC between bins is calculated.
#' However, the brackets of the time bin (window) are then moved forward by a
#' selected amount of time (Z), levels are selected again (subset into bins),
#' and RoC calculated for the new set of WUs. This is repeated X times
#' (where X is the bin size divided by Z) while retaining all the results.
#'
#' R-Ratepol currently provides several options for selecting WU, namely as i
#' ndividual levels (`Working_Units` = `"levels"`), selective binning of levels
#' (`Working_Units` = `"bins"`), and our new method of binning with a moving
#' window (`Working_Units` = `"MW"`)
#'
#' ## Randomisation
#' Due to the inherent statistical errors in uncertainties in the age estimates
#' from age-depth and the assemblage datasets (e.g. pollen counts in each level;
#' Birks and Gordon, 1985), R-Ratepol can be run several times and the results
#' summarised (Steps 5-6). Therefore, two optional settings are available by
#' using age uncertainties and assemblage data standardisation.
#'
#' ## Age uncertainties
#' For each run, a single age sequence from the age uncertainties is randomly
#' selected. The calculation between two consecutive WUs (i.e. one working-unit
#' combination) results in a RoC score and a time position (which is calculated
#' as the mean age position of the two WUs). However, due to random sampling
#' of the age sequence, each WU combination will result in multiple RoC values.
#' The final RoC value for a single WU combination is calculated as the median
#' of the scores from all randomisations. In addition, the 95th quantile from all
#' randomisations is calculated as an error estimate.

#' ## Data standardisation (Step 4b)
#' Taxa in the assemblage dataset can be standardised to a certain count
#' (e.g. number of pollen grains in each WU) by rarefaction. Random sampling
#' without replacement is used to draw a selected number of individuals from
#' each WU (e.g. 150 pollen grains).

#' @references
#' Birks, H.J.B., Gordon, A.D., 1985. Numerical Methods in Quaternary Pollen
#' Analysis. Academic Press, London.
#'
#' Davis, J.C., 1986. Statistics and Data Analysis in Geology, 2nd edn. ed.
#' J. Wiley & Sons, New York.
#'
#' Gower, J.C., 1971. A general coefficient of similarity and some of its
#' properties. Biometrics 27, 857–871.
#'
#' Grimm, E.C., Jacobson, G.L., 1992. Fossil-pollen evidence for abrupt
#' climate changes during the past 18000 years in eastern North America.
#' Clim. Dyn. 6, 179–184.
#'
#' Prentice, I.C., 1980. Multidimensional scaling as a research tool in
#' Quaternary palynology: A review of theory and methods. Rev. Palaeobot.
#' Palynol. 31, 71–104. https://doi.org/10.1016/0034-6667(80)90023-8
#'
#' Wilkinson, L., 2005. The Grammar of Graphics. Springer-Verlag, New York,
#' USA 37. https://doi.org/10.2307/2669493
#' @export
#' @examples
#' \dontrun{
#' example_data <- RRatepol::example_data
#'
#' sequence_01 <-
#'   fc_estimate_RoC(
#'     data_source_community = example_data$pollen_data[[1]],
#'     data_source_age = example_data$sample_age[[1]],
#'     age_uncertainty = FALSE,
#'     smooth_method = "shep",
#'     Working_Units = "MW",
#'     rand = 1e3,
#'     use_parallel = TRUE,
#'     DC = "chisq"
#'   )
#'
#' fc_plot_RoC_sequence(
#'   sequence_01,
#'   age_threshold = 8e3,
#'   Roc_threshold = 1
#' )
#' }
fc_estimate_RoC <-
  function(data_source_community,
           data_source_age,
           age_uncertainty = NULL,
           smooth_method = c("none", "m.avg", "grim", "age.w", "shep"),
           smooth_N_points = 5,
           smooth_age_range = 500,
           smooth_N_max = 9,
           Working_Units = c("levels", "bins", "MW"),
           bin_size = 500,
           Number_of_shifts = 5,
           bin_selection = c("random", "first"),
           standardise = FALSE,
           N_individuals = 150,
           DC = c("euc", "euc.sd", "chord", "chisq", "gower", "bray"),
           tranform_to_proportions = TRUE,
           rand = NULL,
           treads = NULL,
           use_parallel = TRUE,
           interest_threshold = NULL,
           only_subsequent = TRUE,
           time_standardisation = bin_size,
           verbose = FALSE,
           Debug = NULL) {

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

    util_check_class("data_source_community", "data.frame")

    util_check_class("data_source_age", "data.frame")

    util_check_class("age_uncertainty", c("NULL", "matrix"))

    util_check_class("Working_Units", "character")

    util_check_vector_values("Working_Units", c("levels", "bins", "MW"))

    Working_Units <- match.arg(Working_Units)

    util_check_class("time_standardisation", "numeric")

    util_check_if_integer("time_standardisation")

    if (
      Working_Units != "levels"
    ) {
      util_check_class("bin_size", "numeric")

      util_check_if_integer("bin_size")

      util_check_class("bin_selection", "character")

      util_check_vector_values("bin_selection", c("first", "random"))

      util_check_class("only_subsequent", "logical")

      if (
        only_subsequent == FALSE
      ) {
        lifecycle::deprecate_stop("1.0.0", "fc_estimate_RoC(only_subsequent)")
      }

      if (
        Working_Units == "MW"
      ) {
        util_check_class("Number_of_shifts", "numeric")

        util_check_if_integer("Number_of_shifts")
      }
    }

    util_check_class("standardise", "logical")

    if (
      standardise == TRUE
    ) {
      util_check_class("N_individuals", "numeric")

      util_check_if_integer("N_individuals")
    }

    util_check_class("tranform_to_proportions", "logical")

    util_check_class("interest_threshold", c("NULL", "numeric"))

    util_check_class("smooth_method", "character")

    util_check_vector_values(
      "smooth_method",
      c("none", "m.avg", "grim", "age.w", "shep")
    )

    smooth_method <- match.arg(smooth_method)

    if (
      !smooth_method %in% c("none", "shep")
    ) {
      assertthat::assert_that(
        smooth_N_points %% 2 != 0,
        msg = "'smooth_N_points' must be an odd number"
      )

      if (
        smooth_method != "m.avg"
      ) {
        util_check_class("smooth_age_range", "numeric")

        if (
          smooth_method == "grim"
        ) {
          assertthat::assert_that(
            smooth_N_max %% 2 != 0,
            msg = "'smooth_N_max' must be an odd number"
          )

          assertthat::assert_that(
            smooth_N_points < smooth_N_max,
            msg = "'smooth_N_max' must be bigger than 'smooth_N_points"
          )
        }
      }
    }

    util_check_class("DC", "character")

    util_check_vector_values(
      "DC",
      c("euc", "euc.sd", "chord", "chisq", "gower", "bray")
    )

    DC <- match.arg(DC)

    util_check_class("rand", c("NULL", "numeric"))

    if (
      is.null(rand) == FALSE
    ) {
      util_check_if_integer("rand")
    }

    if (
      is.null(treads) == FALSE
    ) {
      lifecycle::deprecate_warn(
        "1.0.0",
        "fc_estimate_RoC(treads)",
        "fc_estimate_RoC(use_parallel)"
      )

      use_parallel <- treads
    }
    util_check_class("use_parallel", c("logical", "numeric"))

    if (
      is.numeric(use_parallel)
    ) {
      util_check_if_integer("use_parallel")
    }

    if (
      is.null(Debug) == FALSE
    ) {
      lifecycle::deprecate_warn(
        "1.0.0",
        "fc_estimate_RoC(Debug)",
        "fc_estimate_RoC(verbose)"
      )

      verbose <- Debug
    }

    util_check_class("verbose", "logical")

    #--------------------------------------------------#
    # 0.1. Report to user -----
    #--------------------------------------------------#

    start_time <- Sys.time()

    util_output_heading(
      paste("RRatepol started", start_time),
      size = "h1"
    )

    if (
      is.null(age_uncertainty == FALSE)
    ) {
      util_output_comment(
        "'age_uncertainty' will be used for in the RoC estimation"
      )

      if (rand < 100) {
        util_output_warning(
          paste(
            "'age_uncertainty' was selected to be used with low number",
            "of replication. Recommend to increase 'rand'"
          )
        )
      }
    }

    switch(Working_Units,
      "levels" = {
        util_output_comment(
          "RoC will be estimated between individual subsequent levels"
        )
      },
      "bins" = {
        util_output_comment(
          paste(
            "RoC will be estimated using selective binning with", bin_size,
            "yr time bin"
          )
        )
      },
      "MW" = {
        util_output_comment(
          paste(
            "RoC will be estimated using 'binning with the mowing window' of",
            bin_size, "yr time bin over", Number_of_shifts, "number of window shifts"
          ),
          "\n",
          fill = TRUE
        )
      }
    )

    if (
      Working_Units != "levels"
    ) {
      if (
        bin_selection == "random"
      ) {
        util_output_comment(
          "Sample will randomly selected for each bin",
          "\n",
          fill = TRUE
        )

        if (
          rand < 100
        ) {
          util_output_warning(
            paste(
              "'bin_selection' was selected as 'random' with low number",
              "of replication. Recommend to increase 'rand'"
            )
          )
        }
      } else {
        util_output_comment(
          "First sample of each time bin will selected"
        )
      }
    }

    util_output_comment(
      paste(
        "'time_standardisation' =", time_standardisation, ":",
        "RoC values will be reported as disimilarity per", time_standardisation,
        "years."
      )
    )

    if (
      Working_Units != "levels" && time_standardisation != bin_size
    ) {
      util_output_comment(
        paste(
          "RoC values will be reported in different units than size of bin.",
          "Recommend to keep 'time_standardisation'",
          "and 'bin_size' as same values"
        )
      )
    }

    if (
      standardise == TRUE
    ) {
      util_output_comment(
        paste(
          "Data will be standardise in each Working unit to", N_individuals,
          "or the lowest number detected in dataset"
        )
      )

      if (rand < 100) {
        util_output_warning(
          paste(
            "'standardise' was selected as 'TRUE' with low number of replication.",
            "Recommend to increase 'rand'"
          )
        )
      }
    }


    #----------------------------------------------------------#
    # 1. Data extraction -----
    #----------------------------------------------------------#

    # extract data into working format
    # already include data check
    data_extract <-
      fc_extract_data(
        data_community_extract = data_source_community,
        data_age_extract = data_source_age,
        age_uncertainty = age_uncertainty,
        verbose = verbose
      )


    #----------------------------------------------------------#
    # 2. Data smoothing -----
    #----------------------------------------------------------#

    if (
      smooth_method != "none"
    ) {
      # smooth data by selected smoothing type
      data_smooth <-
        fc_smooth_community_data(
          data_source_smooth = data_extract,
          smooth_method = smooth_method,
          smooth_N_points = smooth_N_points,
          smooth_N_max = smooth_N_max,
          smooth_age_range = smooth_age_range,
          round_results = standardise,
          verbose = verbose
        )
    } else {
      data_smooth <- data_extract
    }

    # reduce data dimentions
    data_work <-
      fc_reduce(
        data_smooth
      )

    if (
      verbose == TRUE
    ) {
      fc_check_data(data_work)
    }


    #----------------------------------------------------------#
    # 3. Crete datasets to use -----
    #----------------------------------------------------------#

    data_to_run <-
      fc_prepare_data(
        data_source_prep = data_work,
        Working_Units = Working_Units,
        bin_size = bin_size,
        Number_of_shifts = Number_of_shifts,
        rand = rand
      )

    #----------------------------------------------------------#
    # 4. Estimation -----
    #----------------------------------------------------------#

    # select the prefetred number of cores for of cores for parallel computation
    if (
      use_parallel == FALSE
    ) {
      result_table <-
        lapply(
          X = data_to_run,
          FUN = ~ fc_run_iteration()
        )
    } else {
      if (
        class(use_parallel) == "numeric"
      ) {
        n_cores <-
          as.numeric(use_parallel) # set value
      } else {
        n_cores <-
          parallel::detectCores() # detect number
      }

      # create cluster
      cl <- parallel::makeCluster(n_cores)
      parallel::clusterEvalQ(cl, {
        library("RRatepol")
        library("tidyverse")
      })

      result_table <-
        parallel::parLapply(
          X = data_to_run,
          fun = ~ fc_run_iteration(),
          cl = cl
        )

      # close progress bar and cluster
      if (!is.null(cl)) {
        parallel::stopCluster(cl)
        cl <- c()
      }
      gc()
    }

    #----------------------------------------------------------#
    # 5. Results Summary -----
    #----------------------------------------------------------#

    # create new dataframe with summary of randomisation results

    # extract results and match them by bin
    results_full <-
      dplyr::right_join(
        fc_extract_result(
          result_table,
          "RoC",
          rand
        ),
        fc_extract_result(
          result_table,
          "age_position",
          rand
        ),
        by = c("sample_id", "shift", "age_distance")
      )

    # reduce results by the focus age time
    if (interest_threshold != FALSE) {
      results_full <-
        dplyr::filter(
          results_full,
          age_position <= interest_threshold
        )
    }

    # final tibble (sort samples by age and select variables)
    results_full_fin <-
      results_full %>%
      dplyr::arrange(age_position) %>%
      dplyr::select(
        sample_id,
        age_position,
        RoC,
        RoC_95q,
        RoC_05q
      )

    names(results_full_fin) <- c("Working_Unit", "Age", "ROC", "ROC_up", "ROC_dw")

    end_time <- Sys.time()
    time_duration <- end_time - start_time
    cat(paste(
      "R-RATEPOL finished", end_time, "taking", time_duration, units(time_duration)
    ),
    fill = TRUE
    )

    return(results_full_fin)
  }
# end of code
