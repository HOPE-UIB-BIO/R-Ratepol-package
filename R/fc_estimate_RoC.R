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
#'  \item `FALSE` - Age uncertainties are not available and, therefore, will not be used.
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
#' \itemize {
#' \item `"first"` - sample closest to the beginning of the bin is selected 
#' as a representation.
#' \item `"random"` - a random sample is selected as a representation.
#' }
#' @param rand 
#' Numeric. Number of runs used in randomisation.
#' @param treads 
#' Preference of usage of parallel computation of randomisation
#' \itemize{
#' \item `[value]` - selected number of cores
#' \item `TRUE` - automatically selected number of cores
#' \item `FALSE` - does not use parallel computation (only single core)
#' }
#' @param standardise 
#' Logical. If `standardise` == `TRUE`, then standardise 
#' each Working Unit to certain number of individuals (using random resampling 
#' without repetition)
#' @param N_individuals 
#' Numeric. Number of grain to perform standardisation to. 
#' The `N_individual` is automatically adjusted to the smallest number 
#' of pollen grains in sequence. 
#' @param tranform_to_proportions 
#' Logical. Should the community data be transformed to a 
#' proportion during calculations?
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
#' @param interest_threshold 
#' Numeric. Age, after which all results of RoC are excluded 
#' before detection of peak points 
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
#' dissimilarity per 100 yr. If `time_standardisation` = `"auto"` (not recommended), 
#' RoC values will be reported as  standardised by the average distance between 
#' Working Units (levels/bins)" 
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
fc_estimate_RoC <- 
  function(
    data_source_community,
    data_source_age,
    age_uncertainty = FALSE,
    smooth_method = c("none", "m.avg", "grim", "age.w", "shep"),
    smooth_N_points = 5,
    smooth_age_range = 500,
    smooth_N_max = 9,
    Working_Units = c("levels", "bins", "MW"),
    bin_size = 500,
    Number_of_shifts = 5,
    bin_selection = c( "first", "random"),
    rand = 1,
    treads = TRUE,
    standardise = FALSE,
    N_individuals = 150,
    tranform_to_proportions = TRUE,
    DC = c("euc", "euc.sd", "chord", "chisq", "gower", "bray"),
    interest_threshold = FALSE,
    only_subsequent = TRUE,
    time_standardisation = 100,
    verbose = FALSE,
    Debug = NULL){
    
    # Start of the code
    
    #----------------------------------------------------------# 
    # 0. Arguments check -----
    #----------------------------------------------------------#
    
    assertthat::assert_that(
      !missing(data_source_community),
      msg = "Object 'data_source_community' must be included as a 'data.frame'")
    
    assertthat::assert_that(
      !missing(data_source_age),
      msg = "Object 'data_source_age' must be included as a 'data.frame'")
    
    if(missing(age_uncertainty)){
      age_uncertainty <- FALSE
    }
    
    assertthat::assert_that(
      is.numeric(rand),
      msg = "'rand' must be a 'numeric'")
    
    assertthat::assert_that(
      round(rand) == rand,
      msg = "'rand' must be a whole number")
    
    Working_Units <- match.arg(Working_Units)
    
    assertthat::assert_that(
      Working_Units %in% c("levels", "bins", "MW"),
      msg = "'Working_Units' must be a 'levels' 'bins' or 'MW'")
    
    assertthat::assert_that(
      is.numeric(time_standardisation) | time_standardisation == "auto",
      msg = "'time_standardisation' must be a 'numeric' or 'auto'")
    
    if( is.numeric(time_standardisation)){
      assertthat::assert_that(
        round(time_standardisation) == time_standardisation,
        msg = "'time_standardisation' must be a whole number")  
    }
    
    if(Working_Units != "levels"){
      
      assertthat::assert_that(
        is.numeric(bin_size),
        msg = "'bin_size' must be a 'numeric'")
      
      assertthat::assert_that(
        round(bin_size) == bin_size,
        msg = "'bin_size' must be a whole number")
      
      bin_selection <- match.arg(bin_selection)
      
      assertthat::assert_that(
        bin_selection == "first" | bin_selection == "random",
        msg = "'bin_selection' must be a 'first' or 'random'")
      
      assertthat::assert_that(
        is.logical(only_subsequent),
        msg = "'only_subsequent' must be a 'TRUE' or 'FALSE'")
      
      if (Working_Units == "MW"){
        assertthat::assert_that(
          is.numeric(Number_of_shifts),
          msg = "'Number_of_shifts' must be a 'numeric'")
        
        assertthat::assert_that(
          round(Number_of_shifts) == Number_of_shifts,
          msg = "'Number_of_shifts' must be a whole number")
      }
    }  
    
    assertthat::assert_that(
      is.logical(standardise),
      msg = "'standardise' must be a 'TRUE' or 'FALSE'")
    
    if(standardise == TRUE){
      assertthat::assert_that(
        is.numeric(N_individuals),
        msg = "'N_individuals' must be a 'numeric'")
      
      assertthat::assert_that(
        round(N_individuals) == N_individuals,
        msg = "'N_individuals' must be a whole number")
    }
    
    assertthat::assert_that(
      is.logical(tranform_to_proportions),
      msg = "'tranform_to_proportions' must be a 'TRUE' or 'FALSE'")
    
    assertthat::assert_that(
      is.numeric(interest_threshold) | interest_threshold == FALSE,
      msg = "'interest_threshold' must be a 'numeric' or 'FALSE'")
    
    smooth_method <-  match.arg(smooth_method)
    
    assertthat::assert_that(
      any(smooth_method == c("none", "m.avg", "grim", "age.w", "shep")),
      msg = "'smooth_method' must be one of the following:
    'none', 'm.avg', 'grim', 'age.w', 'shep'")
    
    if(!smooth_method %in% c("none", "shep")){
      
      assertthat::assert_that(
        smooth_N_points%%2 != 0,
        msg = "'smooth_N_points' must be an odd number")
      
      if(smooth_method != "m.avg"){
        assertthat::assert_that(
          is.numeric(smooth_age_range),
          msg = "'smooth_age_range' must be 'numeric")
        
        if(smooth_method == "grim"){
          assertthat::assert_that(
            smooth_N_max%%2 != 0,
            msg = "'smooth_N_max' must be an odd number")
          
          assertthat::assert_that(
            smooth_N_points < smooth_N_max,
            msg = "'smooth_N_max' must be bigger than 'smooth_N_points")
        }
      }
    }
    
    DC <- match.arg(DC)
    
    assertthat::assert_that(
      any(DC == c("euc", "euc.sd", "chord", "chisq", "gower", "bray")),
      msg = "'DC' must be one of the following:
    'euc', 'euc.sd', 'chord', 'chisq', 'gower', 'bray'")
    
    assertthat::assert_that(
      is.numeric(treads) | is.logical(treads),
      msg = "'treads' must be a 'numeric' or 'TRUE'/'FALSE'")
    
    if(is.numeric(treads)){
      assertthat::assert_that(
        round(treads) == treads,
        msg = "'treads' must be a whole number")
    }
    
    if(is.null(Debug) == FALSE){
      lifecycle::deprecate_warn(
        "0.0.7",
        "fc_estimate_RoC(Debug)",
        "fc_estimate_RoC(verbose)")
      
      verbose <- Debug
      
    }
    
    assertthat::assert_that(
      is.logical(verbose) == TRUE,
      msg = "'verbose' must be a logical"
    )
    
    #--------------------------------------------------#
    # 0.1. Report to user -----
    #--------------------------------------------------#
    cat("\n")
    
    start_time <- Sys.time()
    cat(paste("R-RATEPOL started", start_time),
        "\n", fill = TRUE)
    
    if (!all(age_uncertainty == FALSE)){
      cat(
        "'age_uncertainty' will be used for in the RoC estimation",
        "\n", fill = TRUE)
      
      if (rand < 100){
        cat(
          paste(
            "'age_uncertainty' was selected to be used with low number",
            "of replication. Recommend to increase 'rand'"),
          "\n", fill = TRUE)
      }
    }
    
    if(smooth_method == "m.avg"){
      cat(
        paste(
          "Data will be smoothed by 'moving average' over", smooth_N_points,
          "points"),
        "\n",fill = TRUE)  
    }
    
    if(smooth_method == "grim"){
      cat(
        paste(
          "Data will be smoothed by 'Grimm method' with min samples", smooth_N_points,
          "max samples", smooth_N_max, "and max age range of", smooth_age_range),
        "\n", fill = TRUE)  
    }
    
    if(smooth_method == "age.w"){
      cat(
        paste(
          "Data will be smoothed by 'age-weighed average' over", smooth_N_points,
          "points with a threshold of", smooth_age_range),
        "\n", fill = TRUE)  
    }
    
    if(smooth_method == "shep"){
      cat(
        paste(
          "Data will be smoothed by 'Shepard's 5-term filter'"),
        "\n", fill = TRUE)  
    }
    
    if(Working_Units == "levels"){
      cat(
        "RoC will be estimated between individual subsequent levels",
        "\n", fill = TRUE)  
    }
    
    if(Working_Units != "levels"){
      
      if(Working_Units == "bins"){
        cat(
          paste(
            "RoC will be estimated using selective binning with", bin_size,
            "yr time bin"),
          "\n", fill = TRUE)  
      }
      
      if(Working_Units == "MW"){
        cat(
          paste(
            "RoC will be estimated using 'binning with the mowing window' of",
            bin_size, "yr time bin over", Number_of_shifts, "number of window shifts"),
          "\n", fill = TRUE)  
      }
      
      if(bin_selection == "random"){
        cat(
          "Sample will randomly selected for each bin",
          "\n", fill = TRUE)  
        
        if(rand < 100){
          cat(
            paste(
              "'bin_selection' was selected as 'random' with low number",
              "of replication. Recommend to increase 'rand'"),
            "\n", fill = TRUE)  
        }
        
      } else {
        cat(
          "First sample of each time bin will selected",
          "\n", fill = TRUE)  
      }
      
      if(only_subsequent == FALSE & Working_Units == "MW"){
        cat(
          paste(
            "'only_subsequent == FALSE' and 'Working_Units == MW'.",
            "This is not a recommended setting.",
            "Please use 'only_subsequent == TRUE' for 'Working_Units == MW'",
            "see '?fc_estimate_ROC' for more information"), 
          "\n", fill = TRUE)
      } else if(only_subsequent == FALSE){
        cat(
          paste(
            "'only_subsequent' was selected as 'FALSE'.",
            "This is not a recommended setting. Results will be affected",
            "see '?fc_estimate_ROC' for more information"),
          "\n", fill = TRUE)
      }
      
    }
    
    if(is.numeric(time_standardisation)){
      cat(
        paste(
          "'time_standardisation' =", time_standardisation,":",
          "RoC values will be reported as disimilarity per", time_standardisation,
          "years."),
        "\n", fill = TRUE)
      
      if(Working_Units != "levels" & time_standardisation != bin_size){
        cat(
          paste(
            "RoC values will be reported in different units than size of bin.",
            "Recommend to keep 'time_standardisation'",
            "and 'bin_size' as same values"),
          "\n", fill = TRUE)
      }
    }
    
    if(time_standardisation == "auto"){
      cat(
        paste(
          "'time_standardisation' = 'auto' is not recomended setting.", 
          "RoC values will be reported as standardised by the average distance", 
          "between Working Units (levels/ bins)"),
        "\n", fill = TRUE)
    }
    
    if(standardise == TRUE ){
      cat(
        paste("Data will be standardise in each Working unit to", N_individuals,
              "or the lowest number detected in dataset"),
        "\n", fill=TRUE)
      
      if(rand < 100){
        cat(
          paste(
            "'standardise' was selected as 'TRUE' with low number of replication.",
            "Recommend to increase 'rand'"),
          "\n", fill = TRUE) 
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
        verbose = verbose)
    
    
    #----------------------------------------------------------#
    # 2. Data smoothing ----- 
    #----------------------------------------------------------#
    data_smooth <- data_extract
    
    # smoothdata by selected smoothing type
    data_smooth <- 
      fc_smooth_community_data(
        data_smooth,
        smooth_method = smooth_method,
        smooth_N_points = smooth_N_points,
        smooth_N_max = smooth_N_max,
        smooth_age_range = smooth_age_range,
        round_results = standardise,
        verbose = verbose)
    
    #check data and reduce data dimentions
    data_work <- 
      fc_check_data(
        data_smooth, 
        proportion = FALSE,
        verbose = verbose)
    
    
    #----------------------------------------------------------#
    # 3. Working Unit selection -----
    #----------------------------------------------------------#
    
    if(Working_Units == "levels"){
      Number_of_shifts <-  1
    } else if(Working_Units == "bins"){
      Number_of_shifts <-  1
      
      bin_sizes <- 
        fc_create_bins(
          data_work,
          shift_value = bin_size,
          Number_of_shifts = 1)
    } else if(Working_Units == "MW"){
      
      shift_value <-  bin_size/Number_of_shifts
      
      bin_sizes <- 
        fc_create_bins(
          data_work,
          shift_value = shift_value,
          Number_of_shifts = Number_of_shifts)
    }
    
    
    #----------------------------------------------------------#
    # 4. Randomisation ----- 
    #----------------------------------------------------------#
    
    # create template for result tibble
    shift_tibble_template <-  tibble::tibble()
    
    # select the prefetred number of cores for of cores for parallel computation
    if(class(treads) == "numeric"){
      Ncores <-  treads # set value
    } else {
      if (treads == TRUE) {
        Ncores <-  parallel::detectCores() # detect number
      } else {
        Ncores <-  1
      }
    }
    
    if(rand > 1){
      cat(
        paste(
          "Starting the randomisation. Number of randomisations set to", rand),
        "\n", fill = TRUE)  
    }
    
    doFuture::registerDoFuture()
    future::plan(multisession, workers = Ncores)
    progressr::handlers(global = TRUE)
    progressr::handlers("progress")
    
    `%dorng%` <- doRNG::`%dorng%`# so %dorng% doesn't need to be attached
    
    # progress bar
    p <- progressr::progressor(steps = rand)
    
    result_tibble <-  
      foreach::foreach(
        l = 1:rand,
        .export = c(
          "data_work",
          "N_individuals"),
        .combine = rbind) %dorng% {
          
          # TIME SAMPLING
          # sample random time sequence from time uncern.
          data_work@Age$newage <- 
            as.numeric(data_work@Age.un[sample(c(1:max(1, nrow(data_work@Age.un))), 1), ])
          
          # create result tible
          shift_tibble <-  shift_tibble_template
          
          # repeat for number of shifts
          for(k in 1:Number_of_shifts){
            
            #----------------------------------------------------------#
            # 4.1 Data subsetting ----- 
            #----------------------------------------------------------#
            data_subset <-  data_work
            
            # select one sample for each bin based on the age of the samples. 
            # Sample is choses if it is the closes one to the upper end of the bin
            if (Working_Units != "levels"){
              
              # select bin for this shift
              selected_bins <-  bin_sizes[bin_sizes$shift == k, ]
              
              #subset data
              data_subset <-  
                fc_subset_samples(
                  data_subset,
                  bins = selected_bins,
                  WU = Working_Units,
                  bin_selection = bin_selection)
              
              data_subset <-  
                fc_check_data(
                  data_subset,
                  proportion = FALSE)
            }
            
            #----------------------------------------------------------#
            # 4.2 Data Standardisation ----- 
            #----------------------------------------------------------#
            data_sd <-  data_subset
            
            # standardisation of community data to X(N_individuals) number of individuals
            if(standardise == TRUE){
              
              # adjust the value by the minimal Community or to a minimal of presected values
              N_individuals <-  min(c(rowSums(data_sd@Community), N_individuals) )
              
              # check if all samples has N_individuals of individuals
              data_sd@Age <-  
                data_sd@Age[rowSums(data_sd@Community, na.rm = TRUE) >= N_individuals, ]
              
              data_sd@Age.un <-  
                data_sd@Age.un[ ,rowSums(data_sd@Community, na.rm = TRUE) >= N_individuals]
              
              data_sd@Community <-  
                data_sd@Community[rowSums(data_sd@Community, na.rm = TRUE) >= N_individuals, , drop = FALSE]
              
              data_sd <-  
                fc_check_data(
                  data_sd,
                  proportion = FALSE,
                  Samples = TRUE,
                  verbose = verbose)
              
              # standardisation
              data_sd <-  
                fc_standardise_community_data(
                  data_source = data_sd, 
                  N_individuals = N_individuals,
                  verbose = verbose)
              
              assertthat::assert_that(
                any(rowSums(data_sd@Community, na.rm = TRUE) == N_individuals),
                msg = "Data standardisation was unsuccesfull, try 'standardise' = FALSE")
            }
            
            # data check with proportioning
            data_sd_check <-  
              fc_check_data(
                data_source_check = data_sd,
                proportion = tranform_to_proportions,
                Samples = FALSE,
                verbose = verbose)
            
            
            #----------------------------------------------------------#
            # 4.3 DC Calculation ----- 
            #----------------------------------------------------------#
            
            # calculate DC between each subsequent samples/bins
            DC_res <-  
              fc_calculate_DC(
                data_source_DC = data_sd_check,
                DC = DC,
                verbose = verbose)
            
            
            #----------------------------------------------------------#
            # 4.4 Age Standardisation ----- 
            #----------------------------------------------------------#
            
            # create empty tible with size = number of samples-1
            shift_tibble_res <- 
              tibble::tibble(
                DC = DC_res) 
            
            # create empty vectors for age difference calcualtion
            shift_tibble_res$age_diff <-  
              vector(
                mode = "numeric",
                length = nrow(shift_tibble_res))
            
            shift_tibble_res$bin <-  
              vector(
                mode = "character",
                length = nrow(shift_tibble_res))
            
            shift_tibble_res$age_distance <- 
              vector(
                mode = "numeric",
                length = nrow(shift_tibble_res))
            
            shift_tibble_res$age_position <-  
              vector(
                mode = "numeric",
                length = nrow(shift_tibble_res))
            
            for (i in 1:nrow(shift_tibble_res)){ # for each RoC
              
              # calcualte the age difference between subsequesnt samples
              shift_tibble_res$age_diff[i] <-  
                data_sd_check@Age$newage[i + 1] - data_sd_check@Age$newage[i]
              
              # Set age difference as 1, if age difference between samples is 
              #   smaller than 1
              if(shift_tibble_res$age_diff[i] < 1){shift_tibble_res$age_diff[i] <-  1}
              
              #calculate the average position of RoC
              shift_tibble_res$age_position[i] <-  
                mean(c(data_sd_check@Age$age[i+1],
                       data_sd_check@Age$age[i]))
              
              # create vector with WU names
              shift_tibble_res$bin[i] <-
                paste(
                  row.names(data_sd_check@Age)[i],
                  "-",
                  row.names(data_sd_check@Age)[i + 1])
              
              if(Working_Units != "levels"){
                shift_tibble_res$age_distance[i] <- 
                  as.numeric(row.names(data_sd_check@Age)[i + 1]) -
                  as.numeric(row.names(data_sd_check@Age)[i])
              } else {
                shift_tibble_res$age_distance[i] <- NA
              }
            }
            
            # remove the non-subsequent levels.
            if(Working_Units != "levels" & only_subsequent == TRUE){
              shift_tibble_res <-
                shift_tibble_res %>% 
                dplyr::filter(age_distance <= bin_size)
            }
            
            if(time_standardisation == "auto"){
              if(Working_Units != "levels"){
                time_standardisation_unit <- bin_size
              } else {
                time_standardisation_unit <- mean(shift_tibble_res$age_diff)
              }
            } else {
              time_standardisation_unit <- time_standardisation
            }
            
            if (verbose == TRUE){
              cat("", fill = TRUE)
              cat(paste("The time standardisation unit (TSU) is",
                        round(time_standardisation_unit,2)), fill=TRUE)
            }
            
            #  calculate DC standardise by time
            shift_tibble_res <-
              shift_tibble_res %>% 
              dplyr::mutate(
                age_diff_st = age_diff / time_standardisation_unit,
                RoC = DC / age_diff_st,
                shift = k)
            
            
            #----------------------------------------------------------#
            # 4.5 Result of a single window shift -----
            #----------------------------------------------------------#
            
            # add the results from this shift into the result tibble
            shift_tibble <- 
              rbind(
                shift_tibble,
                shift_tibble_res)
          }
          
          
          #----------------------------------------------------------#
          # 4.6 Result of a single randomisation run -----
          #----------------------------------------------------------#
          
          if(nrow(shift_tibble)<1 & Working_Units != "levels" & only_subsequent == TRUE){
            stop("Estimation not succesfull, try increase the bin size")
          }
          
          shift_tibble <- 
            shift_tibble %>% 
            dplyr::mutate(
              ID = l)
          
          p(sprintf("l=%g", l))
          
          return(shift_tibble)
        }# end of the randomization
    
    # close progress bar and cluster
    # parallel::stopCluster(cl)
    
    
    #----------------------------------------------------------#
    # 5. Results Summary -----
    #----------------------------------------------------------#
    
    # create new dataframe with summary of randomisation results
    
    # extract results and match them by bin
    results_full <- 
      dplyr::right_join(
        fc_extract_result(
          result_tibble,
          "RoC",
          rand),
        fc_extract_result(
          result_tibble,
          "age_position",
          rand),
        by = c("sample_id","shift","age_distance"))
    
    # reduce results by the focus age time
    if(interest_threshold != FALSE){
      results_full <- 
        dplyr::filter(
          results_full,
          age_position <= interest_threshold)
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
        RoC_05q) 
    
    names(results_full_fin) <-  c("Working_Unit","Age", "ROC", "ROC_up", "ROC_dw")
    
    end_time <-  Sys.time()
    time_duration <-  end_time - start_time
    cat(paste(
      "R-RATEPOL finished", end_time, "taking", time_duration, units(time_duration)),
      fill = TRUE)
    
    return(results_full_fin)
    
  }
# end of code
