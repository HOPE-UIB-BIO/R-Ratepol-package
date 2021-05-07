fc_smooth_community_data <- function(data_source,
                                     smooth_method = "none",
                                     smooth_N_points = 5,
                                     smooth_N_max = 9,
                                     smooth_age_range = 500,
                                     round_results = TRUE,
                                     Debug = FALSE){
  
  # ----------------------------------------------
  # SETUP -----
  # ----------------------------------------------
  
  if(class(data_source) != "RRatepolList")
    stop("Data is not in RRatepolList format")
  
  # split data into 2 datasets
  dat_community <-  as.data.frame(data_source@Community)
  age <- as.data.frame(data_source@Age)
  focus_par <- matrix(data = NA,nrow = nrow(age), ncol = 2)
  
  
  assertthat::assert_that(
    any(smooth_method == c("none", "m.avg", "grim", "age.w", "shep")),
    msg = "`smooth_method` must be one of the following:
    `none`, `m.avg`, `grim`, `age.w`, `shep`")
  
  
  # ----------------------------------------------
  # NONE SMOOTHING -----
  # ----------------------------------------------
  
  if (Debug == TRUE & smooth_method == "none"){
    cat("data will not be smoothed", fill = TRUE)}
  
  if(smooth_method == "none"){
    return(
      RRatepolList(
        Community = dat_community,
        Age = age, 
        Age.un = data_source@Age.un,
        Dim.val = data_source@Dim.val))
  }
  
  
  # ----------------------------------------------
  # Additional information -----  
  # ----------------------------------------------
  
  if(Debug == TRUE){
    
    if(smooth_method == "m.avg"){
      cat(paste(
        "Data will be smoothed by moving average over", smooth_N_points,
        "points"), fill = TRUE)  
    }
    
    if(smooth_method == "grim"){
      cat(paste(
        "Data will be smoothed by Grimm method with min samples", smooth_N_points,
        "max samples", smooth_N_max, "and max age range of", smooth_age_range),
        fill = TRUE)  
    }
    
    if(smooth_method == "age.w"){
      cat(paste(
        "Data will be smoothed by age-weighed average over",smooth_N_points,"points",
        "with a threshod of",smooth_age_range),
        fill = TRUE)  
    }
    
    if(smooth_method == "shep"){
      cat(paste(
        "Data will be smoothed by Shepard's 5-term filter"),
        fill = TRUE)  
    }
  }
  
  
  # crete support fucntion for GRIMM smoothing
  search.parameter <- function(A, B, smooth_age_range){
    
    # test if this increase does not invalidate rules.
    # 1) seach parameter cannot go outside of the sample size (up or down)
    # 2) seach parameter cannot be biger than selected maximum sample sizes
    # 3) the age diference between samples selected by the seach paramated cannot be higher than
    #  defined max age range
    # if all of those ARE TRUE then increase the real search parameter
    
    for (k in 1:(smooth_N_max-smooth_N_points)){
      # create new search parameter that is lower by 1
      A_test <- A - 1
      if( A_test > 0 &  B - A_test < smooth_N_max){ # i+N.active.test < nrow(dat_community) &
        if (abs(age$newage[A_test] - age$newage[B]) < smooth_age_range){
          A <- A_test}
      }
      
      # create new search parameter that higher by 1
      B_test <- B + 1
      if( B_test < nrow(dat_community) &  B - A_test < smooth_N_max){
        if (abs(age$newage[A] - age$newage[B_test]) < smooth_age_range){
          B <- B_test}
      }
    }
    return(c(A, B))
  }
  
  # ----------------------------------------------
  # CALCULATION -----
  # ----------------------------------------------
  
  
  for(j in 1:ncol(dat_community)){ # for every species
    
    col_work <- .subset2(dat_community, j) # select the species
    col_res <- rep(0, length(col_work)) # create empty vector of same lengt for values to be saved
    
    for(i in 1:nrow(dat_community)){ # for each sample
      
      # ----------------------------------------------
      # MOVING AVERAGE SMOOTHING -----
      # ----------------------------------------------
      if(smooth_method == "m.avg"){
        
        # check if smooth_N_points is and odd number
        assertthat::assert_that(
          smooth_N_points%%2 != 0,
          msg = "`smooth_N_points` must be an odd number")
        
        if( i < round(0.5 * (smooth_N_points)) + 1 ){  # Samples near beginning (moving window truncated)
          focus_par[i, ] = c(1, ( i + round(0.5 * (smooth_N_points)) ))
        } else {
          if( i > nrow(age)-round(0.5 * (smooth_N_points)) ){ # Samples near end
            focus_par[i, ] = c( (i - round(0.5 * (smooth_N_points))), nrow(age) )
          } else {
            focus_par[i, ] = c( (i - round(0.5 * (smooth_N_points))), (i + round(0.5 * (smooth_N_points))) )
          }
        }
        col_res[i] <-  mean(col_work[focus_par[i, 1]:focus_par[i, 2]])
      }
      
      
      # ----------------------------------------------
      # GRIMMM SMOOTHING -----
      # ----------------------------------------------
      if(smooth_method == "grim"){
        
        # check if smooth_N_points is and odd number
        assertthat::assert_that(
          smooth_N_points%%2 != 0,
          msg = "`smooth_N_points` must be an odd number")
        
        # check if smooth_N_max is an odd numbers
        assertthat::assert_that(
          smooth_N_max%%2 != 0,
          msg = "`smooth_N_max` must be an odd number")
        
        # Check if miminal number of points in not bigger than maximum
        assertthat::assert_that(
          smooth_N_points < smooth_N_max,
          msg = "`smooth_N_max` must be bigger than `smooth_N_points")
        
        assertthat::assert_that(
          is.numeric(smooth_age_range),
          msg = "`smooth_age_range` must be `numeric")
        
        if( i < round(0.5 * (smooth_N_max)) + 1 ) {  # Samples near beginning (moving window truncated)
          focus_par[i, 1] <-  1
          focus_par[i, 2] <-  ( i + round(0.5 * (smooth_N_points)) )
          focus_par[i, ] <-  search.parameter(focus_par[i, 1], focus_par[i, 2], smooth_age_range)
        } else {
          if( i > nrow(age) - round(0.5 *(smooth_N_points)) ) { # Samples near end
            focus_par[i, 1] <-  (i - round(0.5 * (smooth_N_points)))
            focus_par[i, 2] <-  nrow(age)
            focus_par[i, ] <-  search.parameter(focus_par[i, 1], focus_par[i, 2], smooth_age_range)
            
          } else {
            focus_par[i, 1] <-  (i - round(0.5 * (smooth_N_points)))
            focus_par[i, 2] <-  (i + round(0.5 * (smooth_N_points)))
            focus_par[i, ] <-  search.parameter(focus_par[i, 1], focus_par[i, 2], smooth_age_range)
          }
        }
        col_res[i] <-  mean(col_work[focus_par[i,1]:focus_par[i,2]])
      }
      
      # ----------------------------------------------
      # AGE-WEIGHTED SMOOTHING -----
      # ----------------------------------------------
      if(smooth_method == "age.w"){
        
        # check if smooth_N_points is and odd number
        assertthat::assert_that(
          smooth_N_points%%2 != 0,
          msg = "`smooth_N_points` must be an odd number")
        
        assertthat::assert_that(
          is.numeric(smooth_age_range),
          msg = "`smooth_age_range` must be `numeric")
        
        if( i < round(0.5 * (smooth_N_points)) + 1 ) {  # Samples near beginning (moving window truncated)
          focus_par[i, ] <-  c(1, ( i + round(0.5 * (smooth_N_points)) ))
        } else {
          if( i > nrow(age) - round(0.5 * (smooth_N_points)) ) { # Samples near end
            focus_par[i, ] <-  c( (i - round(0.5 * (smooth_N_points))), nrow(age) )
          } else {
            focus_par[i, ] <-  c( (i - round(0.5 * (smooth_N_points))), (i + round(0.5 * (smooth_N_points))) )
          }
        }
        
        # create small df with values around observed sample (in range of offset)
        df_work <-  
          data.frame(values = col_work[focus_par[i, 1]:focus_par[i, 2]],
                     age = age$newage[focus_par[i, 1]:focus_par[i, 2]],
                     Weight = 1)
        
        # Weith of points is calculated as smooth_age_range / distance bewtween oldest and youngest points.
        # If cannot be smaller than 1. Values very far away from the point
        F_age_dist <-  abs(df_work$age - age$newage[i])
        const <-  smooth_age_range/F_age_dist
        const[const > 1] <-  1
        df_work$Weight <-  const
        
        col_res[i] <-  stats::weighted.mean(df_work$values,df_work$Weight)
        
      }
      
      # ----------------------------------------------
      # Shepard's 5-term filter -----
      # ----------------------------------------------
      if(smooth_method == "shep"){
        
        if(i < round(0.5 * (smooth_N_points)) + 1) {
          col_res[i] <-  col_work[i]
        } else {
          if (i > nrow(age) - round(0.5 * (smooth_N_points))) {
            col_res[i] <-  col_work[i]
          } else {
            w.value  <- (17 * .subset(col_work, i) + 12 * (.subset(col_work, i + 1)+.subset(col_work, i - 1)) - 3 * (.subset(col_work,i + 2) + .subset(col_work, i - 2))) / 35
            if(w.value < 0){w.value <-  0}
            col_res[i] <-  w.value
          }
        }
        
      }
      
    }
    dat_community[ ,j] <-  col_res
  }
  
  if (round_results == TRUE){
    dat_community <-  round(dat_community)
  }
  
  final_list <-  
    RRatepol:::RRatepolList(
      Community = dat_community,
      Age = age,
      Age.un = data_source@Age.un,
      Dim.val = data_source@Dim.val)
  
  return(final_list)
}
