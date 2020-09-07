fc_smooth_pollen_data <- function(data.source,
                            smooth_method="age.w",
                            smooth_N_points = 5,
                            smooth_N_max = 9,
                            smooth_age_range = 500,
                            round_results = T,
                            Debug = F)
{
  # imput variables:
  # data.source - data prepared by the function of fn_extract_data
  # smooth_method = type of smoothing applied for the each of the pollen type
  #     "none"    = None: Pollen data is not smoothed
  #     "m.avg"   = Moving average: Each focus value is calculated as average over N(smooth_N_points) number of levels
  #                   (preferably ½ N before and ½ after focus level, those values are adjusted in the beginning
  #                   and end of the core. N must be odd number). Note that each calculation is done from scratch
  #                   and results are saved separately in order to avoid cumulative rounding errors.
  #     "grim"    = 	Grimm’s smoothing: Similar to moving average but N is not fixed. For each level, N is selected
  #                   as odd number between N_a (smooth_N_points) and N_b (smooth_N_max), while the maintaining the maximum age
  #                   difference from the selected levels as smooth_age_range.
  #     "age.w"   = 	Age weighted average:  Similar to moving average but average is weighted by the age difference
  #                   from the focus level and multiplied by 1/smooth_age_range. To avoid up-weighting levels,
  #                   if smooth_age_range/AGEDIFF exceeds 1, it is saved as 1. This means that levels closer than
  #                   smooth_age_range to the target age are given full weighting, but those farther away are downweighed
  #                   by an amount increasing with age difference.
  #     "shep"    =   Shepard's 5-term filter: Smoothing over 5 points following equation:
  #                   V_NEW=(17*V + 12*(V_((+1) )  + V_((-1) ) )-3*(V_((+2) )  + V_((-2) )))/35 ,
  #                   where V is focal level value. All values that result smaller than zero are saved as zero
  # smooth_N_points = Number of points for (need to be an odd number). Used for moving average, Grimm and Age-Weighted
  # smooth_N_max = maximal number of samples to look in Grimm smoothing
  # smooth_age_range = maximal age range for both Grimm and Age-weight smoothing

  # ----------------------------------------------
  #                     SETUP
  # ----------------------------------------------



  if(class(data.source) != "RRatepolList")
    stop("Data is not in RRatepolList format")

  # split data into 2 datasets
  p.counts <-  as.data.frame(data.source$Pollen)
  age <- as.data.frame(data.source$Age)
  focus.par <- matrix(data=NA,nrow=nrow(age),ncol=2)



  # ----------------------------------------------
  #               NONE SMOOTHING
  # ----------------------------------------------

  if(smooth_method=="none")
  {
    return(list(Pollen=p.counts, Age=age, Age.un=data.source$Age.un, Dim.val= data.source$Dim.val))
  }

  # ----------------------------------------------
  #                 SMOOTHING
  # ----------------------------------------------

  # check if smooth_N_points is and odd number
  if(smooth_N_points%%2 ==0)
    stop("smooth_N_points has to be an odd number")

  # check if smooth_N_max is an odd numbers
  if(smooth_N_max%%2 ==0)
    stop("smooth_N_max has to be an odd number")

  # Check if miminal number of points in not gigger than maximum
  if(smooth_N_points>smooth_N_max)
    stop("smooth_N_max has to be biger than smooth_N_points")


  if (Debug==T & smooth_method == "none")
    {cat("data will not be smoothed",fill=T)}
  if(Debug==T & smooth_method == "m.avg")
    {cat(paste("data will be smoothed by moving average over",smooth_N_points,"points"),fill=T)}
  if (Debug==T & smooth_method == "grim")
    {print(cat("data will be smoothed by Grimm method with min samples",smooth_N_points,
              "max samples",smooth_N_max,"and max age range of",smooth_age_range),fill=T)}
  if(Debug==T & smooth_method == "age.w")
  {cat(paste("data will be smoothed by age-weighed average over",smooth_N_points,"points"),fill=T)}
  if(Debug==T & smooth_method == "shep"){cat(paste("data will be smoothed by Shepard's 5-term filter"),fill=T)}

  # crete support fucntion for GRIMM smoothing
  search.parameter <- function(A, B, smooth_age_range)
    {
      # test if this increase does not invalidate rules.
      # 1) seach parameter cannot go outside of the sample size (up or down)
      # 2) seach parameter cannot be biger than selected maximum sample sizes
      # 3) the age diference between samples selected by the seach paramated cannot be higher than
      #  defined max age range
      # if all of those ARE TRUE then increase the real search parameter


      for (k in 1:(smooth_N_max-smooth_N_points))
      {
        # create new search parameter that is lower by 1
        A.test <- A-1
        if( A.test > 0 &  B-A.test < smooth_N_max) # i+N.active.test < nrow(p.counts) &
        { if (abs(age$newage[A.test]-age$newage[B])<smooth_age_range)
        {A <- A.test}
        }

        # create new search parameter that higher by 1
        B.test <- B+1
        if( B.test < nrow(p.counts) &  B-A.test < smooth_N_max)
        { if (abs(age$newage[A]-age$newage[B.test])<smooth_age_range)
        {B <- B.test}
        }
      }
      return(c(A,B))
    }

 # ----------------------------------------------
 #                   CALCULATION
 # ----------------------------------------------

  for(j in 1:ncol(p.counts)) # for every species
  {
    col.work <- .subset2(p.counts,j) # select the species
    col.res <- rep(0,length(col.work)) # create empty vector of same lengt for values to be saved

    for(i in 1:nrow(p.counts)) # for each sample
    {

      # ----------------------------------------------
      #           MOVING AVERAGE SMOOTHING
      # ----------------------------------------------
      if(smooth_method=="m.avg")
      {
        if( i < round(0.5*(smooth_N_points))+1 ) {  # Samples near beginning (moving window truncated)
          focus.par[i,] = c(1, ( i + round(0.5*(smooth_N_points)) ))
        } else {
          if( i > nrow(age)-round(0.5*(smooth_N_points)) ) { # Samples near end
            focus.par[i,] = c( (i - round(0.5*(smooth_N_points))), nrow(age) )
          } else {
            focus.par[i,] = c( (i - round(0.5*(smooth_N_points))), (i+round(0.5*(smooth_N_points))) )
          }
        }
        col.res[i]<- mean (col.work[focus.par[i,1]:focus.par[i,2]])
      }


      # ----------------------------------------------
      #               GRIMm SMOOTHING
      # ----------------------------------------------
      if(smooth_method == "grim")
      {
        if( i < round(0.5*(smooth_N_max))+1 ) {  # Samples near beginning (moving window truncated)
          focus.par[i,1] <- 1
          focus.par[i,2] <- ( i + round(0.5*(smooth_N_points)) )
          focus.par[i,] <- search.parameter(focus.par[i,1] ,focus.par[i,2], smooth_age_range  )
        } else {
          if( i > nrow(age)-round(0.5*(smooth_N_points)) ) { # Samples near end
            focus.par[i,1] <- (i - round(0.5*(smooth_N_points)))
            focus.par[i,2] <- nrow(age)
            focus.par[i,] <- search.parameter(focus.par[i,1] ,focus.par[i,2], smooth_age_range  )

          } else {
            focus.par[i,1] <- (i - round(0.5*(smooth_N_points)))
            focus.par[i,2] <- (i + round(0.5*(smooth_N_points)))
            focus.par[i,] <- search.parameter(focus.par[i,1] ,focus.par[i,2], smooth_age_range  )
          }
        }
        col.res[i]<- mean(col.work[focus.par[i,1]:focus.par[i,2]])
      }

      # ----------------------------------------------
      #           AGE-WEIGHTED SMOOTHING
      # ----------------------------------------------
      if(smooth_method=="age.w")
      {

        if( i < round(0.5*(smooth_N_points))+1 ) {  # Samples near beginning (moving window truncated)
          focus.par[i,] = c(1, ( i + round(0.5*(smooth_N_points)) ))
        } else {
          if( i > nrow(age)-round(0.5*(smooth_N_points)) ) { # Samples near end
            focus.par[i,] = c( (i - round(0.5*(smooth_N_points))), nrow(age) )
          } else {
            focus.par[i,] = c( (i - round(0.5*(smooth_N_points))), (i+round(0.5*(smooth_N_points))) )
          }
        }

      # create small df with values around observed sample (in range of offset)
      df.work <-  data.frame(values= col.work[focus.par[i,1]:focus.par[i,2]],
                             age = age$newage[focus.par[i,1]:focus.par[i,2]],
                             Weight=1)

      # Weith of points is calculated as smooth_age_range / distance bewtween oldest and youngest points.
      # If cannot be smaller than 1. Values very far away from the point
      F.age.dist <- abs(df.work$age-age$newage[i])
      const <-  smooth_age_range/F.age.dist
      const[const>1] <- 1
      df.work$Weight <- const

      col.res[i]<-weighted.mean(df.work$values,df.work$Weight)

      }

      # ----------------------------------------------
      #             Shepard's 5-term filter
      # ----------------------------------------------
      if(smooth_method=="shep")
      {
        if(i < round(0.5*(smooth_N_points))+1) {
          col.res[i] <- col.work[i]
        } else {
          if (i > nrow(age)-round(0.5*(smooth_N_points)) ) {
            col.res[i] <- col.work[i]
          } else {
            w.value  <- (17*.subset(col.work,i) + 12*(.subset(col.work,i+1)+.subset(col.work,i-1)) - 3*(.subset(col.work,i+2)+.subset(col.work,i-2))) / 35
            if(w.value<0){w.value<-0}
            col.res[i] <- w.value
          }
        }

      }

    }
    p.counts[,j]<-col.res
  }

  if (round_results == T){
    p.counts <- round(p.counts)
  }

  final_list <- list(Pollen=p.counts, Age=age, Age.un=data.source$Age.un, Dim.val= data.source$Dim.val)
  class(final_list) = "RRatepolList"

  return(final_list)
}
