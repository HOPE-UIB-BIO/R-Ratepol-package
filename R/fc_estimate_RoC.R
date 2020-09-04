fc_estimate_RoC <- function (data_source_pollen,
                             data_source_age,
                             smooth_method = "shep",
                             smooth_N_points = 5,
                             smooth_age_range = 500,
                             smooth_N_max = 9,
                             Working_Units = "MW",
                             bin_size = 500,
                             Number_of_shifts = 5,
                             rand = 10000,
                             treads = T,
                             age_uncertainty = T,
                             standardise = T,
                             N_pollen_grains = 150,
                             DC = "chord",
                             interest_threshold = F,
                             Peak_detection = "GAM",
                             Debug = F)
{
  # CODE
  start.time <- Sys.time()
  if (Debug == T){
    cat(paste("RATEPOL started", start.time), fill=T)
  }

  # ----------------------------------------------
  #               DATA EXTRACTION
  # ----------------------------------------------

  # extract data into working format
  # already include data check
  data.extract <- fc_extract_data(data_source_pollen,data_source_age, Debug=Debug, age_uncertainty = age_uncertainty)

  # ----------------------------------------------
  #               POLLEN SMOOTHING
  # ----------------------------------------------
  data.smooth <- data.extract

  # smooth pollen data by selected smoothing type
  data.smooth <- fc_smooth_pollen_data(data.smooth,
                                       smooth_method = smooth_method,
                                       smooth_N_points = smooth_N_points,
                                       smooth_N_max = smooth_N_max,
                                       smooth_age_range = smooth_age_range,
                                       round_results = standardise,
                                       Debug=Debug)

  #check data and reduce data dimentions
  data.work <- fc_check_data(data.smooth, proportion = F, Debug=Debug)

  # ----------------------------------------------
  #               Working Unit selection
  # ----------------------------------------------

  if(Working_Units == "levels"){
    Number_of_shifts <- 1
  }

  if(Working_Units == "BINs"){
    Number_of_shifts <- 1
    BIN_SIZES <- fc_create_BINs(data.work,shift_value = bin_size, Number_of_shifts = 1)
  }

  if(Working_Units == "MW"){
    shift_value <- bin_size/Number_of_shifts
    BIN_SIZES <- fc_create_BINs(data.work,shift_value = shift_value, Number_of_shifts = Number_of_shifts)
  }


  # ----------------------------------------------
  #             RANDOMOMIZATION
  # ----------------------------------------------

  # select the prefetred number of cores for of cores for parallel computation
  if(class(treads) == "numeric"){
    Ncores <- treads # set value
  } else {
    if (treads == T) {
      Ncores <- parallel::detectCores() # detect number
    } else {
      Ncores <- 1
    }
  }


  # create cluster
  cl <- parallel::makeCluster(Ncores)
  doSNOW::registerDoSNOW(cl)

  # add all functions to the cluster
  envir <- environment(fc_estimate_RoC)
  parallel::clusterExport(cl, varlist = c(ls(envir)))

  # create progress bar based os the number of replication
  pb <- utils::txtProgressBar(max = rand, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  # create template for result tibble
  SHIFT.tibble_template <- tibble::tibble()

  result.tibble <-   foreach::`%dopar%`(foreach::foreach(l=1:rand,.combine = rbind, .options.snow=opts), {

    # TIME SAMPLING
    # sample random time sequence from time uncern.
    data.work$Age$newage <- as.numeric(data.work$Age.un[sample(c(1:max(1,nrow(data.work$Age.un))),1),])

    # create result tible
    SHIFT.tibble <- SHIFT.tibble_template

    # repeat for number of shifts
    for(k in 1: Number_of_shifts)
    {
      # ----------------------------------------------
      #             DATA SUBSETTING
      # ----------------------------------------------
      data.subset <- data.work

      # select one sample for each bin based on the age of the samples. Sample is chones if it is the closes one to the upper end of the BIN
      if (Working_Units != "levels")
      {
        # select BIN for this shift
        SELECTED.BINS <- BIN_SIZES[BIN_SIZES$SHIFT==k,]

        #subset data
        data.subset <- fc_subset_samples(data.subset,SELECTED.BINS)
      }

      # ----------------------------------------------
      #         DATA STANDARDISATION
      # ----------------------------------------------
      data.sd <-data.subset

      # standardisation of pollen data to X(N_pollen_grains) number of pollen grains
      if(standardise==T) #
      {

        # adjust the Svalue by the minimal Pollen or to a minimal of presected values
        N_pollen_grains <-  min(c( rowSums(data.subset$Pollen),N_pollen_grains))

        # check if all samples has N_pollen_grains of pollen grains
        data.sd$Age <- data.sd$Age[rowSums(data.sd$Pollen, na.rm = T)>=N_pollen_grains,]
        data.sd$Age.un <- data.sd$Age.un[,rowSums(data.sd$Pollen, na.rm = T)>=N_pollen_grains]
        data.sd$Pollen <- data.sd$Pollen[rowSums(data.sd$Pollen, na.rm = T)>=N_pollen_grains,]
        data.sd<- fc_check_data(data.sd, proportion = F, Samples = T, Debug=Debug)

        # standardisation
        data.sd <- fc_standardise_pollen_data(data.sd, N_pollen_grains, Debug=Debug)

        if(any(rowSums(data.sd$Pollen, na.rm = T)!=N_pollen_grains))
          stop("standardisation was unsuccesfull")
      }

      # data check with proportioning
      data.sd.check <- fc_check_data(data.sd, proportion = T, Samples = F, Debug=Debug)

      # ----------------------------------------------
      #               DC CALCULATION
      # ----------------------------------------------

      # calculate DC between each subsequent samples/BINs
      DC.res <- fc_calculate_DC(data.sd.check,DC=DC, Debug=Debug)

      # ----------------------------------------------
      #             AGE STANDARDISATION
      # ----------------------------------------------

      # create empty vector with size = numeber of samples-1
      sample.size.work <- data.sd.check$Dim.val[2]-1

      # create empty vectors for age difference calcualtion
      age.diff <- vector(mode = "numeric", length = sample.size.work )
      age.diff.names <- vector(mode = "character", length = sample.size.work )
      age.mean <- age.diff

      for (i in 1:sample.size.work) # for each RoC
      {
        # calcualte the age difference between subsequesnt samples
        age.diff[i] <- data.sd.check$Age$newage[i+1]-data.sd.check$Age$newage[i]

        # Set age difference as 1, if age difference between samples is smaller than 1
        if(age.diff[i]<1)
        {age.diff[i]<-1}

        #calculate the average position of RoC
        age.mean[i] <- mean(c(data.sd.check$Age$age[i+1],data.sd.check$Age$age[i]))

        # create vector with bin names
        age.diff.names[i] <- paste(row.names(data.sd.check$Age)[i],"-",row.names(data.sd.check$Age)[i+1])
      }

      if (Debug ==T)
      {
        cat("", fill=T)
        cat(paste("The time standardisation unit (TSU) is",round(mean(age.diff),2)), fill=T)
      }

      #  calculate DC standardise by time
      DC.res.s <- vector(mode = "numeric", length = sample.size.work)
      DC.res.s <- (DC.res*mean(age.diff))/age.diff

      # ----------------------------------------------
      #             Result of single SHIFT
      # ----------------------------------------------

      # add the results from this shift into the result tibble
      SHIFT.tibble <- rbind(SHIFT.tibble,
                            data.frame(BIN= age.diff.names,
                                       DC = DC.res,
                                       Age.Pos = age.mean,
                                       Age.Diff = age.diff,
                                       RoC = DC.res.s,
                                       SHIFT = rep(k,sample.size.work))
      )
    }

    # ----------------------------------------------
    #         RESULT OF SINGLE RAND RUN
    # ----------------------------------------------

    # save result from single randomisation into data.frame with number of randomisation as ID.
    data.result.temp <- as.data.frame(list(ID=l,RUN=SHIFT.tibble))

    return(data.result.temp)
  })# end of the randomization

  # close progress bar and cluster
  close(pb)
  parallel::stopCluster(cl)


  # ----------------------------------------------
  #             RESULTs SUMMARY
  # ----------------------------------------------

  # create new dataframe with summary of randomisation results

  # extract results and match them by BIN
  r.m.full <- dplyr::right_join(fc_extract_result(result.tibble,"RUN.RoC", rand),
                                dplyr::select(fc_extract_result(result.tibble,"RUN.Age.Pos", rand),-shift),
                                by="sample.id")


  # reduce results by the focus age time
  if (interest_threshold!=F){
    r.m.full <- dplyr::filter(r.m.full,RUN.Age.Pos<=interest_threshold)
  }

  # sort samples by age and add smoothing to avoid "waiving"
  r.m.full <- r.m.full[order(r.m.full$RUN.Age.Pos),]
  r.m.full$RUN.RoC.sm = stats::lowess(r.m.full$RUN.Age.Pos,r.m.full$RUN.RoC,f=.1,iter=100)$y
  r.m.full$RUN.RoC.95q.sm = stats::lowess(r.m.full$RUN.Age.Pos,r.m.full$RUN.RoC.95q,f=.1,iter=100)$y
  r.m.full$RUN.RoC.05q.sm = stats::lowess(r.m.full$RUN.Age.Pos,r.m.full$RUN.RoC.05q,f=.1,iter=100)$y
  r.m.full$RUN.RoC.sm = ifelse(r.m.full$RUN.RoC.sm<=0,0.0001,r.m.full$RUN.RoC.sm)
  r.m.full$RUN.RoC.05q.sm = ifelse(r.m.full$RUN.RoC.05q.sm<0,0.0001,r.m.full$RUN.RoC.05q.sm)


  # ----------------------------------------------
  #             PEAK DETECTION
  # ----------------------------------------------

  # Median peak treshold
  if(Peak_detection == "Threshold"){
    # treshold for RoC peaks is set as median of all RoC in dataset
    r.treshold <- median(r.m.full$RUN.RoC.sm)
    # mark peaks which have 95% quantile above the treshold as Peak
    r.m.full$Peak = r.m.full$RUN.RoC.05q.sm > r.treshold
  }

  # GAM
  if(Peak_detection == "GAM"){
    # mark points that are abowe the GAM model (exactly 1.5 SD higher than GAM prediction)
    r.m.full$pred.gam = mgcv::predict.gam(mgcv::gam(RUN.RoC.sm~s(RUN.Age.Pos,k=3), data = r.m.full, family = "Gamma()",method = "REML"), type="response")
    r.m.full$pred.gam.diff = r.m.full$RUN.RoC.sm - r.m.full$pred.gam
    r.m.full$Peak = (r.m.full$pred.gam.diff) > 1.5*sd(r.m.full$pred.gam.diff)

  }

  # SNI
  if (Peak_detection == "SNI"){
    # set moving window of 5 times higher than average distance between samples
    mean.age.window <- 5 * mean( diff(r.m.full$RUN.Age.Pos) )
    # create GAM
    pred.gam <-  mgcv::predict.gam(mgcv::gam(RUN.RoC.sm~s(RUN.Age.Pos,k=3), data = r.m.full, family = "Gamma()",
                                             method = "REML"), type="response")
    # calculate SNI (singal to noise ratio)
    SNI.calc <- fc_CharSNI(data.frame(r.m.full$RUN.Age.Pos, r.m.full$RUN.RoC.sm, pred.gam),mean.age.window)
    # mark points with SNI higher than 3
    r.m.full$Peak = SNI.calc$SNI > 3 & r.m.full$RUN.RoC.sm > pred.gam
  }

  # outro
  r.m.full.fin <- dplyr::select(r.m.full, sample.id,RUN.RoC.sm,RUN.RoC.95q.sm,RUN.RoC.05q.sm,RUN.Age.Pos,Peak)
  names(r.m.full.fin) <- c("Working_Unit","ROC","ROC.up","ROC.dw","AGE","PEAK")

  end.time <- Sys.time()
  time.length <- end.time - start.time
  if(Debug == T){
    cat("", fill=T)
    cat(paste("RATEPOL finished", end.time,"taking",time.length, units(time.length)), fill=T)
  }


  return(r.m.full.fin)

}
#end of CODE
