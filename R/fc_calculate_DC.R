fc_calculate_DC <- function (data_source_DC, DC = "chord", Debug = F)
{
  # DC = disimilarity coeficient. Type of calculation of differences between samples/BINs
  #   "euc"     = 	Euclidean distance: √(∑_(i=1)^n〖(A_i-B_i)〗^2 ), where A_i and B_i are values for Working Unit A and B given for species i.
  #   "euc.sd"  = 	Standardised Euclidean distance: √(∑_(i=1)^n〖((A_i-B_i)/〖SD〗_i )〗^2 ) ,
  #                 where A_i and B_i are values for Working Unit A and B given for species i,and 〖SD〗_i is a standard deviation for species i,
  #                 calculated from whole sequence.
  #   "chord"   = 	Chord distance:√(∑_(i=1)^n〖(√(A_i )-√(B_i ))〗^2 ), where A_i and B_i are values for Working Unit A and B given for species i.
  #   "chisq    = 	Chi-squared coefficient:√(∑_(i=1)^n〖(A_i-B_i)〗^2/((A_i+B_i))), where A_i and B_i are values for Working Unit A and B
  #                 given for species i
  #
  # Comment 1:
  # Note that DC is calculated between each subsequent Working Units and then standardise by the time difference between them. Age of each RoC values
  #   is determined as mean between values of Wokring Units.

  dat.res <- vector(mode="numeric",length = data_source_DC$Dim.val[2]-1)

  # ----------------------------------------------
  #               EUCLIDAN DISTANCE
  # ----------------------------------------------

  if (DC == "euc")
  {
    if (Debug==T){cat("Euclidan distance will be used as DC",fill=T)}
    for (i in 1:(data_source_DC$Dim.val[2]-1)) # for each sample (except the last)
    {
      df.work<- data_source_DC$Pollen[c(i,i+1),] # select only 2 samples (observed + 1 after)

      df.work<-as.data.frame(df.work[,colSums(df.work)>0]) # get rid of "empty species"

      vector.work <- vector(mode="numeric", length = ncol(df.work)) # vector for result for each species

      for( j in 1:ncol(df.work)) # for each species
      {
        a<- .subset2(df.work,j)[1]
        b<- .subset2(df.work,j)[2]
        vector.work[j] <- (a-b)**2 # calculate the diference
      }

      dat.res[i]<- sqrt(sum(vector.work)) # save the square root of sum of all dufereces

    }

  }

  # ----------------------------------------------
  #       STANDARDISED EUCLIDAN DISTANCE
  # ----------------------------------------------

  if (DC=="euc.sd")
  {
    if(Debug==T){cat("Standardised Euclidan distance will be used as DC",fill=T)}

    # calculation of standard deviation for each species

    # vector for standar deviation for each species
    df.sp.supp <- vector(mode="numeric",length = data_source_DC$Dim.val[1])
    # calculate the SD for each species
    df.sp.supp <- apply(data_source_DC$Pollen,2,sd)

    # calculation of the DC
    for (i in 1:(data_source_DC$Dim.val[2]-1)) # for each sample (except the last)
    {
      #print(paste("i",i))
      df.work<- data_source_DC$Pollen[c(i,i+1),] # select only 2 samples (observed + 1 after)

      # get rid of "empty species" in data & in sp.std
      df.sp.supp.work<- df.sp.supp[colSums(df.work)>0]
      df.work<-as.data.frame(df.work[,colSums(df.work)>0])

      vector.work <- vector(mode="numeric", length = ncol(df.work)) # vector for result for each species

      for( j in 1:ncol(df.work)) # for each species
      {
        #print(paste("j",j))
        if (df.sp.supp.work[j]!=0) # check if the standard deviation is not equal zero
        {
          a<- .subset2(df.work,j)[1]
          b<- .subset2(df.work,j)[2]
          vector.work[j] <- ((a-b)/df.sp.supp.work[j])**2 # calculate the diference
        }

      }

      dat.res[i]<- sqrt(sum(vector.work)) # save the square root of sum of all dufereces

    }

  }


  # ----------------------------------------------
  #               CHORD DISTANCE
  # ----------------------------------------------

  if (DC == "chord")
  {
    if(Debug==T){cat("Chord distance will be used as DC",fill=T)}

    for (i in 1:(data_source_DC$Dim.val[2]-1)) # for each sample (except the last)
    {
      df.work<- data_source_DC$Pollen[c(i,i+1),] # select only 2 samples (observed + 1 after)

      df.work<-as.data.frame(df.work[,colSums(df.work)>0]) # get rid of "empty species"

      vector.work <- vector(mode="numeric", length = ncol(df.work)) # vector for result for each species

      for( j in 1:ncol(df.work)) # for each species
      {
        a<- .subset2(df.work,j)[1]
        b<- .subset2(df.work,j)[2]
        vector.work[j] <- (sqrt(a)-sqrt(b))**2 # calculate the diference
      }

      dat.res[i]<- sqrt(sum(vector.work)) # save the square root of sum of all dufereces

    }

  }


  # ----------------------------------------------
  #           CHI-SQUARED COEFICIENT
  # ----------------------------------------------

  if (DC == "chisq")
  {
    if(Debug==T){cat("Chi-squared coeficient will be used as DC",fill=T)}

    for (i in 1:(data_source_DC$Dim.val[2]-1)) # for each sample (except the last)
    {
      df.work<- data_source_DC$Pollen[c(i,i+1),] # select only 2 samples (observed + 1 after)

      df.work<-as.data.frame(df.work[,colSums(df.work)>0]) # get rid of "empty species"

      vector.work <- vector(mode="numeric", length = ncol(df.work)) # vector for result for each species

      for( j in 1:ncol(df.work)) # for each species
      {
        a<- .subset2(df.work,j)[1]
        b<- .subset2(df.work,j)[2]
        vector.work[j] <- ((a-b)**2) / (a+b) # calculate the diference
      }

      dat.res[i]<- sqrt(sum(vector.work)) # save the square root of sum of all dufereces

    }

  }

  # ----------------------------------------------
  #                   RESULT
  # ----------------------------------------------

  return(dat.res)


}
