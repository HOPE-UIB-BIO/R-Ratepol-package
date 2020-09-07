fc_subset_samples <- function(data_subset, BINS)
{
  # create empty template with length of number of BINS for both age and pollen

  data.res.age <- data.frame(matrix(ncol = ncol(data_subset$Age), nrow = nrow(BINS)))
  names(data.res.age) <- names(data_subset$Age)

  data.res.pollen <- data.frame(matrix(ncol=ncol(data_subset$Pollen), nrow = nrow(BINS)))
  names(data.res.pollen) <- names(data_subset$Pollen)

  row.names(data.res.age) <- BINS$NAME
  row.names(data.res.pollen) <- BINS$NAME

  # calculate bin size
  BIN.size <- BINS$NAME[2]-BINS$NAME[1]

  for(i in 1:nrow(BINS)) # for each bin
  {
    selected.BIN <- BINS$NAME[i] #select the bin

    # subset age data so it selected all samples which has higher values than the BIN itself but
    # still small then selected bin + calculated BIN size
    subset.w <- data_subset$Age[data_subset$Age$newage < BINS$NAME[i]+BIN.size &
                                  data_subset$Age$newage > BINS$NAME[i],]

    if (nrow(subset.w)>0) # If selected subset has at least one sample
    {

      subset.w$diff <- abs(subset.w$newage-selected.BIN)
      suppressWarnings(data.res.age[i,] <- subset.w[subset.w$diff==min(subset.w$diff),c(1:4)] )

      data.res.pollen[i,]<-  data_subset$Pollen[row.names(data_subset$Pollen) %in% data.res.age$sample.id[i],]
    }
  }

  list.res <- list(Pollen = data.res.pollen, Age=data.res.age,Dim.val=data_subset$Dim.val )
  list.res <- fc_check_data(list.res, proportion = F, Debug = F)
  class(list.res) <- "RRatepolList"

  return(list.res)
}
