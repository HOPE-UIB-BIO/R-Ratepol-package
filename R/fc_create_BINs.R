fc_create_BINs <-function(data_source_bin, shift_value, Number_of_shifts)
{
  BIN <- shift_value*Number_of_shifts

  BIN.last <- ceiling(max(data_source_bin$Age$newage))
  BIN.breaks <- seq(from=0, to=BIN.last, by=BIN)

  BIN.breaks.temp <- BIN.breaks
  BIN.breaks.fin <- vector(mode = "numeric")


  for (j in 1:Number_of_shifts)
  {
    vector.w <- BIN.breaks.temp+shift_value*(j-1)
    BIN.breaks.fin <-c(BIN.breaks.fin,vector.w)
  }

  DF.sample.names <- data.frame(NAME = BIN.breaks.fin,
                                SHIFT = sort(rep(c(1:Number_of_shifts),length(BIN.breaks))))

  return(DF.sample.names)

}

