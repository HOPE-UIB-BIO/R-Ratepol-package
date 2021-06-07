fc_create_bins <-function(data_source_bin, shift_value, Number_of_shifts){
  
  bin <-  shift_value * Number_of_shifts

  bin_last <-  ceiling(max(data_source_bin@Age$newage))
  bin_first <- min(data_source_bin@Age$newage)
  bin_breaks <-  seq(
    from = bin_first,
    to = bin_last,
    by = bin)

  bin_breaks_temp <-  bin_breaks
  bin_breaks_fin <-  vector(mode = "numeric")

  for (j in 1:Number_of_shifts){
    
    vector_w <-  bin_breaks_temp + shift_value*(j - 1)
    bin_breaks_fin <-  c(bin_breaks_fin, vector_w)
  }

  df_sample_names <- 
    data.frame(
      name = bin_breaks_fin,
      shift = sort(rep(c(1:Number_of_shifts), length(bin_breaks))))

  return(df_sample_names)

}

