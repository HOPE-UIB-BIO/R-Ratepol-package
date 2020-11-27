#function for extracting results
fc_extract_result <- function (x, sel_var, rand){

  # cretae pivot table randomisation ID by BIN code and order it by BIN codes
  r_m <-
    x %>%
    dplyr::mutate(bin_shift = paste0(shift, "_", bin)) %>%
    dplyr::select(c(ID, bin, shift, sel_var)) %>%
    tidyr::pivot_wider(names_from = "ID", values_from = sel_var) %>%
    dplyr::mutate(bin_order = sub(" -.*","", bin) %>%
                    as.numeric()) %>%
    dplyr::arrange(.,bin_order) %>%
    dplyr::select(-c(bin_order))

  # calculate number of NOT NA values for each BIN code
  N_not_NA <-
    r_m  %>%
    dplyr::select(-c(bin)) %>%
    apply(., 1, FUN = function(x){
      y <-  is.na(x) %>%
        table ()
      return(y[1])
    })

  # include only bin codes if they have ast least 10% of number of randomisation
  r_m_sel <- 
    r_m %>%
    dplyr::filter(N_not_NA > 0.1 * rand)

  # calculate median and 95% quantile for selected variable for each sample (BIN)
  r_m_sum <- 
    r_m_sel %>%
    dplyr::mutate(
      sample_id = bin,
      shift = shift,
      !!sel_var := dplyr::select(., -c(bin,shift)) %>%
        apply(., 1, FUN = function(x) median(x, na.rm = T)),
      !!paste0(sel_var,"_05q") := dplyr::select(., -c(bin,shift))
      %>% apply(., 1, FUN = function(x) quantile(x, 0.025, na.rm = T)),
      !!paste0(sel_var,"_95q") := dplyr::select(., -c(bin,shift))
      %>% apply(., 1, FUN = function(x) quantile(x, 0.975, na.rm = T))
    ) %>%
    dplyr::select(
      sample_id,
      shift, 
      sel_var, 
      paste0(sel_var,"_05q"), 
      paste0(sel_var,"_95q"))

  return(r_m_sum)
}
