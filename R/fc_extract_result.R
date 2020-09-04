#function for extracting results
fc_extract_result <- function (x,sel.var, rand)
{

  # cretae pivot table randomisation ID by BIN code and order it by BIN codes
  r.m<-  x %>%
    dplyr::mutate(BIN_shift = paste0(RUN.SHIFT,"_",RUN.BIN)) %>%
    dplyr::select(c(ID,RUN.BIN,RUN.SHIFT,sel.var)) %>%
    tidyr::pivot_wider(names_from = "ID", values_from = sel.var) %>%
    dplyr::mutate(BIN = sub(" -.*","",RUN.BIN) %>% as.numeric()) %>%
    dplyr::arrange(.,BIN) %>%
    dplyr::select(-c(BIN))

  # calculate number of NOT NA values for each BIN code
  N.notNA<-  r.m  %>%
    dplyr::select(-c(RUN.BIN)) %>%
    apply(.,1, FUN = function(x){
      y<- is.na(x) %>% table ()
      return(y[1])
    })

  # include only bin codes if they have ast least 10% of number of randomisation
  r.m.sel <- r.m %>%
    dplyr::filter(N.notNA>0.1*rand)

  # calculate median and 95% quantile for selected variable for each sample (BIN)
  r.m.sum <- r.m.sel %>%
    dplyr::mutate(
      sample.id = RUN.BIN,
      shift = RUN.SHIFT,
      !!sel.var := dplyr::select(.,-c(RUN.BIN,RUN.SHIFT)) %>%
        apply(.,1, FUN = function(x) median(x, na.rm = T)),
      !!paste0(sel.var,".05q") := dplyr::select(.,-c(RUN.BIN,RUN.SHIFT))
      %>% apply(.,1, FUN = function(x) quantile(x,0.025, na.rm = T)),
      !!paste0(sel.var,".95q") := dplyr::select(.,-c(RUN.BIN,RUN.SHIFT))
      %>% apply(.,1, FUN = function(x) quantile(x,0.975, na.rm = T))
    ) %>%
    dplyr::select(sample.id, shift, sel.var, paste0(sel.var,".05q"), paste0(sel.var,".95q"))

  return(r.m.sum)
}
