
#' @title Check the data
#'
#' @param data_source_check
#' List with `community` and `age`
#' @description Output summary information about the data
fc_check_data <-
  function(data_source_check) {
    RUtilpol::check_class("data_source_check", "list")

    RUtilpol::output_comment(
      paste(
        "Community data have", ncol(data_source_check$community),
        "taxa and", nrow(data_source_check$community), "samples.",
        " Age data have", nrow(data_source_check$age), "samples"
      )
    )

    RUtilpol::output_comment(
      paste0(
        "Age data has values of min ", round(min(data_source_check$age$age)),
        ", max ", round(max(data_source_check$age$age)),
        ", mean ", round(mean(data_source_check$age$age)),
        ", and median ", round(stats::median(data_source_check$age$age))
      )
    )
  }
