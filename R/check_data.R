#' @title Check the data
#'
#' @param data_source_check
#' List with `community` and `age`
#' @description Output summary information about the data
#' @keywords internal
check_data <-
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
        "Age data has values of min ", round(min(data_source_check$age$age, na.rm = TRUE)),
        ", max ", round(max(data_source_check$age$age, na.rm = TRUE)),
        ", mean ", round(mean(data_source_check$age$age, na.rm = TRUE)),
        ", median ", round(stats::median(data_source_check$age$age, na.rm = TRUE)),
        ", and NAs ", sum(is.na(data_source_check$age$age))
      )
    )
  }
