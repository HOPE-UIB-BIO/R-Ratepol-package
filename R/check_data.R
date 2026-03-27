#' @title Check the data
#' @description
#' Output summary information about the input data to the console.
#' @param data_source_check
#' `list` with `community` (`data.frame`) and `age` (`data.frame`)
#' elements.
#' @param silent
#' `logical`. If `TRUE`, suppress all console output.
#' @return
#' `NULL` invisibly.
#' @keywords internal
check_data <- function(data_source_check, silent = FALSE) {
  util_check_class(data_source_check, "list")

  assertthat::assert_that(
    nrow(data_source_check$community) > 0,
    ncol(data_source_check$community) >= 1,
    nrow(data_source_check$age) > 0,
    ncol(data_source_check$age) == 1,
    msg = "Error: Object 'data_source_check' was supplied with empty elements: 'community' and 'age'"
  )
  ## 1. dimensions
  dims <-
    lapply(data_source_check, dim)

  if (isFALSE(silent)) {
    if (length(dims) == 3) {
      util_output_comment(
        paste(
          "Community data have",
          dims$community[2],
          "taxa and",
          dims$community[1],
          "samples.",
          "Age data have",
          dims$age[1],
          "samples.",
          "Age uncertainty data have",
          dims$age_un[2],
          "samples."
        )
      )
    } else if (length(dims) == 2) {
      util_output_comment(
        paste(
          "Community data have",
          dims$community[2],
          "taxa and",
          dims$community[1],
          "samples.",
          "Age data have",
          dims$age[1],
          "samples.",
          "Age uncertainty was not provided."
        )
      )
    }
  }

  # 2. Missing data
  nas <-
    lapply(data_source_check, function(x) {
      sum(is.na(x))
    })

  if (isFALSE(silent)) {
    if (length(nas) == 3) {
      util_output_comment(
        paste(
          "Community data has",
          nas$community,
          "NAs.",
          "Age data has",
          nas$age,
          "NAs.",
          "Age uncertainty data has",
          nas$age_un,
          "NAs."
        )
      )
    } else if (length(nas) == 2) {
      util_output_comment(
        paste(
          "Community data has",
          nas$community,
          "NAs.",
          "Age data has",
          nas$age,
          "NAs.",
          "Age uncertainty was not provided."
        )
      )
    }

    # 3. Summary of data
    util_output_comment(
      paste0(
        "Community data has a value of min ",
        round(min(rowSums(data_source_check$community, na.rm = TRUE))),
        ", max ",
        round(max(rowSums(data_source_check$community, na.rm = TRUE))),
        ", mean ",
        round(mean(rowSums(data_source_check$community, na.rm = TRUE))),
        ", and median ",
        round(stats::median(rowSums(
          data_source_check$community,
          na.rm = TRUE
        ))),
        " observations."
      )
    )

    util_output_comment(
      paste0(
        "Age data has values of min ",
        round(min(data_source_check$age$age, na.rm = TRUE)),
        ", max ",
        round(max(data_source_check$age$age, na.rm = TRUE)),
        ", mean ",
        round(mean(data_source_check$age$age, na.rm = TRUE)),
        ", and median ",
        round(stats::median(data_source_check$age$age, na.rm = TRUE))
      )
    )
  }
}
