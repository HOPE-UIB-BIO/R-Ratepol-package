#' @title Subsetting levels in each Working units (WU)
#'
#' @param data_source_subset
#' Data.frame with community data
#' @param data_source_bins
#' Data.frame with individual WU to use
#' @inheritParams fc_estimate_RoC
#' @keywords internal
fc_subset_samples <-
  function(data_source_subset,
           data_source_bins,
           bin_selection = "first") {
    if (
      is.character(data_source_bins$start)
    ) {
      res <-
        data_source_bins %>%
        dplyr::select("label", "age_diff", "res_age", "start") %>%
        dplyr::inner_join(
          data_source_subset %>%
            tibble::rownames_to_column("start"),
          by = "start"
        ) %>%
        dplyr::mutate(
          age_diff = c(diff(.data$age), Inf),
          res_age = .data$age
        ) %>%
        dplyr::select(-c("start", "age"))

      return(res)
    }

    res_com <-
      as.data.frame(
        matrix(
          nrow = nrow(data_source_bins),
          ncol = ncol(data_source_subset) - 1,
          dimnames = list(
            NULL,
            names(data_source_subset)[2:ncol(data_source_subset)]
          )
        )
      )

    # for each bin
    for (i in 1:nrow(data_source_bins)) {

      # subset age data so it selected all samples which has higher values
      # than the BIN itself but
      # still small then selected bin + calculated BIN size
      subset_w <-
        data_source_subset[data_source_subset$age >= data_source_bins$start[i] &
          data_source_subset$age < data_source_bins$end[i], ]

      # If selected subset has at least one sample
      if (nrow(subset_w) > 0) {
        if (
          bin_selection == "random"
        ) {
          # select one random sample from the bin
          random_row <-
            sample(1:nrow(subset_w), 1)

          res_com[i, ] <-
            subset_w[random_row, -1]
        }

        if (
          bin_selection == "first"
        ) {
          # select the sample which is the closest to the beggining of the bin
          res_com[i, ] <-
            subset_w[1, -1]
        }
      }
    }

    res <-
      dplyr::bind_cols(
        data_source_bins %>%
          dplyr::select("label", "age_diff", "res_age"),
        res_com
      )

    return(res)
  }
