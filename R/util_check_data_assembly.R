#' @title Check the data assembly
#' @param data_source Data.frame to be tested
#' @description Check the dataset for mandatory and non-mandatory
#' variables. Stop the session if mandatory not present. Warn about
#' non-mandatory and fill them with 'NA'
#' @return Data.frame with all variables
util_check_data_assembly <-
  function(data_source) {
    util_check_class("data_source", "data.frame")

    mandatory_vars <-
      c(
        "dataset_id",
        "long",
        "lat",
        "altitude",
        "depositionalenvironment",
        "n_sample_counts",
        "sample_depth",
        "raw_counts",
        "chron_control",
        "n_chron_control"
      )

    util_check_col_names("data_source", c(mandatory_vars))

    non_mandatory_vars <-
      c(
        "handle",
        "siteid",
        "sitename",
        "source_of_data",
        "pollen_percentage",
        "data_publicity"
      )

    non_mandatory_vars_present <-
      names(data_source) %>%
      subset(., !names(data_source) %in% mandatory_vars)

    missing_vars <-
      subset(
        non_mandatory_vars,
        !non_mandatory_vars %in% non_mandatory_vars_present
      )


    if (
      length(missing_vars) > 0
    ) {
      util_output_comment(
        msg = paste(
          "The following variables have been added to dataset with NAs:",
          paste(
            paste0("'", missing_vars, "'"),
            collapse = ", "
          )
        )
      )

      data_cor <-
        data_source %>%
        dplyr::bind_cols(
          missing_vars %>%
            purrr::set_names() %>%
            purrr::map_dfc(
              .x = .,
              .f = ~ rep(NA, nrow(data_source))
            )
        )
    } else {
      util_output_comment(
        msg = "All variables present"
      )

      data_cor <- data_source
    }

    # set correct data classes
    data_res <-
      data_cor %>%
      dplyr::mutate(
        dataset_id = as.character(dataset_id),
        handle = as.character(handle),
        siteid = as.character(siteid),
        sitename = as.character(sitename),
        long = as.double(long),
        lat = as.double(lat),
        altitude = as.double(altitude),
        source_of_data = as.character(source_of_data),
        data_publicity = as.character(data_publicity),
        depositionalenvironment = as.character(depositionalenvironment),
        n_sample_counts = as.integer(n_sample_counts),
        sample_depth = purrr::map(sample_depth, ~ tibble::as_tibble(.x)),
        raw_counts = purrr::map(raw_counts, ~ tibble::as_tibble(.x)),
        pollen_percentage = as.logical(pollen_percentage),
        n_chron_control = as.integer(n_chron_control),
        chron_control = purrr::map(chron_control, ~ tibble::as_tibble(.x))
      ) %>%
      dplyr::relocate(
        dataset_id, handle, siteid, sitename, long, lat, altitude,
        source_of_data, data_publicity, depositionalenvironment,
        n_sample_counts, sample_depth, raw_counts, pollen_percentage,
        n_chron_control, chron_control
      )

    util_check_col_names(
      "data_res",
      c(
        "dataset_id",
        "handle",
        "siteid",
        "sitename",
        "long",
        "lat",
        "altitude",
        "source_of_data",
        "data_publicity",
        "depositionalenvironment",
        "n_sample_counts",
        "sample_depth",
        "raw_counts",
        "pollen_percentage",
        "n_chron_control",
        "chron_control"
      )
    )


    return(data_res)
  }
