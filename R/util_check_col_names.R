#' @title Check the presence of variables in data.frame
#' @param data_source Name of the data.frame to check
#' @param var_list Vector with names
#' @description Function will test the presence of ALL names in `var_list`
#' within the `data_source` data.frame and return error message
#' @export
util_check_col_names <-
  function(data_source, var_list) {
    parent_frame <- sys.parent()

    parent_env <- sys.frame(which = parent_frame)

    data_source_obj <- get(data_source, envir = parent_env)

    util_check_class("data_source_obj", "data.frame")

    util_check_class("var_list", "character")

    assertthat::assert_that(
      all(var_list %in% names(data_source_obj)),
      msg = paste0(
        "'", data_source, "' must contains following columns: ",
        util_paste_as_vector(var_list)
      )
    )
  }
