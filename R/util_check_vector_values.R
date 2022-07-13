#' @title Check the presence of variables in vector
#' @param data_source Name of the vector to check
#' @param var_list Vector with names
#' @description Function will test the presence of ANY values from `var_list`
#' within the `data_source` vector and return error message
#' @export
util_check_vector_values <-
  function(data_source, var_list) {
    parent_frame <- sys.parent()

    parent_env <- sys.frame(which = parent_frame)

    data_source_obj <- get(data_source, envir = parent_env)

    util_check_class("data_source_obj", "character")

    util_check_class("var_list", "character")

    assertthat::assert_that(
      any(var_list %in% data_source_obj),
      msg = paste0(
        "'", data_source, "' must contains one of the following values: ",
        util_paste_as_vector(var_list)
      )
    )
  }
