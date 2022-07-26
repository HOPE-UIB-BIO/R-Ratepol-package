#' @title Assert the class of object
#' @param data_source The name of the object in quotes
#' @description The function will evaluate the object of `data_source` name
#' and test if ANY of the classes in equal to `sel_class`
#' @return `TRUE` if class match or error message
util_check_if_integer <-
  function(data_source) {
    parent_frame <- sys.parent()

    parent_env <- sys.frame(which = parent_frame)

    data_source_value <- get(data_source, envir = parent_env)

    assertthat::assert_that(
      round(data_source_value) == data_source_value,
      msg = paste(
        "'", data_source, "' must be a an integer"
      )
    )
  }
