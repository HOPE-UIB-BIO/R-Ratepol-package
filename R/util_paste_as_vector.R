#' @title Paste any list of values as a vector
#' @param var_list List of values
#' @param sep Separator which should be used to wrap each object in
#' the `var_list`
#' @export
util_paste_as_vector <-
  function(var_list, sep = "'") {
    paste(
      paste0(sep, var_list, sep),
      collapse = ", "
    ) %>%
      return()
  }
