#' @keywords internal
#' @noRd
util_paste_as_vector <- function(var_list, sep = "'") {
  paste(
    paste0(sep, var_list, sep),
    collapse = ", "
  )
}

#' @keywords internal
#' @noRd
util_check_class <- function(x, sel_class) {
  var_name <- deparse(substitute(x))
  assertthat::assert_that(
    any(class(x) %in% sel_class),
    msg = paste0(
      "'", var_name, "' must be one of the following: ",
      util_paste_as_vector(sel_class)
    )
  )
}

#' @keywords internal
#' @noRd
util_check_col_names <- function(x, var_list) {
  var_name <- deparse(substitute(x))
  assertthat::assert_that(
    all(var_list %in% names(x)),
    msg = paste0(
      "'", var_name,
      "' must contains following columns: ",
      util_paste_as_vector(var_list)
    )
  )
}

#' @keywords internal
#' @noRd
util_check_vector_values <- function(x, var_list) {
  var_name <- deparse(substitute(x))
  assertthat::assert_that(
    any(var_list %in% x),
    msg = paste0(
      "'", var_name,
      "' must contains one of the following values: ",
      util_paste_as_vector(var_list)
    )
  )
}

#' @keywords internal
#' @noRd
util_check_if_integer <- function(x) {
  var_name <- deparse(substitute(x))
  assertthat::assert_that(
    round(x) == x,
    msg = paste0("'", var_name, "' must be an integer")
  )
}

#' @keywords internal
#' @noRd
util_output_comment <- function(msg = "") {
  message(msg)
}

#' @keywords internal
#' @noRd
util_output_heading <- function(msg = "", size = c("h1", "h2", "h3")) {
  size <- match.arg(size)
  sep_line <-
    switch(
      size,
      h1 = "#----------------------------------------------------------#",
      h2 = "#--------------------------------------------------#",
      h3 = "#----------------------------------------#"
    )
  message(sep_line)
  message(msg)
  message(sep_line)
}

#' @keywords internal
#' @noRd
util_output_warning <- function(msg = "") {
  warning(msg, call. = FALSE)
}

#' @title Flatten the list by one level
#' @param data_source A named list
#' @description
#' Returns a flat list whose names are the concatenation of the
#' parent and child names separated by `"-"`.
#' @keywords internal
#' @noRd
util_flatten_list_by_one <- function(data_source) {
  assertthat::assert_that(
    is.list(data_source),
    msg = "'data_source' must be a list"
  )
  names_high <- names(data_source)
  res <- vector("list", 0)
  for (i in seq_along(data_source)) {
    sel_item <- data_source[[i]]
    names(sel_item) <-
      paste(names_high[i], names(sel_item), sep = "-")
    res <- c(res, sel_item)
  }
  return(res)
}
