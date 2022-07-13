#' @title Output boxed message in console
#' @param msg String with the message that should be printed
#' @return NULL
#' @export
util_output_message <-
  function(msg = "") {
    assertthat::assert_that(
      assertthat::is.string(msg),
      msg = "'msg' must be a 'string'"
    )

    sep_line <-
      "#----------------------------------------------------------#"

    cat("\n")
    usethis::ui_line(sep_line)
    util_output_comment(msg = msg)
    usethis::ui_line(sep_line)
    cat("\n")
  }
