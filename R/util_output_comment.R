#' @title Output message in console
#' @param msg String with the message that should be printed
#' @export
util_output_comment <-
  function(msg = "") {
    assertthat::assert_that(
      assertthat::is.string(msg),
      msg = "'msg' must be a 'string'"
    )

    cat("\n")
    usethis::ui_info(msg)
    cat("\n")
  }
