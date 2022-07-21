#' @title Output warning in console
#' @param msg String with the message that should be printed
util_output_warning <-
  function(msg = "") {
    assertthat::assert_that(
      assertthat::is.string(msg),
      msg = "'msg' must be a 'string'"
    )

    cat("\n")
    usethis::ui_oops(msg)
    cat("\n")
  }
