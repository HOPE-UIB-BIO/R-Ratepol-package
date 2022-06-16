#' @title An S4 class for R-Ratepol data structure
#'
#' @slot Community data.frame. 
#' @slot Age data.frame. 
#' @slot Age.un data.frame. 
#' @slot Dim.val integer. 
#'
#' @export
RRatepolList <-
  methods::setClass(
    "RRatepolList",
    slots = c(
      Community = "data.frame",
      Age = "data.frame",
      Age.un = "data.frame",
      Dim.val = "integer"
    ))
