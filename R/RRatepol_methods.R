#' An S4 class for R- Ratepol data structure .

RRatepolList <-
  setClass(
    "RRatepolList",
    slots = c(
      Community = "data.frame",
      Age = "data.frame",
      Age.un = "data.frame",
      Dim.val = "integer"
    ))
