#' @title Prepare data for single RoC estimation
#'
#' @inheritParams estimate_roc
#' @param data_source_prep
#' List with `community` and `age`
#' @description
#' Create a single list with all information needed to estimate RoC.
#' This list can then be evaluated in parallel across randomisation runs.
#' @return
#' A `list` with preprocessed `data` and `bins` elements.
#' @keywords internal
prepare_data <- function(
  data_source_prep,
  working_units = c("levels", "bins", "MW"),
  bin_size = 500,
  number_of_shifts = 5,
  rand = NULL
) {
  util_check_class(data_source_prep, "list")

  assertthat::assert_that(
    !base::is.null(data_source_prep$community) &&
      base::nrow(data_source_prep$community) > 0,
    msg = "'data_source_prep$community' must not be empty"
  )

  assertthat::assert_that(
    !base::is.null(data_source_prep$age) &&
      base::nrow(data_source_prep$age) > 0,
    msg = "'data_source_prep$age' must not be empty"
  )

  util_check_class(working_units, "character")

  util_check_vector_values(working_units, c("levels", "bins", "MW"))

  working_units <- match.arg(working_units)

  if (working_units != "levels") {
    util_check_class(bin_size, "numeric")

    util_check_if_integer(bin_size)
  }

  util_check_class(rand, c("NULL", "numeric"))

  # check the condition
  is_shift_present <-
    working_units == "MW" && (number_of_shifts != 0)

  if (isFALSE(is_shift_present)) {
    number_of_shifts <- 1
  } else {
    util_check_class(number_of_shifts, "numeric")
    util_check_if_integer(number_of_shifts)
  }

  is_rand_present <-
    isFALSE(is.null(rand))

  if (isTRUE(is_rand_present)) {
    util_check_if_integer(rand)
  } else {
    rand <- 1
  }

  if (working_units == "levels") {
    bin_dummy <-
      make_bins(
        data_source_bins = data_source_prep,
        working_units = "levels"
      )
  } else if (isTRUE(is_shift_present)) {
    bin_dummy <-
      make_bins(
        data_source_bins = data_source_prep,
        working_units = "MW",
        bin_size = bin_size,
        number_of_shifts = number_of_shifts
      )
  } else {
    bin_dummy <-
      make_bins(
        data_source_bins = data_source_prep,
        working_units = "bins",
        bin_size = bin_size
      )
  }

  is_uncertit_present <-
    isFALSE(is.null(data_source_prep$age_un))

  if (isTRUE(is_uncertit_present)) {
    if (isTRUE(is_rand_present)) {
      random_value <-
        sample(
          c(
            1:max(1, nrow(data_source_prep$age_un))
          ),
          rand,
          replace = TRUE
        )

      random_age <-
        purrr::map(
          .x = seq_along(random_value),
          .f = ~ as.numeric(
            data_source_prep$age_un[random_value[.x], ]
          )
        )
    } else {
      random_age <-
        list(data_source_prep$age$age)
    }
  }

  data_merge <-
    data_source_prep$community %>%
    tibble::rownames_to_column("row_name") %>%
    dplyr::left_join(
      data_source_prep$age %>%
        tibble::rownames_to_column("row_name"),
      by = "row_name"
    ) %>%
    dplyr::relocate("age") %>%
    tibble::column_to_rownames("row_name")

  shift_vec <-
    c(1:number_of_shifts) %>%
    rlang::set_names(
      nm = as.character(1:number_of_shifts)
    )

  rand_vec <-
    c(1:rand) %>%
    rlang::set_names(
      nm = as.character(1:rand)
    )

  if (
    # is_rand_present is TRUE
    isFALSE(is_uncertit_present)
  ) {
    rand_vec %>%
      purrr::map(
        .f = ~ purrr::map(
          .x = shift_vec,
          .f = ~ list(
            data = data_merge,
            bins = bin_dummy %>%
              dplyr::filter(shift == .x)
          )
        )
      ) %>%
      return()
  } else {
    # is_uncertit_present is TRUE
    rand_vec %>%
      purrr::map(
        .f = ~ {
          data_with_random_age <-
            data_merge %>%
            dplyr::mutate(
              age = random_age[[.x]]
            )

          purrr::map(
            .x = shift_vec,
            .f = ~ list(
              data = data_with_random_age,
              bins = bin_dummy %>%
                dplyr::filter(shift == .x)
            ) %>%
              return()
          ) %>%
            return()
        }
      ) %>%
      return()
  }
  # stop("Not a valid setting for creating bins")
}
