data_source_smooth <-
  extract_data(
    RRatepol::example_data$pollen_data[[1]],
    RRatepol::example_data$sample_age[[1]]
  )

result <-
  smooth_community_data(
    data_source_smooth,
    smooth_method = c("m.avg", "grim", "age.w", "shep"),
    smooth_n_points = 5,
    smooth_n_max = 9,
    smooth_age_range = 500,
    round_results = FALSE,
    verbose = FALSE
  )

str(result)
