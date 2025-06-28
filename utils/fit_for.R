fitted_drm <- function(drm, ci_mass) {
  tails <- 1 - ci_mass
  lower_tail <- tails * .5
  upper_tail <- 1 - .5 * tails
  drm$stanfit$draws(variables = "y_pp",
                    format = "draws_df") |>
    tidyr::pivot_longer(cols = starts_with("y_pp"),
                        names_to = "pair",
                        values_to = "expected") |>
    dplyr::group_by(pair) |>
    dplyr::summarise(l = quantile(expected, probs = lower_tail),
                     m = median(expected),
                     u = quantile(expected, probs = upper_tail)) |>
    dplyr::ungroup() |>
    dplyr::mutate(pair = gsub("\\D", "", pair)) |>
    dplyr::mutate(pair = as.integer(pair)) |>
    dplyr::arrange(pair)
}

forecast_drm <- function(for_drm, ci_mass) {
  tails <- 1 - ci_mass
  lower_tail <- tails * .5
  upper_tail <- 1 - .5 * tails
  for_drm$draws(variables = "y_proj",
                format = "draws_df") |>
    tidyr::pivot_longer(cols = starts_with("y_proj"),
                        names_to = "pair",
                        values_to = "expected") |>
    dplyr::group_by(pair) |>
    dplyr::summarise(l = quantile(expected, probs = lower_tail),
                     m = median(expected),
                     u = quantile(expected, probs = upper_tail)) |>
    dplyr::ungroup() |>
    dplyr::mutate(pair = gsub("\\D", "", pair)) |>
    dplyr::mutate(pair = as.integer(pair)) |>
    dplyr::arrange(pair)
}
