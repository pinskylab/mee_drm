library(dplyr)
library(ggplot2)
library(bayesplot)
library(sf)
library(drmr)
library(patchwork)
library(cmdstanr)
library(arrow)
library(geoarrow)

## loading data
my_dt <- open_dataset("data/birds/processed.parquet") |>
  st_as_sf()

my_dt <- my_dt |>
  mutate(id = as.integer(factor(id)),
         lon = st_coordinates(st_centroid(geometry))[, 1],
         lat = st_coordinates(st_centroid(geometry))[, 2]) |>
  arrange(id, year)

polygons <- my_dt |>
  filter(year == min(year)) |>
  st_geometry()

polygons |>
  st_area() |>
  units::set_units("km^2") |>
  summary()

my_dt <- st_drop_geometry(my_dt)

##--- splitting data for validation ----

## reserving 5 years for forecast assessment
first_year_forecast <- max(my_dt$year) - 4

## "year to id"
first_id_forecast <-
  first_year_forecast - min(my_dt$year) + 1

years_all <- order(unique(my_dt$year))
years_train <- years_all[years_all < first_id_forecast]
years_test <- years_all[years_all >= first_id_forecast]

dat_test <- my_dt |>
  filter(year >= first_year_forecast)

dat_train <- my_dt |>
  filter(year < first_year_forecast)

##--- centering covariates (for improved mcmc efficiency) ---

avgs <- c("tavg" = mean(dat_train$tavg),
          "lon" = mean(dat_train$lon),
          "lat" = mean(dat_train$lat))

min_year <- dat_train$year |>
  min()

## centering covariates
dat_train <- dat_train |>
  mutate(c_tavg = tavg - avgs["tavg"],
         c_lat  = lat - avgs["lat"],
         c_lon  = lon - avgs["lon"],
         time   = year - min_year)

dat_test <- dat_test |>
  mutate(c_tavg = tavg - avgs["tavg"],
         c_lat  = lat - avgs["lat"],
         c_lon  = lon - avgs["lon"],
         time   = year - min_year)

##--- turning response into density: 1k individuals per km2 ----

dat_train <- dat_train |>
  mutate(dens = 1000 * y,
         .before = y)

dat_test <- dat_test |>
  mutate(dens = 1000 * y,
         .before = y)

chains <- 4
cores <- 4

##--- fitting SDMs ----

## fix pp_sim for SDM

my_sdm <-
  fit_sdm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "id",
          seed = 202505,
          family = "lognormal",
          formula_zero = ~ 1 + n_routes +
            c_tavg + c_lon + c_lat,
          formula_dens = ~ 1 + c_tavg + I(c_tavg * c_tavg),
          .toggles = list(time_ar = 0),
          .priors = list(pr_phi_a = .1, pr_phi_b = .1),
          init = "pathfinder")


sdm$stanfit$summary(variables = c("beta_t", "beta_r",
                                  "tau", "alpha", "phi"))

##--- fitting DRMs ----

my_drm <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "id",
          seed = 202505,
          family = "loglogistic",
          formula_zero = ~ 1 + n_routes +
            c_tavg + c_lon + c_lat,
          formula_rec = ~ 1 + c_tavg + I(c_tavg * c_tavg),
          formula_surv = ~ 1,
          .toggles = list(est_surv = 1,
                          time_ar = 1))


##--- forecasting ----

forecast_drm <- predict_drm(drm = my_drm,
                            new_data = dat_test,
                            past_data = filter(dat_train,
                                               year == max(year)),
                            seed = 125,
                            cores = 4)

##--- * predicted vs observed time series ----

fitted_summary <-
  my_drm$stanfit$draws(variables = "y_pp",
                       format = "draws_df") |>
  tidyr::pivot_longer(cols = starts_with("y_pp"),
                      names_to = "pair",
                      values_to = "expected") |>
  group_by(pair) |>
  summarise(ll = quantile(expected, probs = .05),
            l = quantile(expected, probs = .1),
            m = median(expected),
            u = quantile(expected, probs = .9),
            uu = quantile(expected, probs = .95)) |>
  ungroup() |>
  mutate(pair = gsub("\\D", "", pair)) |>
  mutate(pair = as.integer(pair)) |>
  arrange(pair)

forecasts_summary <-
  forecast_drm$draws(variables = "y_proj",
                     format = "draws_df") |>
  tidyr::pivot_longer(cols = starts_with("y_proj"),
                      names_to = "pair",
                      values_to = "expected") |>
  group_by(pair) |>
  summarise(ll = quantile(expected, probs = .05),
            l = quantile(expected, probs = .1),
            m = median(expected),
            u = quantile(expected, probs = .9),
            uu = quantile(expected, probs = .95)) |>
  ungroup() |>
  mutate(pair = gsub("\\D", "", pair)) |>
  mutate(pair = as.integer(pair)) |>
  arrange(pair)

out_forecast <-
  dat_test |>
  select(year, id, lat, lon, dens) |>
  mutate(pair = row_number(), .before = 1) |>
  left_join(forecasts_summary, by = "pair") |>
  select(- pair)

out_fitted <-
  dat_train |>
  select(year, id, lat, lon, dens) |>
  bind_cols(select(fitted_summary, -pair))

ggplot(data = out_fitted,
       aes(x = year)) +
  geom_ribbon(aes(ymin = l, ymax = u),
              fill = 2,
              alpha = .4) +
  geom_line(aes(y = m)) +
  geom_point(aes(y = dens),
             color = 4) +
  facet_wrap(~ id,
             scales = "free_y") +
  scale_y_continuous(breaks = scales::trans_breaks(identity, identity,
                                                   n = 3)) +
  theme_bw()

ggplot(data = out_forecast,
       aes(x = year)) +
  geom_ribbon(aes(ymin = l, ymax = u),
              fill = 2,
              alpha = .4) +
  geom_line(aes(y = m)) +
  geom_point(aes(y = dens),
             color = 4) +
  facet_wrap(~ id,
             scales = "free_y") +
  scale_y_continuous(breaks = scales::trans_breaks(identity, identity,
                                                   n = 3)) +
  theme_bw()

set.seed(125)
ids <- sample(seq_len(max(out_fitted$id)),
              size = 5)

bind_rows(out_fitted, out_forecast) |>
  filter(id %in% ids) |>
  filter(year != first_year_forecast) |>
  ggplot(data = _) +
  geom_vline(xintercept = first_year_forecast,
             lty = 2) +
  geom_ribbon(aes(x = year,
                  ymin = l, ymax = u),
              fill = 2,
              alpha = .4) +
  geom_line(aes(x = year, y = m)) +
  geom_point(aes(x = year, y = dens),
             color = 4) +
  facet_wrap(~ id) +
  scale_y_continuous(breaks = scales::trans_breaks(identity, identity,
                                                   n = 3),
                     trans = "log1p") +
  theme_bw()

##--- forcasts SDM ----

forecast_sdm <-
  predict_sdm(sdm = my_sdm,
              new_data = dat_test,
              seed = 125,
              cores = 4)

##--- * predicted vs observed time series ----

fitted_sdm <-
  my_sdm$stanfit$draws(variables = "y_pp",
                       format = "draws_df") |>
  tidyr::pivot_longer(cols = starts_with("y_pp"),
                      names_to = "pair",
                      values_to = "expected") |>
  group_by(pair) |>
  summarise(ll = quantile(expected, probs = .05),
            l = quantile(expected, probs = .1),
            m = median(expected),
            u = quantile(expected, probs = .9),
            uu = quantile(expected, probs = .95)) |>
  ungroup() |>
  mutate(pair = gsub("\\D", "", pair)) |>
  mutate(pair = as.integer(pair)) |>
  arrange(pair)

for_sdm <-
  forecast_sdm$draws(variables = "y_proj",
                     format = "draws_df") |>
  tidyr::pivot_longer(cols = starts_with("y_proj"),
                      names_to = "pair",
                      values_to = "expected") |>
  group_by(pair) |>
  summarise(ll = quantile(expected, probs = .05),
            l = quantile(expected, probs = .1),
            m = median(expected),
            u = quantile(expected, probs = .9),
            uu = quantile(expected, probs = .95)) |>
  ungroup() |>
  mutate(pair = gsub("\\D", "", pair)) |>
  mutate(pair = as.integer(pair)) |>
  arrange(pair)

sdm_forecast <-
  dat_test |>
  select(year, id, lat, lon, dens) |>
  mutate(pair = row_number(), .before = 1) |>
  left_join(for_sdm, by = "pair") |>
  select(- pair)

sdm_fitted <-
  dat_train |>
  select(year, id, lat, lon, dens) |>
  bind_cols(select(fitted_sdm, -pair))

ggplot(data = sdm_fitted,
       aes(x = year)) +
  geom_ribbon(aes(ymin = l, ymax = u),
              fill = 2,
              alpha = .4) +
  geom_line(aes(y = m)) +
  geom_point(aes(y = dens),
             color = 4) +
  facet_wrap(~ id,
             scales = "free_y") +
  scale_y_continuous(breaks = scales::trans_breaks(identity, identity,
                                                   n = 3)) +
  theme_bw()

ggplot(data = sdm_forecast,
       aes(x = year)) +
  geom_ribbon(aes(ymin = l, ymax = u),
              fill = 2,
              alpha = .4) +
  geom_line(aes(y = m)) +
  geom_point(aes(y = dens),
             color = 4) +
  facet_wrap(~ id,
             scales = "free_y") +
  scale_y_continuous(breaks = scales::trans_breaks(identity, identity,
                                                   n = 3)) +
  theme_bw()

bind_rows(sdm_fitted, sdm_forecast) |>
  ## filter(patch == 6) |>
  ggplot(data = _) +
  geom_vline(xintercept = first_year_forecast,
             lty = 2) +
  geom_ribbon(aes(x = year,
                  ymin = l, ymax = u),
              fill = 2,
              alpha = .4) +
  geom_line(aes(x = year, y = m)) +
  geom_point(aes(x = year, y = dens),
             color = 4) +
  facet_wrap(~ id) +
  scale_y_continuous(breaks = scales::trans_breaks(identity, identity,
                                                   n = 3),
                     trans = "log1p") +
  theme_bw()

##--- comparing both ----

bind_rows(sdm_fitted, sdm_forecast) |>
  mutate(model = "SDM") |>
  bind_rows(
      bind_rows(out_fitted, out_forecast) |>
      mutate(model = "DRM")
  ) |>
  filter(id %in% ids) |>
  filter(year != first_year_forecast) |>
  ggplot(data = _) +
  geom_vline(xintercept = first_year_forecast,
             lty = 2) +
  geom_ribbon(aes(x = year,
                  ymin = l, ymax = u,
                  fill = model,
                  color = model),
              alpha = .4) +
  geom_line(aes(x = year, y = m, color = model)) +
  geom_point(aes(x = year, y = dens)) +
  facet_wrap(model ~ id) +
  scale_y_continuous(breaks = scales::trans_breaks(identity, identity,
                                                   n = 3),
                     trans = "log1p") +
  theme_bw()

##--- metric of comparison ----

out_forecast |>
  mutate(model = "DRM", .before = "year") |>
  bind_rows(mutate(sdm_forecast, model = "SDM", .before = "year")) |>
  filter(year != first_year_forecast) |>
  mutate(bias = dens - m) |>
  mutate(rmse = bias * bias) |>
  mutate(is = int_score(dens, l = l, u = u, alpha = .2)) |>
  mutate(cvg = 100 * data.table::between(dens, l, u)) |>
  ungroup() |>
  group_by(model) |>
  summarise(across(rmse:cvg, mean)) |>
  ungroup() |>
  rename_all(toupper) |>
  rename("Model" = "MODEL",
         "IS (80%)" = "IS",
         "PIC (80%)" = "CVG") |>
  ## left_join(aux_qt,
  ##           by = "Model") |>
  ## arrange(desc(LOOIC)) |>
  ## relocate(LOOIC, .after = "Model") |>
  print() ## |>
  ## xtable::xtable(caption = "Forecasting skill according to different metrics",
  ##                digits = 2) |>
  ## print(include.rownames = FALSE)

wood_out <-
  bind_rows(sdm_fitted, sdm_forecast) |>
  mutate(model = "SDM") |>
  bind_rows(
      bind_rows(out_fitted, out_forecast) |>
      mutate(model = "DRM")
  )


wood_out <- open_dataset("data/birds/processed.parquet") |>
  st_as_sf() |>
  mutate(id = as.integer(factor(id)),
         lon = st_coordinates(st_centroid(geometry))[, 1],
         lat = st_coordinates(st_centroid(geometry))[, 2]) |>
  arrange(id, year) |>
  filter(year == min(year)) |>
  select(id) |>
  left_join(wood_out, by = "id")

wood_out

write_dataset(wood_out,
              path = "~/git-projects/lcgodoy.github.io/slides/2025-frcheck/data/wood.parquet")

wood_out |>
  filter(year %in% c((first_year_forecast + 1):max(year))) |>
  tidyr::pivot_longer(cols = c("dens", "m")) |>
  mutate(model = if_else(name == "dens", "Observed", model)) |>
  ggplot(data = _) +
  geom_sf(aes(fill = value)) +
  scale_fill_viridis_c(option = "H", trans = "log1p") +
  facet_grid(model ~ year) +
  theme_bw() +
  theme()
