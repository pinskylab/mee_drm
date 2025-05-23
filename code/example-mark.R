library(dplyr)
library(ggplot2)
library(bayesplot)
library(sf)
library(drmr)
library(patchwork)
library(cmdstanr)
library(mgcv)

## loading data
data(sum_fl)

## loading map
map_name <- system.file("maps/sum_fl.shp", package = "drmr")

polygons <- st_read(map_name)

polygons |>
  st_area() |>
  units::set_units("km^2") |>
  summary()

##--- splitting data for validation ----

## reserving 5 years for forecast assessment
first_year_forecast <- max(sum_fl$year) - 4

## "year to id"
first_id_forecast <-
  first_year_forecast - min(sum_fl$year) + 1

years_all <- order(unique(sum_fl$year))
years_train <- years_all[years_all < first_id_forecast]
years_test <- years_all[years_all >= first_id_forecast]

dat_test <- sum_fl |>
  filter(year >= first_year_forecast)

dat_train <- sum_fl |>
  filter(year < first_year_forecast)

##--- centering covariates (for improved mcmc efficiency) ---

avgs <- c("stemp" = mean(dat_train$stemp),
          "btemp" = mean(dat_train$btemp),
          "depth" = mean(dat_train$depth),
          "n_hauls" = mean(dat_train$n_hauls),
          "lat" = mean(dat_train$lat),
          "lon" = mean(dat_train$lon))

min_year <- dat_train$year |>
  min()

## centering covariates
dat_train <- dat_train |>
  mutate(c_stemp = stemp - avgs["stemp"],
         c_btemp = btemp - avgs["btemp"],
         c_hauls = n_hauls - avgs["n_hauls"],
         ## depth = depth - avgs["depth"],
         c_lat   = lat - avgs["lat"],
         c_lon   = lon - avgs["lon"],
         time  = year - min_year)

dat_test <- dat_test |>
  mutate(c_stemp = stemp - avgs["stemp"],
         c_btemp = btemp - avgs["btemp"],
         c_hauls = n_hauls - avgs["n_hauls"],
         ## depth = depth - avgs["depth"],
         c_lat   = lat - avgs["lat"],
         c_lon   = lon - avgs["lon"],
         time  = year - min_year)

##--- turning response into density: 1k individuals per km2 ----

dat_train <- dat_train |>
  mutate(dens = 1000 * y / area_km2,
         .before = y)

dat_test <- dat_test |>
  mutate(dens = 1000 * y / area_km2,
         .before = y)

chains <- 4
cores <- 4

##--- fitting DRMs ----

my_drm <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "patch",
          n_ages = 12,
          seed = 202505,
          formula_zero = ~ 1 + c_hauls + c_btemp + c_stemp,
          formula_rec = ~ 1 + c_stemp + I(c_stemp * c_stemp),
          formula_surv = ~ 1,
          .toggles = list(est_surv = 1,
                          time_ar = 1),
          init = "pathfinder")

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
  select(year, patch, lat_floor, dens) |>
  mutate(pair = row_number(), .before = 1) |>
  left_join(forecasts_summary, by = "pair") |>
  select(- pair)

out_fitted <-
  dat_train |>
  select(year, patch, lat_floor, dens) |>
  bind_cols(select(fitted_summary, -pair))

ggplot(data = out_fitted,
       aes(x = year)) +
  geom_ribbon(aes(ymin = l, ymax = u),
              fill = 2,
              alpha = .4) +
  geom_line(aes(y = m)) +
  geom_point(aes(y = dens),
             color = 4) +
  facet_wrap(~ rev(patch),
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
  facet_wrap(~ rev(patch),
             scales = "free_y") +
  scale_y_continuous(breaks = scales::trans_breaks(identity, identity,
                                                   n = 3)) +
  theme_bw()

count(out_fitted, patch)
count(out_forecast, patch)

count(out_fitted, year)
count(out_forecast, year)

bind_rows(out_fitted, out_forecast) |>
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
  facet_wrap(~ patch) +
  scale_y_continuous(breaks = scales::trans_breaks(identity, identity,
                                                   n = 3),
                     trans = "log1p") +
  theme_bw()

##--- LN SDM ----

my_sdm <-
  fit_sdm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "patch",
          family = "lognormal",
          seed = 202505,
          formula_zero = ~ 1 + c_hauls + c_btemp + c_stemp,
          formula_dens = ~ 1 + c_stemp + I(c_stemp * c_stemp),
          init = "pathfinder")

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
  select(year, patch, lat_floor, dens) |>
  mutate(pair = row_number(), .before = 1) |>
  left_join(for_sdm, by = "pair") |>
  select(- pair)

sdm_fitted <-
  dat_train |>
  select(year, patch, lat_floor, dens) |>
  bind_cols(select(fitted_sdm, -pair))

ggplot(data = sdm_fitted,
       aes(x = year)) +
  geom_ribbon(aes(ymin = l, ymax = u),
              fill = 2,
              alpha = .4) +
  geom_line(aes(y = m)) +
  geom_point(aes(y = dens),
             color = 4) +
  facet_wrap(~ rev(patch),
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
  facet_wrap(~ rev(patch),
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
  facet_wrap(~ patch) +
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
  facet_wrap(model ~ patch) +
  scale_y_continuous(breaks = scales::trans_breaks(identity, identity,
                                                   n = 3),
                     trans = "log1p") +
  theme_bw()

##--- metric of comparison ----

out_forecast |>
  mutate(model = "DRM", .before = "year") |>
  bind_rows(mutate(sdm_forecast, model = "SDM", .before = "year")) |>
  ## filter(year != first_year_forecast) |>
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

flounder_out <-
  bind_rows(sdm_fitted, sdm_forecast) |>
  mutate(model = "SDM") |>
  bind_rows(
      bind_rows(out_fitted, out_forecast) |>
      mutate(model = "DRM")
  )

arrow::write_dataset(flounder_out,
                     path = "~/git-projects/lcgodoy.github.io/slides/2025-frcheck/data/sumf.parquet")

##--- taking a look at the Z ----

zt_summary <-
  my_drm$stanfit$draws(variables = "z_t",
                       format = "draws_df") |>
  tidyr::pivot_longer(cols = starts_with("z_t"),
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

ztp_summary <-
  forecast_drm$draws(variables = "z_tp",
                     format = "draws_df") |>
  tidyr::pivot_longer(cols = starts_with("z_tp"),
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

ztp_summary <- mutate(ztp_summary, pair = pair + max(zt_summary$pair))

zt_summary <- zt_summary |>
  bind_rows(ztp_summary)

ggplot(data = zt_summary,
       aes(x = pair,
           y = m)) +
  geom_ribbon(aes(ymin = l, ymax = u)) +
  geom_line(color = "white", lwd = 2) +
  geom_vline(xintercept = 31, color = 2, lty = 2,
             lwd = 2) +
  geom_hline(yintercept = 0, color = 2, lty = 3,
             lwd = 2) +
  theme_bw()

##--- taking a look at lambda ----

lbd_drm <- drmr:::lambda_drm(my_drm, cores = 1)

## age, time, patch

lbd_summary <-
  lbd_drm$draws(variables = "lambda",
                format = "draws_df") |>
  tidyr::pivot_longer(cols = starts_with("lambda"),
                      names_to = "pair",
                      values_to = "expected") |>
  group_by(pair) |>
  summarise(ll = quantile(expected, probs = .05),
            l = quantile(expected, probs = .1),
            m = median(expected),
            u = quantile(expected, probs = .9),
            uu = quantile(expected, probs = .95)) |>
  ungroup() |>
  mutate(ages = gsub("[a-z]", "", pair)) |>
  mutate(patch = sapply(strsplit(ages, split = ","),
                        \(x) gsub("\\D", "", x[3]))) |>
  mutate(time  = sapply(strsplit(ages, split = ","),
                        \(x) gsub("\\D", "", x[2]))) |>
  mutate(age   = sapply(strsplit(ages, split = ","),
                        \(x) gsub("\\D", "", x[1]))) |>
  mutate(across(patch:age, as.integer)) |>
  select(-ages)

lbd_fct <-
  forecast_drm$draws(variables = "lambda_proj",
                     format = "draws_df") |>
  tidyr::pivot_longer(cols = starts_with("lambda_proj"),
                      names_to = "pair",
                      values_to = "expected") |>
  group_by(pair) |>
  summarise(ll = quantile(expected, probs = .05),
            l = quantile(expected, probs = .1),
            m = median(expected),
            u = quantile(expected, probs = .9),
            uu = quantile(expected, probs = .95)) |>
  ungroup() |>
  mutate(ages = gsub("[a-z]", "", pair)) |>
  mutate(patch = sapply(strsplit(ages, split = ","),
                        \(x) gsub("\\D", "", x[3]))) |>
  mutate(time  = sapply(strsplit(ages, split = ","),
                        \(x) gsub("\\D", "", x[2]))) |>
  mutate(age   = sapply(strsplit(ages, split = ","),
                        \(x) gsub("\\D", "", x[1]))) |>
  mutate(across(patch:age, as.integer)) |>
  select(-ages)

lbd_fct <- lbd_fct |>
  mutate(time = time + max(lbd_summary$time))

lbd_all <- bind_rows(lbd_summary, lbd_fct)

lbd_all |>
  filter(time %in% c(29:max(time))) |>
  ggplot(data = _,
         aes(x = age)) +
  geom_segment(aes(yend = m, y = 0)) +
  geom_point(aes(y = m)) +
  facet_grid(patch ~ time) +
  theme_bw()
