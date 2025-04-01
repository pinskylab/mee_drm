library(dplyr)
library(ggplot2)
library(bayesplot)
library(sf)
library(drmr)
library(cmdstanr)

fix_intercept <- function(beta0, beta1, beta2, offset) {
  beta0 + offset * (beta2 - beta1)
}

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
first_year_forecast <- max(sum_fl$year) - 5

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

##--- turning response into density: individuals pear km2 per haul ----

dat_train <- dat_train |>
  mutate(dens = 1000 * y / area_km2,
         .before = y)

dat_test <- dat_test |>
  mutate(dens = 1000 * y / area_km2,
         .before = y)

chains <- 4
cores <- 4

##--- fitting SDM (with movement and F) ----

drm_0 <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "patch",
          seed = 202505)

drm_1 <-
  fit_drm(.data = dat_train,
          y_col = "dens",
          time_col = "year",
          site_col = "patch",
          seed = 202505,
          formula_zero = ~ 1 + c_hauls + c_btemp + c_stemp,
          formula_rec = ~ 1 + c_stemp + I(c_stemp * c_stemp),
          m = 0.25,
          init = "pathfinder")

drm_2 <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "patch",
          seed = 202505,
          formula_zero = ~ 1 + c_hauls + c_btemp + c_stemp,
          formula_rec = ~ 1 + c_stemp + I(c_stemp * c_stemp),
          formula_surv = ~ 1,
          .toggles = list(est_surv = 1),
          init = "pathfinder")

drm_3 <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "patch",
          seed = 202505,
          formula_zero = ~ 1 + c_hauls + c_btemp + c_stemp,
          formula_rec = ~ 1 + c_stemp + I(c_stemp * c_stemp),
          formula_surv = ~ 1,
          .toggles = list(est_surv = 1,
                          time_ar = 1),
          init = "pathfinder")

loos <- list("drm_0" = drm_0$stanfit$loo(),
             "drm_1" = drm_1$stanfit$loo(),
             "drm_2" = drm_2$stanfit$loo(),
             "drm_3" = drm_3$stanfit$loo())

loo::loo_compare(loos)

adj_mat <- gen_adj(st_buffer(st_geometry(polygons),
                             dist = 2500))

## row-standardized matrix
adj_mat <-
  t(apply(adj_mat, 1, \(x) x / (sum(x))))

drm_4 <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "patch",
          family = "gamma",
          seed = 202505,
          formula_zero = ~ 1 + c_hauls + c_btemp + c_stemp,
          formula_rec = ~ 1 + c_stemp + I(c_stemp * c_stemp),
          formula_surv = ~ 1,
          n_ages = 6,
          adj_mat = adj_mat, ## A matrix for movement routine
          ages_movement = c(0, 0, 1, 1, 1, 0), ## ages allowed to move
          .toggles = list(est_surv = 1,
                          time_ar = 1,
                          movement = 1),
          init = "pathfinder")

drm_4$stanfit$summary(variables = c("phi", "alpha", "zeta",
                                    "beta_r", "beta_s",
                                    "beta_t"))

loos <- c(loos, list("drm_4" = drm_4$stanfit$loo()))

loo::loo_compare(loos)

drm_5 <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "patch",
          family = "gamma",
          seed = 202505,
          formula_zero = ~ 1 + c_hauls + c_btemp + c_stemp,
          formula_rec = ~ 1 + c_stemp + I(c_stemp * c_stemp),
          formula_surv = ~ 1,
          n_ages = 6,
          .toggles = list(est_surv = 1,
                          time_ar = 1),
          init = "pathfinder")

loos <- c(loos, list("drm_5" = drm_5$stanfit$loo()))

loo::loo_compare(loos)

## instantaneous fishing mortality rates
fmat <-
  system.file("fmat.rds", package = "drmr") |>
  readRDS()

f_train <- fmat[, years_train]
f_test  <- fmat[, years_test]

drm_6 <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "patch",
          family = "gamma",
          seed = 202505,
          formula_zero = ~ 1 + c_hauls + c_btemp + c_stemp,
          formula_rec = ~ 1 + c_stemp + I(c_stemp * c_stemp),
          n_ages = NROW(f_train),
          f_mort = f_train,
          m = 0.25,
          .toggles = list(time_ar = 1),
          init = "pathfinder")

loos <- c(loos, list("drm_6" = drm_6$stanfit$loo()))

drm_7 <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "patch",
          family = "gamma",
          seed = 202505,
          formula_zero = ~ 1 + c_hauls + c_btemp + c_stemp,
          formula_rec = ~ 1,
          formula_surv = ~ 1 + c_btemp + I(c_btemp * c_btemp),
          .toggles = list(time_ar = 1,
                          est_surv = 1),
          init = "pathfinder")

loos <- c(loos, list("drm_7" = drm_7$stanfit$loo()))

loo::loo_compare(loos)

##--- * MCMC diagnostics ----

viz_pars <-
  c("tau[1]",
    "alpha[1]",
    "phi[1]")

mcmc_trace(drm_3$stanfit$draws(variables = viz_pars))
mcmc_dens_overlay(drm_3$stanfit$draws(variables = viz_pars))

mcmc_trace(drm_3$stanfit$draws(variables = c("beta_r")))
mcmc_dens_overlay(drm_3$stanfit$draws(variables = c("beta_r")))

mcmc_trace(drm_3$stanfit$draws(variables = c("beta_t")))
mcmc_dens_overlay(drm_3$stanfit$draws(variables = c("beta_t")))

##--- Viz relationships ----

## * make this into a function!
## ** that is far from easy!

## recruitment

newdata_rec <- data.frame(c_stemp =
                            seq(from = quantile(dat_train$c_stemp, .05),
                                to = quantile(dat_train$c_stemp, .95),
                                length.out = 200))

rec_samples <- marg_rec(drm_3, newdata_rec)

rec_samples <- rec_samples |>
  mutate(stemp = c_stemp + avgs["stemp"])

rec_summary <-
  rec_samples |>
  group_by(stemp) |>
  summarise(ll = quantile(recruitment, probs = .05),
            l = quantile(recruitment, probs = .1),
            m = median(recruitment),
            u = quantile(recruitment, probs = .9),
            uu = quantile(recruitment, probs = .95)) |>
  ungroup()

set.seed(1111)
rec_samples |>
  group_by(stemp) |>
  sample_n(size = 200) |>
  ungroup() |>
  ggplot(data = _,
       aes(x = stemp,
           y = recruitment,
           group = iter)) +
  geom_ribbon(data = rec_summary,
              aes(x = stemp,
                  ymin = l, ymax = u),
              inherit.aes = FALSE,
              fill = 2,
              alpha = .2,
              linewidth = 1.2) +
  geom_ribbon(data = rec_summary,
              aes(x = stemp,
                  ymin = ll, ymax = uu),
              inherit.aes = FALSE,
              fill = 2,
              alpha = .2,
              linewidth = 1.2) +
  geom_line(alpha = .05) +
  geom_line(data = rec_summary,
            aes(x = stemp, y = m),
            inherit.aes = FALSE,
            color = 2,
            linewidth = 1.2) +
  theme_bw()

##--- * survival ----

newdata_surv <- data.frame(c_btemp =
                             seq(from = quantile(dat_train$c_btemp, .05),
                                 to = quantile(dat_train$c_btemp, .95),
                                 length.out = 200))

surv_samples <- marg_surv(drm_7, newdata_surv)

surv_samples <- surv_samples |>
  mutate(btemp = c_btemp + avgs["btemp"])

surv_summary <-
  surv_samples |>
  group_by(btemp) |>
  summarise(ll = quantile(survival, probs = .05),
            l = quantile(survival, probs = .1),
            m = median(survival),
            u = quantile(survival, probs = .9),
            uu = quantile(survival, probs = .95)) |>
  ungroup()

set.seed(1111)
surv_samples |>
  group_by(btemp) |>
  sample_n(size = 200) |>
  ungroup() |>
  ggplot(data = _,
       aes(x = btemp,
           y = survival,
           group = iter)) +
  geom_ribbon(data = surv_summary,
              aes(x = btemp,
                  ymin = l, ymax = u),
              inherit.aes = FALSE,
              fill = 2,
              alpha = .2,
              linewidth = 1.2) +
  geom_ribbon(data = surv_summary,
              aes(x = btemp,
                  ymin = ll, ymax = uu),
              inherit.aes = FALSE,
              fill = 2,
              alpha = .2,
              linewidth = 1.2) +
  geom_line(alpha = .05) +
  geom_line(data = surv_summary,
            aes(x = btemp, y = m),
            inherit.aes = FALSE,
            color = 2,
            linewidth = 1.2) +
  theme_bw()

##--- * prob of absence ----

newdata_abs <- expand.grid(c_hauls = 0,
                           c_btemp =
                             seq(from = quantile(dat_train$c_btemp, .05),
                                 to = quantile(dat_train$c_btemp, .95),
                                 length.out = 20),
                           c_stemp =
                             seq(from = quantile(dat_train$c_stemp, .05),
                                 to = quantile(dat_train$c_stemp, .95),
                                 length.out = 20))

abs_samples <- marg_pabs(drm_3, newdata_abs)

abs_samples <- abs_samples |>
  mutate(btemp = c_btemp + avgs["btemp"]) |>
  mutate(stemp = c_stemp + avgs["stemp"])

abs_summary <-
  abs_samples |>
  group_by(btemp, stemp, c_hauls) |>
  summarise(ll = quantile(prob_abs, probs = .05),
            l = quantile(prob_abs, probs = .1),
            m = median(prob_abs),
            u = quantile(prob_abs, probs = .9),
            uu = quantile(prob_abs, probs = .95)) |>
  ungroup()

ggplot(data = abs_summary,
       aes(x = btemp,
           y = stemp,
           fill = m)) +
  geom_raster() +
  scale_fill_viridis_c(option = "H") +
  theme_bw()

##--- forecasting ----

##--- * DRM ----

forecast_0 <- predict_drm(drm = drm_0,
                          ntime_for =
                            length(unique(dat_test$year)),
                          new_data = dat_test,
                          seed = 125,
                          cores = 4)

forecast_1 <- predict_drm(drm = drm_1,
                          ntime_for =
                            length(unique(dat_test$year)),
                          new_data = dat_test,
                          seed = 125,
                          cores = 4)

forecast_2 <- predict_drm(drm = drm_2,
                          ntime_for =
                            length(unique(dat_test$year)),
                          new_data = dat_test,
                          past_data = filter(dat_train,
                                             year == max(year)),
                          seed = 125,
                          cores = 4)

forecast_3 <- predict_drm(drm = drm_3,
                          ntime_for =
                            length(unique(dat_test$year)),
                          new_data = dat_test,
                          past_data = filter(dat_train,
                                             year == max(year)),
                          seed = 125,
                          cores = 4)

forecast_4 <- predict_drm(drm = drm_4,
                          ntime_for =
                            length(unique(dat_test$year)),
                          new_data = dat_test,
                          past_data = filter(dat_train,
                                             year == max(year)),
                          seed = 125,
                          cores = 4)

forecast_5 <- predict_drm(drm = drm_5,
                          ntime_for =
                            length(unique(dat_test$year)),
                          new_data = dat_test,
                          past_data = filter(dat_train,
                                             year == max(year)),
                          seed = 125,
                          cores = 4)

forecast_6 <- predict_drm(drm = drm_6,
                          ntime_for =
                            length(unique(dat_test$year)),
                          new_data = dat_test,
                          past_data = filter(dat_train,
                                             year == max(year)),
                          f_test = f_test,
                          seed = 125,
                          cores = 4)

forecast_7 <- predict_drm(drm = drm_7,
                          ntime_for =
                            length(unique(dat_test$year)),
                          new_data = dat_test,
                          past_data = filter(dat_train,
                                             year == max(year)),
                          seed = 125,
                          cores = 4)

##--- * obtaining the summary for predictions ----

all_forecasts <-
  ls(pattern = "^forecast_")
all_drms <-
  ls(pattern = "^drm_")

forecasts_summary <-
  Map(f = \(x, nm) {
    fct <- get(x)
    fct$draws(variables = "y_proj",
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
      arrange(pair) |>
      mutate(model = nm,
             .before = 1)
  }, x = all_forecasts, nm = all_drms)


forecasts_summary <-
  bind_rows(forecasts_summary)

forecasts_summary <-
  dat_test |>
  select(dens, lat_floor, patch, year) |>
  mutate(pair = row_number()) |>
  left_join(forecasts_summary, by = "pair")

forecasts_summary |>
  mutate(bias = dens - m) |>
  mutate(rmse = bias * bias) |>
  mutate(abias = abs(bias)) |>
  mutate(is1 = int_score(dens, l = ll, u = uu, alpha = .1)) |>
  mutate(is2 = int_score(dens, l = l, u = u, alpha = .2)) |>
  mutate(cvg1 = data.table::between(dens, ll, uu)) |>
  mutate(cvg2 = data.table::between(dens, l, u)) |>
  ungroup() |>
  group_by(model) |>
  summarise(across(bias:cvg2, mean)) |>
  ungroup() |>
  rename_all(toupper) |>
  rename("Model" = "MODEL",
         "IS (90%)" = "IS1",
         "IS (80%)" = "IS2",
         "PIC (90%)" = "CVG1",
         "PIC (80%)" = "CVG2")
