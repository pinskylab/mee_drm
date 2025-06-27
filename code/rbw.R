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

map <- my_dt |>
  filter(year == max(year)) |>
  select(id) |>
  distinct()

polygons <- map |>
  st_geometry()

polygons |>
  st_area() |>
  units::set_units("km^2") |>
  summary()

my_dt <- st_drop_geometry(my_dt)

##--- splitting data for validation ----

## reserving 5 years for forecast assessment
first_year_forecast <- max(my_dt$year) - 7

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
  mutate(dens = 100 * y,
         .before = y)

dat_test <- dat_test |>
  mutate(dens = 100 * y,
         .before = y)

chains <- 4
cores <- 4

##--- fig 1 ----

world <-
  rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

filter <- st_bbox(map) |>
  st_as_sfc()

bk <- st_intersection(world, filter) |>
  st_as_sfc()

map |>
  left_join(filter(dat_train, year == 2011),
            by = "id") |>
  ggplot(data = _,
         aes(fill = dens)) +
  geom_sf(data = bk, inherit.aes = FALSE) +
  geom_sf(alpha = .8) +
  scale_fill_viridis_c(option = "H") +
  theme_bw()

##--- fitting DRMs ----

adj_mat <- gen_adj(st_geometry(polygons))

## row-standardized matrix
adj_mat <-
  t(apply(adj_mat, 1, \(x) x / (sum(x))))

##--- * Recruitment ----

n_ages <- 12

init_ldens <-
  dat_train |>
  filter(year < 1985) |>
  summarise(dens = sum(area * dens) / sum(area)) |>
  pull(dens)

init_ldens <-
  init_ldens *
  seq(.9, .1, len = n_ages) / sum(seq(.9, .1, len = n_ages))

init_ldens <- log(init_ldens[-1])

drm_rec <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "id",
          family = "gamma",
          seed = 202505,
          formula_zero = ~ 1 + n_routes,
          formula_rec = ~ 1 + c_tavg + I(c_tavg * c_tavg),
          formula_surv = ~ 1,
          n_ages = 12,
          adj_mat = adj_mat, ## A matrix for movement routine
          init_data = init_ldens,
          ages_movement = c(0, 0,
                            rep(1, n_ages - 4),
                            0, 0), ## ages allowed to move
          .toggles = list(ar_re = "rec",
                          movement = 1,
                          est_surv = 1,
                          est_init = 0,
                          minit = 0),
          .priors = list(pr_phi_a = 1, pr_phi_b = .1,
                         pr_alpha_a = 4.2, pr_alpha_b = 5.8,
                         pr_zeta_a = 7, pr_zeta_b = 3))

##--- comparing some priors and posteriors ----

drm_rec$stanfit$draws(variables = c("phi")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dgamma(x, shape = 1, rate = .1),
                xlim = c(0, 3),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_rec$stanfit$draws(variables = c("alpha")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dbeta(x, shape1 = 4.2, shape2 = 5.8),
                xlim = c(0, 1),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_rec$stanfit$draws(variables = c("zeta")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dbeta(x, shape1 = 7, shape2 = 3),
                xlim = c(0, 1),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_rec$stanfit$draws(variables = c("xi")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) {
    y <- - x
    dnorm(log(y),
          mean = drm_rec$data$pr_lmxi_mu,
          sd = drm_rec$data$pr_lmxi_sd) / y
  },
  xlim = c(-4, -1e-16),
  n = 501,
  inherit.aes = FALSE,
  color = 2,
  lwd = 1.2)

##--- Viz relationships ----

## * make this into a function!
## ** that is far from easy!

## recruitment

newdata_rec <- data.frame(c_tavg =
                            seq(from = quantile(dat_train$c_tavg, .05),
                                to = quantile(dat_train$c_tavg, .95),
                                length.out = 200),
                          c_lon = 0,
                          c_lat = 0)

rec_samples <- marg_rec(drm_rec, newdata_rec)

rec_samples <- rec_samples |>
  mutate(stemp = c_tavg + avgs["tavg"])

rec_summary <-
  rec_samples |>
  group_by(stemp) |>
  summarise(ll = quantile(recruitment, probs = .05),
            l = quantile(recruitment, probs = .1),
            m = median(recruitment),
            u = quantile(recruitment, probs = .9),
            uu = quantile(recruitment, probs = .95)) |>
  ungroup() |>
  mutate(model = "drm_rec")

rec_fig <-
  ggplot(data = rec_summary,
         aes(x = stemp,
             y = m)) +
  geom_ribbon(aes(ymin = l, ymax = u),
              fill = "gray50",
              color = "transparent",
              linewidth = 1.2) +
  geom_line(color = "white", linewidth = 1.2) +
  theme_bw() +
  guides(fill = "none") +
  labs(color = "Model",
       fill = "Model",
       x = "Air temperature (in Celsius)",
       y = "Est. recruitment (per km2)") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.175, 0.75))

rec_fig

gratio <- 0.5 * (1 + sqrt(5))

ggsave(filename = "overleaf/img/recruitment-rbw.pdf",
       plot = rec_fig,
       width = 6,
       height = 6 / gratio)

##--- * temperature that optimizes recruitment ----

betas_rec <-
  drm_rec$stanfit$draws(variables = "beta_r",
                         format = "matrix")

max_quad_x(betas_rec[, 2], betas_rec[, 3],
           offset = avgs["tavg"]) |>
  apply(2, quantile, probs = c(.1, .5, .9))

##--- Convergence & estimates ----

## all rhat's look good (no larger than 1.01)
drm_rec$stanfit$summary(variables = c("beta_r", "beta_t",
                                      "alpha", "sigma_t",
                                      "zeta", "phi",
                                      "xi"))

## the different chains are in agreement and converging.
drm_rec$stanfit$draws(variables = c("beta_r", "beta_t")) |>
  mcmc_trace()

drm_rec$stanfit$draws(variables = c("beta_r", "beta_t")) |>
  mcmc_dens_overlay()

drm_rec$stanfit$draws(variables = c("alpha", "sigma_t",
                                    "zeta", "phi")) |>
  mcmc_trace(facet_args = list(labeller = ggplot2::label_parsed))

drm_rec$stanfit$draws(variables = c("alpha", "sigma_t",
                                    "zeta", "phi")) |>
  mcmc_dens_overlay(facet_args = list(labeller = ggplot2::label_parsed))

##--- * Survival ----

drm_surv <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "id",
          family = "gamma",
          seed = 202505,
          formula_zero = ~ 1 + n_routes,
          formula_rec = ~ 1,
          formula_surv = ~ 1 + c_tavg + I(c_tavg * c_tavg),
          n_ages = n_ages,
          adj_mat = adj_mat, ## A matrix for movement routine
          init_data = init_ldens,
          ages_movement = c(0, 0,
                            rep(1, n_ages - 4),
                            0, 0), ## ages allowed to move
          .toggles = list(ar_re = "rec",
                          est_surv = 1,
                          movement = 1,
                          est_init = 0,
                          minit = 0),
          .priors = list(pr_phi_a = 1, pr_phi_b = .1,
                         pr_alpha_a = 4.2, pr_alpha_b = 5.8,
                         pr_zeta_a = 7, pr_zeta_b = 3))

##--- Convergence & estimates ----

## all rhat's look good (no larger than 1.01)
drm_surv$stanfit$summary(variables = c("beta_s", "beta_t",
                                       "alpha", "sigma_t",
                                       "zeta", "phi",
                                       "xi"))

## the different chains are in agreement and converging.
drm_surv$stanfit$draws(variables = c("beta_s", "beta_t")) |>
  mcmc_trace()

drm_surv$stanfit$draws(variables = c("beta_s", "beta_t")) |>
  mcmc_dens_overlay()

drm_surv$stanfit$draws(variables = c("alpha", "sigma_t",
                                    "zeta", "phi")) |>
  mcmc_trace(facet_args = list(labeller = ggplot2::label_parsed))

drm_surv$stanfit$draws(variables = c("alpha", "sigma_t",
                                    "zeta", "phi")) |>
  mcmc_dens_overlay(facet_args = list(labeller = ggplot2::label_parsed))

##--- comparing some priors and posteriors ----

drm_surv$stanfit$draws(variables = c("phi")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dgamma(x, shape = 1, rate = .1),
                xlim = c(0, 3),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_surv$stanfit$draws(variables = c("alpha")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dbeta(x, shape1 = 4.2, shape2 = 5.8),
                xlim = c(0, 1),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_surv$stanfit$draws(variables = c("zeta")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dbeta(x, shape1 = 7, shape2 = 3),
                xlim = c(0, 1),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_surv$stanfit$draws(variables = c("beta_r")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dnorm(x),
                xlim = c(-2, 2),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_surv$stanfit$draws(variables = c("xi")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) {
    y <- - x
    dnorm(log(y),
          mean = drm_surv$data$pr_lmxi_mu,
          sd = drm_surv$data$pr_lmxi_sd) / y
  },
  xlim = c(-4, -1e-16),
  n = 501,
  inherit.aes = FALSE,
  color = 2,
  lwd = 1.2)

##--- surv and environment ----

newdata_surv <- data.frame(c_tavg =
                             seq(from = quantile(dat_train$c_tavg, .05),
                                 to = quantile(dat_train$c_tavg, .95),
                                 length.out = 200),
                           c_lon = 0,
                           c_lat = 0)

surv_samples <- marg_surv(drm_surv, newdata_surv)

surv_samples <- surv_samples |>
  mutate(btemp = c_tavg + avgs["tavg"])

surv_summary <-
  surv_samples |>
  group_by(btemp) |>
  summarise(ll = quantile(survival, probs = .05),
            l = quantile(survival, probs = .1),
            m = median(survival),
            u = quantile(survival, probs = .9),
            uu = quantile(survival, probs = .95)) |>
  ungroup()

surv_fig <-
  ggplot(data = surv_summary,
         aes(x = btemp,
             y = m)) +
  geom_ribbon(aes(ymin = l, ymax = u),
              fill = "gray50",
              color = "transparent",
              linewidth = 1.2) +
  geom_line(color = "white", linewidth = 1.2) +
  theme_bw() +
  labs(x = "Air temperature (in Celsius)",
       y = "Est. survival")

surv_fig

ggsave(filename = "overleaf/img/surv-rbw.pdf",
       plot = surv_fig,
       width = 6,
       height = 6 / gratio)

##--- Panel for recruitment and survival ----

rec_fig + surv_fig +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10))

ggsave(filename = "overleaf/img/rec_surv-rbw.pdf",
       width = 7,
       height = .75 * 7 / gratio)

##--- ** obtaining the SBT that maximizes surival ----

betas_surv <-
  drm_surv$stanfit$draws(variables = "beta_s",
                         format = "matrix")

max_quad_x(betas_surv[, 2], betas_surv[, 3],
           offset = avgs["tavg"]) |>
  apply(2, quantile, probs = c(.1, .5, .9))

##--- * SDM ----

sdm <-
  fit_sdm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "id",
          family = "lognormal",
          seed = 202505,
          formula_zero = ~ 1 + n_routes,
          formula_dens = ~ 1 + c_tavg + I(c_tavg * c_tavg),
          .toggles = list(rho_mu = 0))

##--- forecasting ----

##--- * DRM ----

forecast_rec <- predict_drm(drm = drm_rec,
                            new_data = dat_test,
                            past_data = filter(dat_train,
                                               year == max(year)),
                            seed = 125,
                            cores = 4)

forecast_surv <- predict_drm(drm = drm_surv,
                             new_data = dat_test,
                             past_data = filter(dat_train,
                                                year == max(year)),
                             seed = 125,
                             cores = 4)

forecast_sdm <-
  predict_sdm(sdm = sdm,
              new_data = dat_test,
              seed = 125,
              cores = 4)

##--- Viz predicted and observed ----

fitted_rec <-
  drm_rec$stanfit$draws(variables = "y_pp",
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

forecast_rec <-
  forecast_rec$draws(variables = "y_proj",
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

forecast_rec <-
  dat_test |>
  select(year, id, lat, lon, dens) |>
  mutate(pair = row_number(), .before = 1) |>
  left_join(forecast_rec, by = "pair") |>
  select(- pair)

fitted_rec <-
  dat_train |>
  select(year, id, lat, lon, dens) |>
  bind_cols(select(fitted_rec, -pair))

fitted_surv <-
  drm_surv$stanfit$draws(variables = "y_pp",
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

forecasts_surv <-
  forecast_surv$draws(variables = "y_proj",
                      format = "draws_df") |>
  tidyr::pivot_longer(cols = starts_with("y_proj"),
                      names_to = "pair",
                      values_to = "expected") |>
  group_by(pair) |>
  summarise(ll = quantile(expected, probs = .05, na.rm = TRUE),
            l = quantile(expected, probs = .1, na.rm = TRUE),
            m = median(expected, na.rm = TRUE),
            u = quantile(expected, probs = .9, na.rm = TRUE),
            uu = quantile(expected, probs = .95, na.rm = TRUE)) |>
  ungroup() |>
  mutate(pair = gsub("\\D", "", pair)) |>
  mutate(pair = as.integer(pair)) |>
  arrange(pair)

forecast_surv <-
  dat_test |>
  select(year, id, lat, lon, dens) |>
  mutate(pair = row_number(), .before = 1) |>
  left_join(forecasts_surv, by = "pair") |>
  select(- pair)

fitted_surv <-
  dat_train |>
  select(year, id, lat, lon, dens) |>
  bind_cols(select(fitted_surv, -pair))

fitted_sdm <-
  sdm$stanfit$draws(variables = "y_pp",
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

for_sdm <-
  dat_test |>
  select(year, id, lat, lon, dens) |>
  mutate(pair = row_number(), .before = 1) |>
  left_join(for_sdm, by = "pair") |>
  select(- pair)

fitted_sdm <-
  dat_train |>
  select(year, id, lat, lon, dens) |>
  bind_cols(select(fitted_sdm, -pair))

##--- table ----

for_sdm |>
  mutate(model = "SDM") |>
  bind_rows(
      forecast_rec |>
      mutate(model = "DRM (rec)"),
      forecast_surv |>
      mutate(model = "DRM (surv)")
  ) |>
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
  arrange(RMSE) |>
  print() |>
  xtable::xtable(caption = "Forecasting skill according to different metrics",
                 digits = 2) |>
  print(include.rownames = FALSE)

##--- Figure ? ----

set.seed(2025)
ids <- sample(seq_len(max(dat_train$id)),
              size = 5)

bind_rows(fitted_sdm, for_sdm) |>
  mutate(model = "SDM") |>
  bind_rows(
      bind_rows(fitted_rec, forecast_rec) |>
      mutate(model = "DRM (rec)"),
      bind_rows(fitted_surv, forecast_surv) |>
      mutate(model = "DRM (surv)")
  ) |>
  ## filter(model != "DRM (surv)") |>
  filter(id %in% ids) |>
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
  facet_grid(id ~ model) +
  guides(color = "none",
         fill = "none") +
  labs(x = "Year",
       y = "Density") +
  scale_y_continuous(breaks = scales::trans_breaks(identity, identity,
                                                   n = 3),
                     trans = "log1p") +
  theme_bw()

ggsave(filename = "overleaf/img/forecasts_birds.pdf",
       width = 6,
       height = 8)

##--- saving some stuff ----

saveRDS(sdm, file = "data/birds/models/sdm.rds")
saveRDS(drm_rec, file = "data/birds/models/drm_rec.rds")
saveRDS(drm_surv, file = "data/birds/models/drm_surv.rds")

wood_out <-
  bind_rows(fitted_sdm, for_sdm) |>
  mutate(model = "SDM") |>
  bind_rows(
      bind_rows(fitted_rec, forecast_rec) |>
      mutate(model = "DRM (rec)"),
      bind_rows(fitted_surv, forecast_surv) |>
      mutate(model = "DRM (surv)")      
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

write_dataset(wood_out,
              path = "data/birds/models/predfor.parquet")

wood_out |>
  filter(year %in% c(1985, 1995, 2005, 2015)) |>
  tidyr::pivot_longer(cols = c("dens", "m")) |>
  mutate(model = if_else(name == "dens", "Observed", model)) |>
  mutate(model = factor(model,
                        levels = c("Observed", "DRM (rec)",
                                   "DRM (surv)", "SDM"))) |>
  ggplot(data = _) +
  geom_sf(aes(fill = value)) +
  scale_fill_viridis_c(option = "H", trans = "log1p") +
  facet_grid(model ~ year) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

wood_out |>
  tidyr::pivot_longer(cols = c("dens", "m")) |>
  mutate(model = if_else(name == "dens", "Observed", model)) |>
  mutate(model = factor(model,
                        levels = c("Observed", "DRM (rec)",
                                   "DRM (surv)", "SDM"))) |>
  filter(model == "DRM (rec)") |>
  filter(year > 1980) |>
  ggplot(data = _) +
  geom_sf(aes(fill = value)) +
  scale_fill_viridis_c(option = "H", trans = "log1p") +
  facet_wrap(~ year) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

wood_out |>
  filter(model == "DRM (rec)") |>
  mutate(bias = m - dens) |>
  mutate(mape = 100 * abs(bias) / dens) |>
  filter(!is.nan(mape)) |>
  filter(year > 1980) |>
  ggplot(data = _) +
  geom_sf(aes(fill = mape)) +
  scale_fill_viridis_c(option = "H", trans = "log") +
  facet_wrap(~ year) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

wood_out |>
  filter(model == "DRM (rec)") |>
  mutate(bias = m - dens) |>
  filter(bias > -8) |>
  filter(year > 1980) |>
  ggplot(data = _) +
  geom_sf(aes(fill = bias)) +
  scale_fill_viridis_c(option = "H") +
  facet_wrap(~ year) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
