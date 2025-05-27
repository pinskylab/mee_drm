library(dplyr)
library(ggplot2)
library(bayesplot)
library(sf)
library(drmr)
library(patchwork)
library(cmdstanr)

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

adj_mat <- gen_adj(st_buffer(st_geometry(polygons),
                             dist = 2500))

## row-standardized matrix
adj_mat <-
  t(apply(adj_mat, 1, \(x) x / (sum(x))))

## instantaneous fishing mortality rates
fmat <-
  system.file("fmat.rds", package = "drmr") |>
  readRDS()

f_train <- fmat[, years_train]
f_test  <- fmat[, years_test]

drm_rec <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "patch",
          family = "gamma",
          seed = 202505,
          formula_zero = ~ 1 + c_hauls + c_btemp + c_stemp,
          formula_rec = ~ 1 + c_stemp + I(c_stemp * c_stemp),
          formula_surv = ~ 1,
          f_mort = f_train,
          n_ages = NROW(f_train),
          adj_mat = adj_mat, ## A matrix for movement routine
          ages_movement = c(0, 0,
                            rep(1, 12),
                            0, 0), ## ages allowed to move
          .toggles = list(time_ar = 1,
                          est_surv = 1,
                          movement = 1),
          init = "pathfinder")

## estimates

drm_rec$stanfit$summary(variables = c("phi", "beta_r"))

##--- Viz relationships ----

## * make this into a function!
## ** that is far from easy!

## recruitment

newdata_rec <- data.frame(c_stemp =
                            seq(from = quantile(dat_train$c_stemp, .05),
                                to = quantile(dat_train$c_stemp, .95),
                                length.out = 200))

rec_samples_3 <- marg_rec(drm_rec, newdata_rec)

rec_samples_3 <- rec_samples_3 |>
  mutate(stemp = c_stemp + avgs["stemp"])

rec_summary <-
  rec_samples_3 |>
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
       x = "SST (in Celsius)",
       y = "Est. recruitment (per km2)") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.175, 0.75))

drm_surv <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "patch",
          family = "gamma",
          seed = 202505,
          formula_zero = ~ 1 + c_hauls + c_btemp + c_stemp,
          formula_rec = ~ 1,
          formula_surv = ~ 1 + c_btemp + I(c_btemp * c_btemp),
          f_mort = f_train,
          n_ages = NROW(f_train),
          adj_mat = adj_mat, ## A matrix for movement routine
          ages_movement = c(0, 0,
                            rep(1, 12),
                            0, 0), ## ages allowed to move
          .toggles = list(time_ar = 1,
                          est_surv = 1,
                          movement = 1),
          init = "pathfinder")

##--- * SDM ----

sdm <-
  fit_sdm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "patch",
          family = "lognormal",
          seed = 202505,
          formula_zero = ~ 1 + c_hauls + c_btemp + c_stemp,
          formula_dens = ~ 1 + c_stemp + I(c_stemp * c_stemp),
          init = "pathfinder")

##--- * model comparison ----

loos <- list("drec"  = drm_rec$stanfit$loo(),
             "dsurv" = drm_surv$stanfit$loo(),
             "sdm"   = sdm$stanfit$loo())

loos_out <- loo::loo_compare(loos)

##--- * some quantities for model comparison ----

aux_qt <-
  data.frame(Model = names(loos),
             delta_LOOIC = loos_out[order(rownames(loos_out)), 1],
             LOOIC = sapply(loos, \(x) x$estimates[3, 1]))

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
  select(year, patch, lat_floor, dens) |>
  mutate(pair = row_number(), .before = 1) |>
  left_join(forecast_rec, by = "pair") |>
  select(- pair)

fitted_rec <-
  dat_train |>
  select(year, patch, lat_floor, dens) |>
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
  select(year, patch, lat_floor, dens) |>
  mutate(pair = row_number(), .before = 1) |>
  left_join(forecasts_surv, by = "pair") |>
  select(- pair)

fitted_surv <-
  dat_train |>
  select(year, patch, lat_floor, dens) |>
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
  select(year, patch, lat_floor, dens) |>
  mutate(pair = row_number(), .before = 1) |>
  left_join(for_sdm, by = "pair") |>
  select(- pair)

fitted_sdm <-
  dat_train |>
  select(year, patch, lat_floor, dens) |>
  bind_cols(select(fitted_sdm, -pair))

##--- Figure 1 ----

bind_rows(fitted_sdm, for_sdm) |>
  mutate(model = "SDM") |>
  bind_rows(
      bind_rows(fitted_rec, forecast_rec) |>
      mutate(model = "DRM (rec)"),
      bind_rows(fitted_surv, forecast_surv) |>
      mutate(model = "DRM (surv)")      
  ) |>
  filter(year > 1985) |>
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
  facet_grid(patch ~ model) +
  scale_y_continuous(breaks = scales::trans_breaks(identity, identity,
                                                   n = 3),
                     trans = "log1p") +
  theme_bw()

aux_qt <-
  mutate(aux_qt,
         Model = case_when(Model == "sdm"  ~ "SDM",
                           Model == "drec" ~ "DRM (rec)",
                           TRUE            ~ "DRM (surv)"))

sdm_forecast |>
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
  left_join(aux_qt,
            by = "Model") |>
  arrange(LOOIC) |>
  relocate(LOOIC, delta_LOOIC, .after = "Model") |>
  print() |>
  xtable::xtable(caption = "Forecasting skill according to different metrics",
                 digits = 2) |>
  print(include.rownames = FALSE)


##--- * MCMC diagnostics ----

viz_pars <-
  c("tau[1]",
    "alpha[1]",
    "phi[1]")

mcmc_trace(drm_rec$stanfit$draws(variables = viz_pars))
mcmc_dens_overlay(drm_rec$stanfit$draws(variables = viz_pars))

mcmc_trace(drm_rec$stanfit$draws(variables = c("beta_r")))
mcmc_dens_overlay(drm_rec$stanfit$draws(variables = c("beta_r")))

mcmc_trace(drm_rec$stanfit$draws(variables = c("beta_t")))
mcmc_dens_overlay(drm_rec$stanfit$draws(variables = c("beta_t")))

mcmc_trace(drm_surv$stanfit$draws(variables = viz_pars))
mcmc_dens_overlay(drm_surv$stanfit$draws(variables = viz_pars))

mcmc_trace(drm_surv$stanfit$draws(variables = c("beta_s")))
mcmc_dens_overlay(drm_surv$stanfit$draws(variables = c("beta_s")))

mcmc_trace(drm_surv$stanfit$draws(variables = c("beta_t")))
mcmc_dens_overlay(drm_surv$stanfit$draws(variables = c("beta_t")))

##--- Viz relationships ----

## * make this into a function!
## ** that is far from easy!

## recruitment

newdata_rec <- data.frame(c_stemp =
                            seq(from = quantile(dat_train$c_stemp, .05),
                                to = quantile(dat_train$c_stemp, .95),
                                length.out = 200))

rec_samples_3 <- marg_rec(drm_rec, newdata_rec)

rec_samples_3 <- rec_samples_3 |>
  mutate(stemp = c_stemp + avgs["stemp"])

rec_summary <-
  rec_samples_3 |>
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
       x = "SST (in Celsius)",
       y = "Est. recruitment (per km2)") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.175, 0.75))

gratio <- 0.5 * (1 + sqrt(5))

ggsave(filename = "overleaf/img/recruitment.pdf",
       plot = rec_fig,
       width = 6,
       height = 6 / gratio)

##--- * survival ----

newdata_surv <- data.frame(c_btemp =
                             seq(from = quantile(dat_train$c_btemp, .05),
                                 to = quantile(dat_train$c_btemp, .95),
                                 length.out = 200))

surv_samples <- marg_surv(drm_surv, newdata_surv)

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
  labs(x = "SBT (in Celsius)",
       y = "Est. survival")

ggsave(filename = "overleaf/img/surv.pdf",
       plot = surv_fig,
       width = 6,
       height = 6 / gratio)

##--- ** obtaining the SBT that maximizes surival ----

betas_s_7 <-
  drm_surv$stanfit$draws(variables = "beta_s",
                         format = "matrix")

max_quad_x(betas_s_7[, 2], betas_s_7[, 3],
           offset = avgs["stemp"]) |>
  apply(2, quantile, probs = c(.1, .5, .9))

##--- Panel for recruitment and survival ----

rec_fig + surv_fig +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10))

ggsave(filename = "overleaf/img/rec_surv.pdf",
       width = 7,
       height = .75 * 7 / gratio)
