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

sdm <-
  fit_sdm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "id",
          seed = 202505,
          formula_zero = ~ 1 + c_tavg + c_lon + c_lat,
          formula_dens = ~ 1 + c_tavg + I(c_tavg * c_tavg),
          .toggles = list(time_ar = 1),
          .priors = list(pr_phi_a = .1, pr_phi_b = .1))

sdm_ll <-
  fit_sdm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "id",
          family = "loglogistic",
          seed = 202505,
          formula_zero = ~ 1 + c_tavg + c_lon + c_lat,
          formula_dens = ~ 1 + c_tavg + I(c_tavg * c_tavg),
          .toggles = list(time_ar = 1),
          .priors = list(pr_phi_a = .1, pr_phi_b = .1))

sdm_ln <-
  fit_sdm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "id",
          family = "lognormal",
          seed = 202505,
          formula_zero = ~ 1 + c_tavg + c_lon + c_lat,
          formula_dens = ~ 1 + c_tavg + I(c_tavg * c_tavg),
          .toggles = list(time_ar = 1),
          .priors = list(pr_phi_a = .1, pr_phi_b = .1))

sdm$stanfit$summary(variables = c("beta_t", "beta_r",
                                  "tau", "alpha", "phi"))
sdm_ll$stanfit$summary(variables = c("beta_t", "beta_r",
                                     "tau", "alpha", "phi"))
sdm_ln$stanfit$summary(variables = c("beta_t", "beta_r",
                                     "tau", "alpha", "phi"))

loos <- list("sdm (gamma)"   = sdm$stanfit$loo(),
             "sdm (LL)"      = sdm_ll$stanfit$loo())

loo::loo_compare(loos)

##--- fitting DRMs ----

baseline <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "id",
          seed = 202505,
          .priors = list(pr_phi_a = .1, pr_phi_b = .1))

drm_1 <-
  fit_drm(.data = dat_train,
          y_col = "dens",
          time_col = "year",
          site_col = "id",
          seed = 202505,
          formula_zero = ~ 1 + c_tavg + c_lon + c_lat,
          formula_rec = ~ 1 +  c_tavg + I(c_tavg * c_tavg),
          m = 0.25,
          init = "prior",
          .priors = list(pr_phi_a = .1, pr_phi_b = .1))

drm_2 <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "id",
          seed = 202505,
          formula_zero = ~ 1 + c_tavg + c_lon + c_lat,
          formula_rec = ~ 1 + c_tavg + c_tavg * c_tavg,
          formula_surv = ~ 1,
          m = 0,
          .toggles = list(est_surv = 1),
          init = "prior",
          .priors = list(pr_phi_a = .1, pr_phi_b = .1))

drm_3 <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "id",
          seed = 202505,
          formula_zero = ~ 1 + c_tavg + c_lon + c_lat,
          formula_rec = ~ 1 + c_tavg + I(c_tavg * c_tavg),
          formula_surv = ~ 1,
          .toggles = list(est_surv = 1,
                          time_ar = 1))

adj_mat <- gen_adj(st_geometry(polygons))

## row-standardized matrix
adj_mat <-
  t(apply(adj_mat, 1, \(x) x / (sum(x))))

drm_4 <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "id",
          seed = 202505,
          formula_zero = ~ 1 + c_tavg + c_lon + c_lat,
          formula_rec = ~ 1 + c_tavg + I(c_tavg * c_tavg),
          formula_surv = ~ 1,
          n_ages = 12,
          adj_mat = adj_mat, ## A matrix for movement routine
          ages_movement = c(0, 0, rep(1, 8),
                            0, 0), ## ages allowed to move
          .toggles = list(est_surv = 1,
                          time_ar = 1,
                          movement = 1),
          init = "prior")

drm_4$stanfit$summary(variables = c("phi", "alpha", "zeta",
                                    "beta_r", "beta_s",
                                    "beta_t"))

drm_5 <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "id",
          seed = 202505,
          formula_zero = ~ 1 + c_tavg + c_lon + c_lat,
          formula_rec = ~ 1 + c_tavg + I(c_tavg * c_tavg),
          formula_surv = ~ 1,
          n_ages = 12,
          .toggles = list(est_surv = 1,
                          time_ar = 1),
          .priors = list(pr_phi_a = .1, pr_phi_b = .1))

drm_6 <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "id",
          family = "gamma",
          seed = 202505,
          n_ages = 12,
          formula_zero = ~ 1 + c_tavg + c_lon + c_lat,
          formula_rec = ~ 1,
          formula_surv = ~ 1 + c_tavg + I(c_tavg * c_tavg),
          .toggles = list(time_ar = 1,
                          est_surv = 1),
          .priors = list(pr_phi_a = .1, pr_phi_b = .1))

##--- * model comparison ----

loos <- list("baseline" = baseline$stanfit$loo(),
             "drm_1" = drm_1$stanfit$loo(),
             "drm_2" = drm_2$stanfit$loo(),
             "drm_3" = drm_3$stanfit$loo(),
             "drm_4" = drm_4$stanfit$loo(),
             "drm_5" = drm_5$stanfit$loo(),
             "drm_6" = drm_6$stanfit$loo(),
             "sdm"   = sdm$stanfit$loo())

loos_out <- loo::loo_compare(loos)

##--- * some quantities for model comparison ----

times <- ls(pattern = "^(drm|baseline|sdm)") |>
  sapply(\(x) get(x)$stanfit$time()$total)

times <- times[1:8]

aux_qt <-
  data.frame(Model = names(times),
             LOOIC = loos_out[order(rownames(loos_out)), 1],
             time = times)

##--- forecasting ----

##--- * DRM ----

forecast_0 <- predict_drm(drm = baseline,
                          new_data = dat_test,
                          seed = 125,
                          cores = 4)

forecast_1 <- predict_drm(drm = drm_1,
                          new_data = dat_test,
                          seed = 125,
                          cores = 4)

forecast_2 <- predict_drm(drm = drm_2,
                          new_data = dat_test,
                          past_data = filter(dat_train,
                                             year == max(year)),
                          seed = 125,
                          cores = 4)

forecast_3 <- predict_drm(drm = drm_3,
                          new_data = dat_test,
                          past_data = filter(dat_train,
                                             year == max(year)),
                          seed = 125,
                          cores = 4)

forecast_4 <- predict_drm(drm = drm_4,
                          new_data = dat_test,
                          past_data = filter(dat_train,
                                             year == max(year)),
                          seed = 125,
                          cores = 4)

forecast_5 <- predict_drm(drm = drm_5,
                          new_data = dat_test,
                          past_data = filter(dat_train,
                                             year == max(year)),
                          seed = 125,
                          cores = 4)

forecast_6 <- predict_drm(drm = drm_6,
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

##--- * obtaining the summary for predictions ----

all_forecasts <-
  ls(pattern = "^forecast_")
all_drms <-
  ls(pattern = "^(baseline|drm_|sdm)")[1:8]

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
  select(dens, lat, lon, id, year) |>
  mutate(pair = row_number()) |>
  left_join(forecasts_summary, by = "pair")

## Table 5

forecasts_summary |>
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
  arrange(desc(LOOIC)) |>
  relocate(LOOIC, .after = "Model") |>
  print() |>
  xtable::xtable(caption = "Forecasting skill according to different metrics",
                 digits = 2) |>
  print(include.rownames = FALSE)

## Figure 1 with selected models' forecast

set.seed(202505)
sample_ids <- sample(seq_len(max(forecasts_summary$id)),
                     size = 10) |>
  sort()

filter(forecasts_summary,
       model %in% c("drm_1",
                    "drm_3",
                    "drm_4",
                    "sdm")) |>
  filter(id %in% sample_ids) |>
  ggplot(data = _,
       aes(x = year,
           y = m)) +
  geom_ribbon(aes(ymin = l, ymax = u),
              fill = 2,
              alpha = .4) +
  geom_line() +
  geom_point(aes(y = dens),
             color = 4) +
  facet_grid(rev(id) ~ model,
             scales = "free_y") +
  scale_y_continuous(breaks = scales::trans_breaks(identity, identity,
                                                   n = 3)) +
  theme_bw() +
  labs(x = "Year",
       y = "Predicted density") +
  theme(panel.spacing.y = unit(.125, "in"),
        axis.text.x = element_text(angle = 30))

ggsave(filename = "overleaf/img/forecasts_birds.pdf",
       width = 6,
       height = 8)

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

newdata_rec <- data.frame(c_tavg =
                            seq(from = quantile(dat_train$c_tavg, .05),
                                to = quantile(dat_train$c_tavg, .95),
                                length.out = 200))

rec_samples_3 <- marg_rec(drm_3, newdata_rec)

rec_samples_3 <- rec_samples_3 |>
  mutate(tavg = c_tavg + avgs["tavg"])

rec_samples_4 <- marg_rec(drm_4, newdata_rec)

rec_samples_4 <- rec_samples_4 |>
  mutate(tavg = c_tavg + avgs["tavg"])

rec_summary <-
  rec_samples_3 |>
  group_by(tavg) |>
  summarise(ll = quantile(recruitment, probs = .05),
            l = quantile(recruitment, probs = .1),
            m = median(recruitment),
            u = quantile(recruitment, probs = .9),
            uu = quantile(recruitment, probs = .95)) |>
  ungroup() |>
  mutate(model = "drm_3") |>
  bind_rows(
      rec_samples_4 |>
      group_by(tavg) |>
      summarise(ll = quantile(recruitment, probs = .05),
                l = quantile(recruitment, probs = .1),
                m = median(recruitment),
                u = quantile(recruitment, probs = .9),
                uu = quantile(recruitment, probs = .95)) |>
      ungroup() |>
      mutate(model = "drm_4")
  )

rec_fig <-
  ggplot(data = rec_summary,
         aes(x = tavg,
             y = m,
             color = model,
             fill = model)) +
  geom_ribbon(aes(ymin = l, ymax = u),
              alpha = .4,
              color = "transparent",
              linewidth = 1.2) +
  geom_line(linewidth = 1.2) +
  theme_bw() +
  guides(fill = "none") +
  labs(color = "Model",
       fill = "Model",
       x = "Air temperature (in Celsius)",
       y = "Est. recruitment (per km2)") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.175, 0.75))

gratio <- 0.5 * (1 + sqrt(5))

ggsave(filename = "overleaf/img/recruitment_birds.pdf",
       plot = rec_fig,
       width = 6,
       height = 6 / gratio)

##--- ** obtaining the SST that maximizes the recruitment ----

betas_r_3 <-
  drm_3$stanfit$draws(variables = "beta_r",
                      format = "matrix")

max_quad_x(betas_r_3[, 2], betas_r_3[, 3],
           offset = avgs["tavg"]) |>
  apply(2, quantile, probs = c(.1, .5, .9))

betas_r_4 <-
  drm_4$stanfit$draws(variables = "beta_r",
                      format = "matrix")

max_quad_x(betas_r_4[, 2], betas_r_4[, 3],
           offset = avgs["tavg"]) |>
  apply(2, quantile, probs = c(.1, .5, .9))

##--- * survival ----

newdata_surv <- data.frame(c_tavg =
                             seq(from = quantile(dat_train$c_tavg, .05),
                                 to = quantile(dat_train$c_tavg, .95),
                                 length.out = 200))

surv_samples <- marg_surv(drm_6, newdata_surv)

surv_samples <- surv_samples |>
  mutate(tavg = c_tavg + avgs["tavg"])

surv_summary <-
  surv_samples |>
  group_by(tavg) |>
  summarise(ll = quantile(survival, probs = .05),
            l = quantile(survival, probs = .1),
            m = median(survival),
            u = quantile(survival, probs = .9),
            uu = quantile(survival, probs = .95)) |>
  ungroup()

surv_fig <-
  ggplot(data = surv_summary,
         aes(x = tavg,
             y = m)) +
  geom_ribbon(aes(ymin = l, ymax = u),
              alpha = .4,
              color = "transparent",
              fill = 2,
              linewidth = 1.2) +
  geom_line(linewidth = 1.2) +
  theme_bw() +
  labs(x = "Air temperature (in Celsius)",
       y = "Est. survival")

ggsave(filename = "overleaf/img/surv_birds.pdf",
       plot = surv_fig,
       width = 6,
       height = 6 / gratio)

##--- ** obtaining the SBT that maximizes surival ----

betas_s_7 <-
  drm_6$stanfit$draws(variables = "beta_s",
                      format = "matrix")

max_quad_x(betas_s_7[, 2], betas_s_7[, 3],
           offset = avgs["tavg"]) |>
  apply(2, quantile, probs = c(.1, .5, .9))

##--- Panel for recruitment and survival ----

rec_fig + surv_fig +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10))

ggsave(filename = "overleaf/img/rec_surv_birds.pdf",
       width = 7,
       height = .75 * 7 / gratio)

##--- * prob of absence ----

newdata_abs <- expand.grid(c_hauls = 0,
                           c_tavg =
                             seq(from = quantile(dat_train$c_tavg, .05),
                                 to = quantile(dat_train$c_tavg, .95),
                                 length.out = 20),
                           c_tavg =
                             seq(from = quantile(dat_train$c_tavg, .05),
                                 to = quantile(dat_train$c_tavg, .95),
                                 length.out = 20))

abs_samples <- marg_pabs(drm_3, newdata_abs)

abs_samples <- abs_samples |>
  mutate(tavg = c_tavg + avgs["tavg"]) |>
  mutate(tavg = c_tavg + avgs["tavg"])

abs_summary <-
  abs_samples |>
  group_by(tavg, tavg, c_hauls) |>
  summarise(ll = quantile(prob_abs, probs = .05),
            l = quantile(prob_abs, probs = .1),
            m = median(prob_abs),
            u = quantile(prob_abs, probs = .9),
            uu = quantile(prob_abs, probs = .95)) |>
  ungroup()

ggplot(data = abs_summary,
       aes(x = tavg,
           y = tavg,
           fill = m)) +
  geom_raster() +
  scale_fill_viridis_c(option = "H") +
  theme_bw()
