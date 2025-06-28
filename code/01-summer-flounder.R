library(dplyr)
library(ggplot2)
library(bayesplot)
library(sf)
library(drmr)
library(patchwork)
library(cmdstanr)

source("utils/fit_for.R")

bayesplot::color_scheme_set(scheme = "mix-pink-teal")

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
  mutate(dens = 100 * y / area_km2,
         .before = y)

dat_test <- dat_test |>
  mutate(dens = 100 * y / area_km2,
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

## vizualizing different beta priors
mode_zeta <- .4
conc_zeta1 <- 10
conc_zeta2 <- 3.5

alpha1 <- mode_zeta * (conc_zeta1 - 2) + 1
alpha2 <- mode_zeta * (conc_zeta2 - 2) + 1

beta1 <- (1 - mode_zeta) * (conc_zeta1 - 2) + 1
beta2 <- (1 - mode_zeta) * (conc_zeta2 - 2) + 1

ggplot() +
  stat_function(fun =
                  \(x) dbeta(x, shape1 = alpha1, shape2 = beta1),
                color = 2) +
  stat_function(fun =
                  \(x) dbeta(x, shape1 = alpha2, shape2 = beta2),
                color = 4) +
  theme_bw()

mode_zeta <- .75
conc_zeta1 <- 10
conc_zeta2 <- 3.5

alpha1 <- mode_zeta * (conc_zeta1 - 2) + 1
alpha2 <- mode_zeta * (conc_zeta2 - 2) + 1

beta1 <- (1 - mode_zeta) * (conc_zeta1 - 2) + 1
beta2 <- (1 - mode_zeta) * (conc_zeta2 - 2) + 1

ggplot() +
  stat_function(fun =
                  \(x) dbeta(x, shape1 = alpha1, shape2 = beta1),
                color = 2) +
  stat_function(fun =
                  \(x) dbeta(x, shape1 = alpha2, shape2 = beta2),
                color = 4) +
  theme_bw()

##--- DRM recruitment ----

drm_rec <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "patch",
          family = "gamma",
          seed = 202505,
          formula_zero = ~ 1 + c_hauls,
          formula_rec = ~ 1 + c_stemp + I(c_stemp * c_stemp),
          formula_surv = ~ 1,
          f_mort = f_train,
          n_ages = NROW(f_train),
          adj_mat = adj_mat, ## A matrix for movement routine
          ages_movement = c(0, 0,
                            rep(1, 12),
                            0, 0), ## ages allowed to move
          .toggles = list(ar_re = "rec",
                          movement = 1,
                          est_surv = 1,
                          est_init = 0,
                          minit = 1),
          .priors = list(pr_phi_a = 1, pr_phi_b = .1,
                         pr_alpha_a = 4.2, pr_alpha_b = 5.8,
                         pr_zeta_a = 7, pr_zeta_b = 3))

##--- * code 1 ----

drm_rec$stanfit$summary(variables = c("beta_r"))

##--- Convergence check & parameter estimates ----

## all rhat's look good (no larger than 1.01)
drm_rec$stanfit$summary(variables = c("beta_r"))


drm_rec$stanfit$summary(variables = c("beta_r", "beta_t",
                                      "alpha", "sigma_t",
                                      "xi",
                                      "zeta", "phi"))

## the different chains are in agreement and converging.
drm_rec$stanfit$draws(variables = c("beta_r", "beta_t")) |>
  mcmc_trace()

drm_rec$stanfit$draws(variables = c("beta_r", "beta_t")) |>
  mcmc_dens_overlay()

drm_rec$stanfit$draws(variables = c("alpha", "sigma_t",
                                    "xi",
                                    "zeta", "phi")) |>
  mcmc_trace(facet_args = list(labeller = ggplot2::label_parsed))

drm_rec$stanfit$draws(variables = c("alpha", "sigma_t",
                                    "xi",
                                    "zeta", "phi")) |>
  mcmc_dens_overlay(facet_args = list(labeller = ggplot2::label_parsed))

##--- comparing some priors and posteriors ----

drm_rec$stanfit$draws(variables = c("phi")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dgamma(x,
                                  shape = drm_rec$data$pr_phi_a,
                                  rate = drm_rec$data$pr_phi_b),
                xlim = c(0, 3),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_rec$stanfit$draws(variables = c("alpha")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dbeta(x,
                                 shape1 = drm_rec$data$pr_alpha_a,
                                 shape2 = drm_rec$data$pr_alpha_b),
                xlim = c(0, 1),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_rec$stanfit$draws(variables = c("zeta")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dbeta(x,
                                 shape1 = drm_rec$data$pr_zeta_a,
                                 shape2 = drm_rec$data$pr_zeta_b),
                xlim = c(0, 1),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_rec$stanfit$draws(variables = c("sigma_t")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dlnorm(x,
                                  meanlog = drm_rec$data$pr_lsigma_t_mu,
                                  sdlog = drm_rec$data$pr_lsigma_t_sd),
                xlim = c(0, .5),
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
  xlim = c(-5, -1e-16),
  n = 501,
  inherit.aes = FALSE,
  color = 2,
  lwd = 1.2)

##--- DRM Survival ----

drm_surv <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "patch",
          family = "gamma",
          seed = 202505,
          formula_zero = ~ 1 + c_hauls,
          formula_rec = ~ 1,
          formula_surv = ~ 1 + c_btemp + I(c_btemp * c_btemp),
          f_mort = f_train,
          n_ages = NROW(f_train),
          adj_mat = adj_mat, ## A matrix for movement routine
          ages_movement = c(0, 0,
                            rep(1, 12),
                            0, 0), ## ages allowed to move
          .toggles = list(ar_re = "rec",
                          est_surv = 1,
                          movement = 1,
                          est_init = 0,
                          minit = 1),
          .priors = list(pr_phi_a = 1, pr_phi_b = .1,
                         pr_alpha_a = 4.2, pr_alpha_b = 5.8,
                         pr_zeta_a = 7, pr_zeta_b = 3))

##--- Convergence & estimates ----

## r_hat for beta_s indicate convergence issues
drm_surv$stanfit$summary(variables = c("beta_s", "beta_t",
                                       "alpha", "sigma_t",
                                       "xi",
                                       "zeta", "phi"))

## the different chains are in agreement and converging.
drm_surv$stanfit$draws(variables = c("beta_s", "beta_t")) |>
  mcmc_trace()

drm_surv$stanfit$draws(variables = c("beta_s", "beta_t")) |>
  mcmc_dens_overlay()

drm_surv$stanfit$draws(variables = c("alpha", "sigma_t",
                                     "xi",
                                     "zeta", "phi")) |>
  mcmc_trace(facet_args = list(labeller = ggplot2::label_parsed))

drm_surv$stanfit$draws(variables = c("alpha", "sigma_t",
                                     "xi",
                                     "zeta", "phi")) |>
  mcmc_dens_overlay(facet_args = list(labeller = ggplot2::label_parsed))

##--- comparing some priors and posteriors ----

drm_surv$stanfit$draws(variables = c("phi")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dgamma(x,
                                  shape = drm_surv$data$pr_phi_a,
                                  rate = drm_surv$data$pr_phi_b),
                xlim = c(0, 3),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_surv$stanfit$draws(variables = c("alpha")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dbeta(x,
                                 shape1 = drm_surv$data$pr_alpha_a,
                                 shape2 = drm_surv$data$pr_alpha_b),
                xlim = c(0, 1),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_surv$stanfit$draws(variables = c("zeta")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dbeta(x,
                                 shape1 = drm_surv$data$pr_zeta_a,
                                 shape2 = drm_surv$data$pr_zeta_b),
                xlim = c(0, 1),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_surv$stanfit$draws(variables = c("sigma_t")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dlnorm(x,
                                  meanlog = drm_surv$data$pr_lsigma_t_mu,
                                  sdlog = drm_surv$data$pr_lsigma_t_sd),
                xlim = c(0, .5),
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
  xlim = c(-5, -1e-16),
  n = 501,
  inherit.aes = FALSE,
  color = 2,
  lwd = 1.2)

##--- SDM ----

sdm <-
  fit_sdm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "patch",
          family = "lognormal",
          seed = 202505,
          formula_zero = ~ 1 + c_hauls,
          formula_dens = ~ 1 + c_stemp + I(c_stemp * c_stemp),
          .priors = list(pr_phi_a = 1, pr_phi_b = .1,
                         pr_alpha_a = 4.2, pr_alpha_b = 5.8),
          ## the model is not converging with `rho_mu = 1`
          .toggles = list(rho_mu = 0))

##--- forecasting ----

##--- * DRM ----

forecast_rec <- predict_drm(drm = drm_rec,
                            new_data = dat_test,
                            past_data = filter(dat_train,
                                               year == max(year)),
                            seed = 125,
                            f_test = f_test,
                            cores = 4)

forecast_surv <- predict_drm(drm = drm_surv,
                             new_data = dat_test,
                             past_data = filter(dat_train,
                                                year == max(year)),
                             seed = 125,
                             f_test = f_test,
                             cores = 4)

forecast_sdm <-
  predict_sdm(sdm = sdm,
              new_data = dat_test,
              seed = 125,
              cores = 4)

##--- Viz predicted and observed ----

fitted_rec <-
  fitted_drm(drm_rec, ci_mass = .8)

for_rec <-
  forecast_drm(forecast_rec, ci_mass = .8)

for_rec <-
  dat_test |>
  select(year, patch, lat_floor, dens) |>
  mutate(pair = row_number(), .before = 1) |>
  left_join(for_rec, by = "pair") |>
  select(- pair)

fitted_rec <-
  dat_train |>
  select(year, patch, lat_floor, dens) |>
  bind_cols(select(fitted_rec, -pair))

fitted_surv <-
  fitted_drm(drm_surv, ci_mass = .8)

for_surv <-
  forecast_drm(forecast_surv, ci_mass = .8)

for_surv <-
  dat_test |>
  select(year, patch, lat_floor, dens) |>
  mutate(pair = row_number(), .before = 1) |>
  left_join(for_surv, by = "pair") |>
  select(- pair)

fitted_surv <-
  dat_train |>
  select(year, patch, lat_floor, dens) |>
  bind_cols(select(fitted_surv, -pair))

fitted_sdm <-
  fitted_drm(sdm, ci_mass = .8)

for_sdm <-
  forecast_drm(forecast_sdm, ci_mass = .8)

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

##--- Figure 2 ----

bind_rows(fitted_sdm, for_sdm) |>
  mutate(model = "SDM") |>
  bind_rows(
      bind_rows(fitted_rec, for_rec) |>
      mutate(model = "DRM (rec)"),
      bind_rows(fitted_surv, for_surv) |>
      mutate(model = "DRM (surv)")
  ) |>
  ## filter(model != "DRM (surv)") |>
  ggplot(data = _) +
  geom_vline(xintercept = first_year_forecast,
             lty = 2) +
  geom_ribbon(aes(x = year,
                  ymin = l, ymax = u,
                  fill = model,
                  color = model),
              alpha = .4) +
  geom_line(aes(x = year, y = m, color = model)) +
  geom_point(aes(x = year, y = dens), size = .5) +
  facet_grid(patch ~ model, scales = "free_y") +
  scale_y_continuous(breaks = scales::trans_breaks(identity, identity,
                                                   n = 3),
                     trans = "log1p") +
  theme_bw() +
  guides(color = "none",
         fill = "none") +
  labs(y = "Density (in hundreds of individuals per square-km)",
       x = "Year") +
  theme(strip.background = element_rect(fill = "white"))

ggsave(filename = "overleaf/img/forecast_sf.pdf",
       width = 6,
       height = 7)

##--- Viz relationships ----

## * make this into a function!
## ** that is far from easy!

## recruitment

newdata_rec <- data.frame(c_stemp =
                            seq(from = quantile(dat_train$c_stemp, .05),
                                to = quantile(dat_train$c_stemp, .95),
                                length.out = 200))

rec_samples <- marg_rec(drm_rec, newdata_rec)

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

rec_fig

gratio <- 0.5 * (1 + sqrt(5))

ggsave(filename = "overleaf/img/recruitment.pdf",
       plot = rec_fig,
       width = 6,
       height = 6 / gratio)

##--- surv and environment ----

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

surv_fig

ggsave(filename = "overleaf/img/surv.pdf",
       plot = surv_fig,
       width = 6,
       height = 6 / gratio)

##--- Figure 3 ----
##--- Panel for recruitment and survival ----

rec_fig + surv_fig +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10))

ggsave(filename = "overleaf/img/rec_surv.pdf",
       width = 7,
       height = .75 * 7 / gratio)

##--- ** SST that optimizes recruitment ----

betas_rec <-
  drm_rec$stanfit$draws(variables = "beta_r",
                         format = "matrix")

max_quad_x(betas_rec[, 2], betas_rec[, 3],
           offset = avgs["stemp"]) |>
  apply(2, quantile, probs = c(.1, .5, .9))

##--- ** SBT that maximizes surival ----

betas_surv <-
  drm_surv$stanfit$draws(variables = "beta_s",
                         format = "matrix")

max_quad_x(betas_surv[, 2], betas_surv[, 3],
           offset = avgs["btemp"]) |>
  apply(2, quantile, probs = c(.1, .5, .9))

##--- Table 4 ----

for_sdm |>
  mutate(model = "SDM") |>
  bind_rows(
      for_rec |>
      mutate(model = "DRM (rec)"),
      for_surv |>
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

