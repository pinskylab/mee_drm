library(dplyr)
library(ggplot2)
library(bayesplot)
library(sf)
library(drmr)
library(patchwork)
library(cmdstanr)

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

## initializing with data

init_ldens <-
  dat_train |>
  filter(year < 1987) |>
  summarise(dens = sum(area_km2 * dens) / sum(area_km2)) |>
  pull(dens)

init_ldens <-
  init_ldens *
  seq(.9, .1, len = NROW(f_train)) / sum(seq(.9, .1, len = NROW(f_train)))

init_ldens <- log(init_ldens[-1])

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
          init_data = init_ldens,
          adj_mat = adj_mat, ## A matrix for movement routine
          ages_movement = c(0, 0,
                            rep(1, 12),
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

##--- * temperature that optimizes recruitment ----

betas_rec <-
  drm_rec$stanfit$draws(variables = "beta_r",
                         format = "matrix")

max_quad_x(betas_rec[, 2], betas_rec[, 3],
           offset = avgs["stemp"]) |>
  apply(2, quantile, probs = c(.1, .5, .9))

##--- Convergence & estimates ----

## all rhat's look good (no larger than 1.01)
drm_rec$stanfit$summary(variables = c("beta_r", "beta_t",
                                      "alpha", "sigma_t",
                                      ## "sigma_s",
                                      "zeta", "phi"))

## the different chains are in agreement and converging.
drm_rec$stanfit$draws(variables = c("beta_r", "beta_t")) |>
  mcmc_trace()

drm_rec$stanfit$draws(variables = c("beta_r", "beta_t")) |>
  mcmc_dens_overlay()

drm_rec$stanfit$draws(variables = c("alpha", "sigma_t",
                                    ## "sigma_s",
                                    "zeta", "phi")) |>
  mcmc_trace(facet_args = list(labeller = ggplot2::label_parsed))

drm_rec$stanfit$draws(variables = c("alpha", "sigma_t",
                                    ## "sigma_s",
                                    "zeta", "phi")) |>
  mcmc_dens_overlay(facet_args = list(labeller = ggplot2::label_parsed))

##--- formula 2 ----

## estimating initialization

drm_rec2 <-
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
          .toggles = list(ar_re = "rec",
                          movement = 1,
                          est_surv = 1,
                          est_init = 1,
                          minit = 0),
          .priors = list(pr_phi_a = 1, pr_phi_b = .1,
                         pr_alpha_a = 4.2, pr_alpha_b = 5.8,
                         pr_zeta_a = 7, pr_zeta_b = 3))

##--- comparing some priors and posteriors ----

drm_rec2$stanfit$draws(variables = c("phi")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dgamma(x, shape = 1, rate = .1),
                xlim = c(0, 3),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_rec2$stanfit$draws(variables = c("alpha")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dbeta(x,
                                 shape1 = drm_rec2$data$pr_alpha_a,
                                 shape2 = drm_rec2$data$pr_alpha_b),
                xlim = c(0, 1),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_rec2$stanfit$draws(variables = c("zeta")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dbeta(x,
                                 shape1 = drm_rec2$data$pr_zeta_a,
                                 shape2 = drm_rec2$data$pr_zeta_b),
                xlim = c(0, 1),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_rec2$stanfit$draws(variables = c("sigma_t")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dlnorm(x,
                                  meanlog = drm_rec2$data$pr_lsigma_t_mu,
                                  sdlog = drm_rec2$data$pr_lsigma_t_sd),
                xlim = c(0, .5),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

##--- Viz relationships ----

## * make this into a function!
## ** that is far from easy!

## recruitment

newdata_rec <- data.frame(c_stemp =
                            seq(from = quantile(dat_train$c_stemp, .05),
                                to = quantile(dat_train$c_stemp, .95),
                                length.out = 200))

rec_samples <- marg_rec(drm_rec2, newdata_rec)

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
  mutate(model = "drm_rec2")

rec_fig2 <-
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

rec_fig2

##--- * temperature that optimizes recruitment ----

betas_rec <-
  drm_rec2$stanfit$draws(variables = "beta_r",
                         format = "matrix")

max_quad_x(betas_rec[, 2], betas_rec[, 3],
           offset = avgs["stemp"]) |>
  apply(2, quantile, probs = c(.1, .5, .9))

##--- Convergence & estimates ----

## all rhat's look good (no larger than 1.01)
drm_rec2$stanfit$summary(variables = c("beta_r", "beta_t",
                                       "alpha", "sigma_t",
                                       ## "sigma_s",
                                       "zeta", "phi"))

## the different chains are in agreement and converging.
drm_rec2$stanfit$draws(variables = c("beta_r", "beta_t")) |>
  mcmc_trace()

drm_rec2$stanfit$draws(variables = c("beta_r", "beta_t")) |>
  mcmc_dens_overlay()

drm_rec2$stanfit$draws(variables = c("alpha", "sigma_t",
                                    ## "sigma_s",
                                    "zeta", "phi")) |>
  mcmc_trace(facet_args = list(labeller = ggplot2::label_parsed))

drm_rec2$stanfit$draws(variables = c("alpha", "sigma_t",
                                    ## "sigma_s",
                                    "zeta", "phi")) |>
  mcmc_dens_overlay(facet_args = list(labeller = ggplot2::label_parsed))

## estimating initialization with mortality

drm_rec3 <-
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
          .toggles = list(ar_re = "rec",
                          movement = 1,
                          est_surv = 1,
                          minit = 1),
          .priors = list(pr_phi_a = 1, pr_phi_b = .1,
                         pr_alpha_a = 4.2, pr_alpha_b = 5.8,
                         pr_zeta_a = 7, pr_zeta_b = 3))

##--- comparing some priors and posteriors ----

drm_rec3$stanfit$draws(variables = c("phi")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dgamma(x, shape = 1, rate = .1),
                xlim = c(0, 3),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_rec3$stanfit$draws(variables = c("alpha")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dbeta(x,
                                 shape1 = drm_rec3$data$pr_alpha_a,
                                 shape2 = drm_rec3$data$pr_alpha_b),
                xlim = c(0, 1),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_rec3$stanfit$draws(variables = c("zeta")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dbeta(x,
                                 shape1 = drm_rec3$data$pr_zeta_a,
                                 shape2 = drm_rec3$data$pr_zeta_b),
                xlim = c(0, 1),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_rec3$stanfit$draws(variables = c("sigma_t")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dlnorm(x,
                                  meanlog = drm_rec3$data$pr_lsigma_t_mu,
                                  sdlog = drm_rec3$data$pr_lsigma_t_sd),
                xlim = c(0, .5),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

##--- Viz relationships ----

## * make this into a function!
## ** that is far from easy!

## recruitment

newdata_rec <- data.frame(c_stemp =
                            seq(from = quantile(dat_train$c_stemp, .05),
                                to = quantile(dat_train$c_stemp, .95),
                                length.out = 200))

rec_samples <- marg_rec(drm_rec3, newdata_rec)

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
  mutate(model = "drm_rec3")

rec_fig3 <-
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

rec_fig3

##--- * temperature that optimizes recruitment ----

betas_rec <-
  drm_rec3$stanfit$draws(variables = "beta_r",
                         format = "matrix")

max_quad_x(betas_rec[, 2], betas_rec[, 3],
           offset = avgs["stemp"]) |>
  apply(2, quantile, probs = c(.1, .5, .9))

##--- Convergence & estimates ----

## all rhat's look good (no larger than 1.01)
drm_rec3$stanfit$summary(variables = c("beta_r", "beta_t",
                                       "alpha", "sigma_t",
                                       ## "sigma_s",
                                       "zeta", "phi"))

## the different chains are in agreement and converging.
drm_rec3$stanfit$draws(variables = c("beta_r", "beta_t")) |>
  mcmc_trace()

drm_rec3$stanfit$draws(variables = c("beta_r", "beta_t")) |>
  mcmc_dens_overlay()

drm_rec3$stanfit$draws(variables = c("alpha", "sigma_t",
                                    ## "sigma_s",
                                    "zeta", "phi")) |>
  mcmc_trace(facet_args = list(labeller = ggplot2::label_parsed))

drm_rec3$stanfit$draws(variables = c("alpha", "sigma_t",
                                    ## "sigma_s",
                                    "zeta", "phi")) |>
  mcmc_dens_overlay(facet_args = list(labeller = ggplot2::label_parsed))

## checking figs

rec_fig | rec_fig2 | rec_fig3

##--- * model comparison ----

loos <- list("init_data" = drm_rec$stanfit$loo(),
             "init_est"  = drm_rec3$stanfit$loo(),
             "init_mort" = drm_rec3$stanfit$loo())

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

forecast_rec2 <- predict_drm(drm = drm_rec2,
                             new_data = dat_test,
                             past_data = filter(dat_train,
                                                year == max(year)),
                             seed = 125,
                             cores = 4)

forecast_rec3 <- predict_drm(drm = drm_rec3,
                             new_data = dat_test,
                             past_data = filter(dat_train,
                                                year == max(year)),
                             seed = 125,
                             cores = 4)

##--- proc ----

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

##--- f2 ----

##--- Viz predicted and observed ----

fitted_rec2 <-
  drm_rec2$stanfit$draws(variables = "y_pp",
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

forecast_rec2 <-
  forecast_rec2$draws(variables = "y_proj",
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

forecast_rec2 <-
  dat_test |>
  select(year, patch, lat_floor, dens) |>
  mutate(pair = row_number(), .before = 1) |>
  left_join(forecast_rec2, by = "pair") |>
  select(- pair)

fitted_rec2 <-
  dat_train |>
  select(year, patch, lat_floor, dens) |>
  bind_cols(select(fitted_rec2, -pair))

##--- f3 ----

##--- Viz predicted and observed ----

fitted_rec3 <-
  drm_rec3$stanfit$draws(variables = "y_pp",
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

forecast_rec3 <-
  forecast_rec3$draws(variables = "y_proj",
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

forecast_rec3 <-
  dat_test |>
  select(year, patch, lat_floor, dens) |>
  mutate(pair = row_number(), .before = 1) |>
  left_join(forecast_rec3, by = "pair") |>
  select(- pair)

fitted_rec3 <-
  dat_train |>
  select(year, patch, lat_floor, dens) |>
  bind_cols(select(fitted_rec3, -pair))

##--- comparing!! ----

bind_rows(
    bind_rows(fitted_rec, forecast_rec) |>
    mutate(model = "data"),
    bind_rows(fitted_rec3, forecast_rec3) |>
    mutate(model = "est"),
    bind_rows(fitted_rec3, forecast_rec3) |>
    mutate(model = "mort")
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

aux_qt <-
  data.frame(Model = names(loos),
             delta_LOOIC = loos_out[order(rownames(loos_out)), 1],
             LOOIC = sapply(loos, \(x) x$estimates[3, 1])) |>
  mutate(Model = case_when(Model == "init_data" ~ "data",
                           Model == "init_est"  ~ "est",
                           TRUE                 ~ "mort"))

bind_rows(
    forecast_rec |>
    mutate(model = "data"),
    forecast_rec2 |>
    mutate(model = "est"),
    forecast_rec3 |>
    mutate(model = "mort")
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

##--- comparing estimates ----

drm_rec2$stanfit$summary(variables = c("beta_r", "beta_t",
                                      "alpha", "sigma_t",
                                      "zeta", "phi"))

drm_rec3$stanfit$summary(variables = c("beta_r", "beta_t",
                                       "alpha", "sigma_t",
                                       "zeta", "phi"))
