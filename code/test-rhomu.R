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

init_ldens <-
  dat_train |>
  filter(year < 1987) |>
  summarise(dens = sum(area_km2 * dens) / sum(area_km2)) |>
  pull(dens)

init_ldens <-
  init_ldens *
  seq(.9, .1, len = NROW(f_train)) / sum(seq(.9, .1, len = NROW(f_train)))

init_ldens <- log(init_ldens[-1])

drm_rec_old <-
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
          ## init_data = init_ldens,
          adj_mat = adj_mat, ## A matrix for movement routine
          ages_movement = c(0, 0,
                            rep(1, 12),
                            0, 0), ## ages allowed to move
          .toggles = list(ar_re = "rec",
                          rho_mu = 0,
                          movement = 1,
                          est_surv = 1,
                          est_init = 0,
                          minit = 1),
          .priors = list(pr_phi_a = 1, pr_phi_b = .1,
                         pr_alpha_a = 4.2, pr_alpha_b = 5.8,
                         pr_zeta_a = 7, pr_zeta_b = 3))

drm_rec_new <-
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
          ## init_data = init_ldens,
          adj_mat = adj_mat, ## A matrix for movement routine
          ages_movement = c(0, 0,
                            rep(1, 12),
                            0, 0), ## ages allowed to move
          .toggles = list(ar_re = "rec",
                          rho_mu = 1,
                          movement = 1,
                          est_surv = 1,
                          est_init = 0,
                          minit = 1),
          .priors = list(pr_phi_a = 1, pr_phi_b = .1,
                         pr_alpha_a = 4.2, pr_alpha_b = 5.8,
                         pr_zeta_a = 7, pr_zeta_b = 3))

##--- comparing some priors and posteriors ----

drm_rec_old$stanfit$draws(variables = c("phi")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dgamma(x, shape = 1, rate = .1),
                xlim = c(0, 3),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_rec_new$stanfit$draws(variables = c("phi")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dgamma(x, shape = 1, rate = .1),
                xlim = c(0, 3),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_rec_new$stanfit$draws(variables = c("xi")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) {
    y <- - x
    dnorm(log(y),
          mean = drm_rec_new$data$pr_lmxi_mu,
          sd = drm_rec_new$data$pr_lmxi_sd) / y
  },
  xlim = c(-4, -1e-16),
  n = 501,
  inherit.aes = FALSE,
  color = 2,
  lwd = 1.2)

drm_rec_new$stanfit$draws(variables = c("xi")) |>
  mcmc_trace()

drm_rec_old$stanfit$draws(variables = c("alpha")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dbeta(x,
                                 shape1 = drm_rec_old$data$pr_alpha_a,
                                 shape2 = drm_rec_old$data$pr_alpha_b),
                xlim = c(0, 1),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_rec_new$stanfit$draws(variables = c("alpha")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dbeta(x,
                                 shape1 = drm_rec_new$data$pr_alpha_a,
                                 shape2 = drm_rec_new$data$pr_alpha_b),
                xlim = c(0, 1),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_rec_old$stanfit$draws(variables = c("zeta")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dbeta(x,
                                 shape1 = drm_rec_old$data$pr_zeta_a,
                                 shape2 = drm_rec_old$data$pr_zeta_b),
                xlim = c(0, 1),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_rec_new$stanfit$draws(variables = c("zeta")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dbeta(x,
                                 shape1 = drm_rec_new$data$pr_zeta_a,
                                 shape2 = drm_rec_new$data$pr_zeta_b),
                xlim = c(0, 1),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_rec_old$stanfit$draws(variables = c("sigma_t")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dlnorm(x,
                                  meanlog = drm_rec_old$data$pr_lsigma_t_mu,
                                  sdlog = drm_rec_old$data$pr_lsigma_t_sd),
                xlim = c(0, .5),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

drm_rec_new$stanfit$draws(variables = c("sigma_t")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dlnorm(x,
                                  meanlog = drm_rec_new$data$pr_lsigma_t_mu,
                                  sdlog = drm_rec_new$data$pr_lsigma_t_sd),
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

rec_samples <- marg_rec(drm_rec_old, newdata_rec)

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

rec_fold <-
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

rec_samples <- marg_rec(drm_rec_new, newdata_rec)

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

rec_fnew <-
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


rec_fold | rec_fnew

##--- pzero vs density ----

mu_range <- drm_rec_new$stanfit$summary(variables = "mu")
mu_range <- c(min(mu_range$q5), max(mu_range$q95))

mu_range <- seq(from = mu_range[1],
                to = mu_range[2], len = 100)

xis <- drm_rec_new$stanfit$draws("xi", format = "matrix")

beta_0 <- drm_rec_new$stanfit$draws("beta_t[1]", format = "matrix")

aux_mu <-
  sapply(mu_range, \(x) {
    rho <- plogis(as.numeric(beta_0) + as.numeric(xis) * x)
    data.frame(mu = x, rho = rho)
  }, simplify = FALSE) |>
  bind_rows()

summary_mu <-
  group_by(aux_mu, mu) |>
  summarise(rho_m = median(rho),
            rho_l = quantile(rho, .05),
            rho_u = quantile(rho, .95)) |>
  ungroup()

ggplot(data = summary_mu,
       aes(x = mu,
           y = rho_m)) +
  geom_ribbon(aes(ymin = rho_l,
                  ymax = rho_u)) +
  geom_line() +
  theme_bw()

##--- * temperature that optimizes recruitment ----

betas_rold <-
  drm_rec_old$stanfit$draws(variables = "beta_r",
                            format = "matrix")

betas_rnew <-
  drm_rec_new$stanfit$draws(variables = "beta_r",
                            format = "matrix")

cbind(
    "OLD" = max_quad_x(betas_rold[, 2], betas_rold[, 3],
                       offset = avgs["stemp"]) |>
    apply(2, quantile, probs = c(.1, .5, .9)),
    "NEW" = max_quad_x(betas_rnew[, 2], betas_rnew[, 3],
               offset = avgs["stemp"]) |>
    apply(2, quantile, probs = c(.1, .5, .9))
)
    
##--- Convergence & estimates ----

## all rhat's look good (no larger than 1.01)
drm_rec_old$stanfit$summary(variables = c("beta_r", "beta_t",
                                          "alpha", "sigma_t",
                                          ## "sigma_s",
                                          "zeta", "phi"))
drm_rec_new$stanfit$summary(variables = c("beta_r", "beta_t",
                                          "alpha", "sigma_t",
                                          "xi",
                                          "zeta", "phi"))

rec_fold + rec_fnew +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10))

##--- forecasting ----

##--- * DRM ----

forecast_rold <- predict_drm(drm = drm_rec_old,
                             new_data = dat_test,
                             past_data = filter(dat_train,
                                                year == max(year)),
                             seed = 125,
                             f_test = f_test,
                             cores = 4)

forecast_rnew <- predict_drm(drm = drm_rec_new,
                             new_data = dat_test,
                             past_data = filter(dat_train,
                                                year == max(year)),
                             seed = 125,
                             f_test = f_test,
                             cores = 4)

##--- Viz predicted and observed ----

fitted_rec_old <-
  drm_rec_old$stanfit$draws(variables = "y_pp",
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

for_rec_old <-
  forecast_rold$draws(variables = "y_proj",
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

for_rec_old <-
  dat_test |>
  select(year, patch, lat_floor, dens) |>
  mutate(pair = row_number(), .before = 1) |>
  left_join(for_rec_old, by = "pair") |>
  select(- pair)

fitted_rec_old <-
  dat_train |>
  select(year, patch, lat_floor, dens) |>
  bind_cols(select(fitted_rec_old, -pair))

fitted_rec_new <-
  drm_rec_new$stanfit$draws(variables = "y_pp",
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

for_rec_new <-
  forecast_rnew$draws(variables = "y_proj",
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

for_rec_new <-
  dat_test |>
  select(year, patch, lat_floor, dens) |>
  mutate(pair = row_number(), .before = 1) |>
  left_join(for_rec_new, by = "pair") |>
  select(- pair)

fitted_rec_new <-
  dat_train |>
  select(year, patch, lat_floor, dens) |>
  bind_cols(select(fitted_rec_new, -pair))

##--- Figure 1 ----

fitted_all <-
  bind_rows(
      bind_rows(fitted_rec_old, for_rec_old) |>
      mutate(model = "DRM (old)"),
      bind_rows(fitted_rec_new, for_rec_new) |>
      mutate(model = "DRM (new)")
  )

fitted_all |>
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

##--- densities (theoretical means) ----

fitted_rec_old <-
  drm_rec_old$stanfit$draws(variables = "mu",
                        format = "draws_df") |>
  tidyr::pivot_longer(cols = starts_with("mu"),
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

for_rec_old <-
  forecast_rold$draws(variables = "mu_proj",
                     format = "draws_df") |>
  tidyr::pivot_longer(cols = starts_with("mu_proj"),
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

for_rec_old <-
  dat_test |>
  select(year, patch, lat_floor, dens) |>
  mutate(pair = row_number(), .before = 1) |>
  left_join(for_rec_old, by = "pair") |>
  select(- pair)

fitted_rec_old <-
  dat_train |>
  select(year, patch, lat_floor, dens) |>
  bind_cols(select(fitted_rec_old, -pair))

fitted_rec_new <-
  drm_rec_new$stanfit$draws(variables = "mu",
                         format = "draws_df") |>
  tidyr::pivot_longer(cols = starts_with("mu"),
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

for_rec_new <-
  forecast_rnew$draws(variables = "mu_proj",
                      format = "draws_df") |>
  tidyr::pivot_longer(cols = starts_with("mu_proj"),
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

for_rec_new <-
  dat_test |>
  select(year, patch, lat_floor, dens) |>
  mutate(pair = row_number(), .before = 1) |>
  left_join(for_rec_new, by = "pair") |>
  select(- pair)

fitted_rec_new <-
  dat_train |>
  select(year, patch, lat_floor, dens) |>
  bind_cols(select(fitted_rec_new, -pair))

##--- Figure 1 ----

theoretical_all <-
  bind_rows(
      bind_rows(fitted_rec_old, for_rec_old) |>
      mutate(model = "DRM (old)"),
      bind_rows(fitted_rec_new, for_rec_new) |>
      mutate(model = "DRM (new)")
  )

theoretical_all |>
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
  scale_y_continuous(breaks = scales::trans_breaks(identity,
                                                   identity,
                                                   n = 3),
                     trans = "log1p") +
  theme_bw() +
  guides(color = "none",
         fill = "none") +
  labs(y = "Density (in hundreds of individuals per square-km)",
       x = "Year") +
  theme(strip.background = element_rect(fill = "white"))

##--- relationship between mu and rho ----

rho_mu <-
  drm_rec_old$stanfit$summary(variables = c("y_pp",
                                            "mu", "rho")) |>
  mutate(pair = gsub("\\D", "", variable),
         parameter = gsub("\\W", "", substr(variable, 1, 3)),
         .before = variable) |>
  mutate(pair = as.integer(pair)) |>
  arrange(pair)

rho_mu |>
  select(pair, parameter, median, q5, q95) |>
  tidyr::pivot_wider(names_from = parameter,
                     values_from = median:q95) |>
  ggplot(data = _,
         aes(x = median_rho,
             y = median_mu)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1,
              col = 2) +
  ## geom_linerange(aes(ymin = q5_mu,
  ##                    ymax = q95_mu)) +
  ## geom_linerange(aes(xmin = q5_rho,
  ##                    xmax = q95_rho)) +
  theme_bw() +
  labs(y = expression(hat(y)),
       x = expression(mu))

rho_mu |>
  select(pair, parameter, median, q5, q95) |>
  tidyr::pivot_wider(names_from = parameter,
                     values_from = median:q95) |>
  ggplot(data = _,
         aes(x = median_mu,
             y = median_y_p)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1,
              col = 2) +
  ## geom_linerange(aes(ymin = q5_mu,
  ##                    ymax = q95_mu)) +
  ## geom_linerange(aes(xmin = q5_rho,
  ##                    xmax = q95_rho)) +
  theme_bw() +
  labs(y = expression(hat(y)),
       x = expression(mu))
