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

##--- DRM Survival ----

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
          .toggles = list(ar_re = "dens",
                          sp_re = "dens",
                          est_surv = 1,
                          movement = 1),
          .priors = list(pr_phi_a = 1, pr_phi_b = .1,
                         pr_alpha_a = 4.2, pr_alpha_b = 5.8,
                         pr_zeta_a = 7, pr_zeta_b = 3))

##--- Convergence & estimates ----

## r_hat for beta_s indicate convergence issues
drm_surv$stanfit$summary(variables = c("beta_s", "beta_t",
                                       "alpha", "sigma_t",
                                       "sigma_s",
                                       "zeta", "phi"))

## the different chains are in agreement and converging.
drm_surv$stanfit$draws(variables = c("beta_s", "beta_t")) |>
  mcmc_trace()

drm_surv$stanfit$draws(variables = c("beta_s", "beta_t")) |>
  mcmc_dens_overlay()

drm_surv$stanfit$draws(variables = c("alpha", "sigma_s",
                                     "sigma_t",
                                     "zeta", "phi")) |>
  mcmc_trace(facet_args = list(labeller = ggplot2::label_parsed))

drm_surv$stanfit$draws(variables = c("alpha", "sigma_s",
                                     "sigma_t",
                                     "zeta", "phi")) |>
  mcmc_dens_overlay(facet_args = list(labeller = ggplot2::label_parsed))

## the traceplots suggests we may need a reparametrization for the betas or a
## longer warmup period. The former can be achieved through a qr
## reparametrization
## (https://mc-stan.org/docs/stan-users-guide/regression.html#QR-reparameterization.section)

X_aux <- model.matrix(update(drm_surv$formulas$formula_surv, ~. -1 + 0),
                      data = dat_train)

qr_x <- qr(X_aux)

q_qrx <- qr.Q(qr_x)
q_ast <- q_qrx * sqrt(NROW(X_aux) - 1)
r_qrx <- qr.R(qr_x)
r_ast <- r_qrx / sqrt(NROW(X_aux) - 1)
r_inv <- solve(r_ast)

dat_train <-
  mutate(dat_train, aux_x1 = q_ast[, 1],
         aux_x2 = q_ast[, 2])

drm_surv <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "patch",
          family = "gamma",
          seed = 202505,
          formula_zero = ~ 1 + c_hauls + c_btemp + c_stemp,
          formula_rec = ~ 1,
          formula_surv = ~ 1 + aux_x1 + aux_x2,
          f_mort = f_train,
          n_ages = NROW(f_train),
          adj_mat = adj_mat, ## A matrix for movement routine
          ages_movement = c(0, 0,
                            rep(1, 12),
                            0, 0), ## ages allowed to move
          .toggles = list(ar_re = "dens",
                          sp_re = "dens",
                          est_surv = 1,
                          movement = 1),
          .priors = list(pr_phi_a = 1, pr_phi_b = .1,
                         pr_alpha_a = 4.2, pr_alpha_b = 5.8,
                         pr_zeta_a = 7, pr_zeta_b = 3))


##--- convergence again ----

## now all rhat's look good (no larger than 1.01)
drm_surv$stanfit$summary(variables = c("beta_s", "beta_t",
                                       "alpha", "sigma_t",
                                       "sigma_s",
                                       "zeta", "phi"))

## the different chains are in agreement and converging.
drm_surv$stanfit$draws(variables = c("beta_s", "beta_t")) |>
  mcmc_trace()

drm_surv$stanfit$draws(variables = c("beta_s", "beta_t")) |>
  mcmc_dens_overlay()

drm_surv$stanfit$draws(variables = c("alpha", "sigma_s",
                                     "sigma_t",
                                     "zeta", "phi")) |>
  mcmc_trace(facet_args = list(labeller = ggplot2::label_parsed))

drm_surv$stanfit$draws(variables = c("alpha", "sigma_s",
                                     "sigma_t",
                                     "zeta", "phi")) |>
  mcmc_dens_overlay(facet_args = list(labeller = ggplot2::label_parsed))
