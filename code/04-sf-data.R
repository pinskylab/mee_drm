library(dplyr)
library(ggplot2)
library(bayesplot)
library(sf)
library(drmr)
library(cmdstanr)

## loading data
data(sum_fl)

## loading map
map_name <- system.file("maps/sum_fl.shp", package = "drmr")

polygons <- st_read(map_name)

## instantaneous fishing mortality rates
fmat <-
  system.file("fmat.rds", package = "drmr") |>
  readRDS()

##--- splitting data for validation ----

## reserving 5 years for forecast assessment
first_year_forecast <- max(sum_fl$year) - 5

## "year to id"
first_id_forecast <-
  first_year_forecast - min(sum_fl$year) + 1

years_all <- seq_len(NCOL(fmat))
years_train <- years_all[years_all < first_id_forecast]
years_test <- years_all[years_all >= first_id_forecast]

f_train <- fmat[, years_train]
f_test  <- fmat[, years_test]

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

##--- using function to prepare data for stan ----

form_m <- ~ 1 # + btemp + I(btemp * btemp)
form_r <- ~ 1 +
  c_stemp +
  I(c_stemp * c_stemp) ## +
## as.factor(patch)
form_t <- ~ 1 +
  c_hauls + 
  c_btemp +
  I(c_btemp * c_btemp) +
  c_stemp +
  I(c_stemp * c_stemp) +
  I(c_btemp * c_stemp) ## +
## as.factor(patch)

x_m <- model.matrix(form_m,
                    data = dat_train)
x_r <- model.matrix(form_r,
                    data = dat_train)

x_t <- model.matrix(form_t,
                    data = dat_train)

## adj_mat <- gen_adj(st_buffer(st_geometry(filter(polygons,
##                                                 patch <= 8)),
adj_mat <- gen_adj(st_buffer(st_geometry(polygons),
                             dist = 2500))

## row-standardized matrix
adj_mat <-
  t(apply(adj_mat, 1, \(x) x / (sum(x))))

chains <- 4
cores <- 4

##--- fitting SDM (with movement and F) ----

mcmc_drm <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "patch", ## vector of patches
          f_mort = f_train,
          n_ages = NROW(f_train), ## number of age groups
          family = "gamma",
          seed = 2025,
          ## adj_mat = adj_mat,
          formula_zero = form_t,
          formula_rec = form_r,
          formula_surv = form_m,
          ## ages_movement = 3, ## age at which fish are able to move
          iter_sampling = 1000, ## number of samples after warmup
          iter_warmup = 1000, ## number of warmup samples
          parallel_chains = cores,
          init = "pathfinder",
          chains = chains,
          .toggles = list(time_ar = 1,
                          movement = 0,
                          est_mort = 0),
          .priors = list(pr_alpha_a = 5, pr_alpha_b = 5))

##--- * MCMC diagnostics ----

viz_pars <-
  c("tau[1]",
    "alpha[1]",
    "phi[1]")

mcmc_trace(mcmc_drm$draws$draws(variables = viz_pars))
mcmc_dens_overlay(mcmc_drm$draws$draws(variables = viz_pars))

mcmc_trace(mcmc_drm$draws$draws(variables = c("coef_r")))
mcmc_dens_overlay(mcmc_drm$draws$draws(variables = c("coef_r")))

mcmc_trace(mcmc_drm$draws$draws(variables = c("coef_t")))
mcmc_dens_overlay(mcmc_drm$draws$draws(variables = c("coef_t")))

##--- comparing models ----

loo_drm <- mcmc_drm$draws$loo()

##--- estimated relationship between recruitment and temperature ----

## recruitment

rec_coefs <-
  mcmc_drm$draws$draws(variables =
                         c("coef_r"),
                       format = "draws_matrix")

y_max <- max_quad_x(beta1 = rec_coefs[, 2],
                    beta2 = rec_coefs[, 3],
                    offset = avgs["stemp"])

quantile(y_max, probs = c(.05, .5, .95))

rec_coefs <- rec_coefs |>
  as_tibble() |>
  mutate(beta_1 = fix_linbeta(`coef_r[2]`, `coef_r[3]`,
                              offset = avgs["stemp"]),
         beta_2 = `coef_r[3]`) |>
  select(beta_1:beta_2) |>
  posterior::as_draws_matrix()

temps <- seq(from = quantile(dat_train$c_stemp + avgs["stemp"], .05),
             to = quantile(dat_train$c_stemp + avgs["stemp"], .95),
             length.out = 200)
x_aux <- model.matrix(~ 0 + temps + I(temps^2))

dim(x_aux)
dim(rec_coefs)

rec_estimates <-
  tcrossprod(rec_coefs, x_aux)

rec_estimates <- exp(c(rec_estimates))

rec_estimates <-
  tibble(temperature = c(sapply(temps, \(x) rep(x, nrow(rec_coefs)))),
         mortality   = rec_estimates)

to_plot <-
  rec_estimates |>
  group_by(temperature) |>
  summarise(ll = quantile(mortality, probs = .05),
            l = quantile(mortality, probs = .1),
            m = median(mortality),
            u = quantile(mortality, probs = .9),
            uu = quantile(mortality, probs = .95)) |>
  ungroup()

to_plot2 <- to_plot |>
  filter(temperature <= quantile(y_max, probs = .75),
         temperature >= quantile(y_max, probs = .25))

to_plot |>
  ggplot(data = ,
         aes(x = temperature,
             y = m)) +
  ## geom_ribbon(aes(ymin = l, ymax = u),
  ##             fill = 2,
  ##             alpha = .2) +
  ## geom_ribbon(aes(ymin = ll, ymax = uu),
  ##             fill = 2,
  ##             alpha = .2) +
  geom_area(data = to_plot2,
            fill = 4, alpha = .3) +
  annotate("segment",
           x = quantile(y_max, probs = .5),
           y = 0, yend = max(to_plot$m),
           col = 4, lty = 2, lwd = 1.2) +
  geom_line() +
  labs(y = "Recruitment",
       x = "Temperature") +
  theme_bw()


##--- forecasting ----

##--- * DRM ----

x_tt <- model.matrix(form_t,
                     data = dat_test)
x_mt <- model.matrix(form_m,
                     data = dat_test)
x_rt <- model.matrix(form_r,
                     data = dat_test)

x_mpast <-
  model.matrix(form_m,
               data = filter(dat_train, year == max(year)))

forecast_drm <- predict_drm(drm = mcmc_drm$draws,
                            drm_data = mcmc_drm$data,
                            ntime_for =
                              length(unique(dat_test$year)),
                            x_tt = x_tt,
                            x_rt = x_rt,
                            f_test = f_test[, -1],
                            seed = 125,
                            cores = 4)
