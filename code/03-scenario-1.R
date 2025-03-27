library(sf)
library(dplyr)
library(ggplot2)
library(arrow)
library(geoarrow)
library(drmr)

##--- loading "base" dataset ----

my_dt <- open_dataset("data/basemap.parquet") |>
  st_as_sf()

## number of timepoints
n_time <- my_dt |>
  pull(year) |>
  unique() |>
  length()

##--- simulating data for scenario 1 ----

## setting temperature with maximum recruitment to mean(my_dt$avg_sst) + t_max
sim_tempt <- my_dt$cs_avg_sst
t_max <- 0.08
br_0 <- -2
br_2 <- -0.25
br_1 <- - 2 * t_max * br_2

betas_r <- c(br_0, br_1, br_2)

n_ages <- 10

form_m <- ~ 1 ## + btemp + I(btemp * btemp)
form_r <- ~ 1 + stemp + I(stemp * stemp)
form_t <- ~ 1 +
  btemp +
  I(btemp * btemp) +
  stemp +
  I(stemp * stemp) +
  I(btemp * stemp)

x_m <- model.matrix(form_m,
                    data = dat_train)
x_r <- model.matrix(form_r,
                    data = dat_train)

x_t <- model.matrix(form_t,
                    data = dat_train)

data0 <-
  make_data(y = my_dt$min_sst,
            time = dat_train$year,
            site = dat_train$id,
            n_ages = n_ages,
            ## f_mort = f_train,
            ## m = 0.25,
            age_selectivity = rep(1, n_ages),
            ## ages_movement = 3,
            ## adj_mat = adj_mat,
            x_t = x_t,
            x_r = x_r,
            x_m = x_m,
            init_data = rep(1, n_ages - 1),
            family = "gamma",
            .priors = list(pr_phi_mu = -1,
                           pr_phi_sd = .5,
                           pr_coef_r_mu = rep(0, 3),
                           pr_coef_r_sd = c(1, .25, .25)),
            .toggles = list(cloglog = 0,
                            movement = 0,
                            est_mort = 1,
                            time_ar = 1,
                            est_init = 0,
                            qr_t = 0,
                            qr_r = 0,
                            qr_m = 0))

pars <- list("alpha" = alpha_r,
             "tau" = tau,
             "coef_m" = betas_m,
             "coef_r" = betas_r,
             "coef_t" = betas_t)

scen_1 <-
  model_sim(dat = data0,
            model = "drm",
            selectivity = data0$age_selectivity,
            ar_time = TRUE)
