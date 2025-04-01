library(sf)
library(dplyr)
library(ggplot2)
library(arrow)
library(geoarrow)
library(drmr)
library(bayesplot)

##--- loading "base" dataset ----

my_dt <- open_dataset("data/basemap.parquet") |>
  st_as_sf()

##--- lat band instead ----

my_dt <- my_dt |>
  mutate(ct_y = round(st_coordinates(st_centroid(geometry))[, 2], 1)) |>
  group_by(year, time, ct_y) |>
  summarise(id = min(id),
            across(min_sal:cs_med_vo, mean),
            geometry = st_union(geometry)) |>
  ungroup()

## number of timepoints
n_time <- my_dt |>
  pull(year) |>
  unique() |>
  length()

my_dt <- my_dt |>
  group_by(year) |>
  mutate(id = rank(id),
         .before = id) |>
  ungroup()

##--- simulating data for scenario 1 ----

## setting temperature with maximum recruitment to mean(my_dt$avg_sst) + t_max
sim_tempt <- my_dt$cs_avg_sst
t_max <- 0.08
br_0 <- -.5
br_2 <- -4
br_1 <- - 2 * t_max * br_2

betas_r <- c(br_0, br_1, br_2)
betas_t <- c(0, .8, -.2)

n_ages <- 10

form_m <- ~ 1 ## + btemp + I(btemp * btemp)
form_r <- ~ 1 + cs_avg_sst + I(cs_avg_sst * cs_avg_sst)
form_t <- ~ 1 +
  cs_avg_uo +
  I(cs_avg_uo * cs_avg_uo)

x_m <- model.matrix(form_m,
                    data = my_dt)
x_r <- model.matrix(form_r,
                    data = my_dt)
x_t <- model.matrix(form_t,
                    data = my_dt)

data0 <-
  make_data(y = my_dt$min_sst,
            time = my_dt$year,
            site = my_dt$id,
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
                            est_mort = 0,
                            time_ar = 1,
                            est_init = 0,
                            qr_t = 0,
                            qr_r = 0,
                            qr_m = 0))

pars <- list("alpha"  = 0.4,
             "tau"    = 0.4,
             "coef_m" = NA_real_,
             "coef_r" = betas_r,
             "coef_t" = betas_t,
             "phi"    = 1.2)

set.seed(1254)
scen1_dt <-
  model_sim(dat = data0,
            model = "drm",
            selectivity = data0$age_selectivity,
            ar_time = TRUE,
            init_type = data0$est_init,
            pars = pars)


my_dt <- my_dt |>
  mutate(time_date = time, .after = time) |>
  mutate(time = year - min(year) + 1) |>
  left_join(scen1_dt$mu_and_y,
            by = c("id" = "site", "time"))

ggplot(data = my_dt,
       aes(x = avg_sst,
           y = y)) +
  geom_point() +
  theme_bw()

ggplot(data = my_dt,
       aes(x = y)) +
  geom_density() +
  theme_bw()

my_dt |>
  group_by(id) |>
  summarise(nas_y = mean(is.na(y)),
            sst = mean(cs_avg_sst),
            uo = mean(cs_avg_uo)) |>
  ungroup() |>
  ggplot(data = _) +
  geom_sf(aes(fill = nas_y)) +
  scale_fill_viridis_c(option = "H") +
  theme_bw()

##--- splitting data and fitting a DRM ----

train_dt <- my_dt |>
  filter(year <= 2015)

test_dt <- my_dt |>
  filter(year > 2015)

drm_1 <- fit_drm(.data = train_dt,
                 y_col = "y",
                 time_col = "time",
                 site_col = "id",
                 family = "gamma",
                 formula_zero = form_t,
                 formula_rec = form_r,
                 init = "prior",
                 n_ages = n_ages,
                 init_data = rep(1, n_ages - 1),
                 .priors = list(pr_phi_mu = -1,
                                pr_phi_sd = .5,
                                pr_coef_r_mu = rep(0, 3),
                                pr_coef_r_sd = c(1, .25, .25)),
                 .toggles = list(movement = 0,
                                 est_mort = 0,
                                 time_ar = 1,
                                 est_init = 0),
                 seed = 2025)

drm_1$draws$draws(variables = c("phi", "alpha", "tau")) |>
  mcmc_trace()

drm_1$draws$draws(variables = c("phi", "alpha", "tau")) |>
  mcmc_dens_overlay()

truth <-
  Map(\(x, nm) {
    tibble(variable = sprintf("%s[%d]", nm, seq_along(x)),
           truth = x)
  }, pars, names(pars)) |>
  bind_rows() |>
  filter(!is.na(truth))

truth |>
  left_join(drm_1$draws$summary(variables = names(pars)[-3]),
            by = "variable")

##--- forecasting ----

x_tt <- model.matrix(form_t,
                     data = test_dt)

x_mt <- model.matrix(form_m,
                     data = test_dt)
x_rt <- model.matrix(form_r,
                     data = test_dt)

x_mpast <-
  model.matrix(form_m,
               data = filter(train_dt, year == max(year)))

my_pred <- predict_drm(drm = drm_1$draws,
                       drm_data = drm_1$data,
                       ntime_for =
                         length(unique(test_dt$year)),
                       x_tt = x_tt,
                       x_rt = x_rt,
                       f_test = matrix(0,
                                       ncol = length(unique(test_dt$year)),
                                       nrow = n_ages),
                       cores = 4,
                       seed = 2025)

##--- viz forecasts ----

y_pred <-
  my_pred$draws(variables = "y_proj",
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

y_pred <-
  test_dt |>
  select(id, year, time, y) |>
  mutate(pair = row_number()) |>
  left_join(y_pred, by = "pair")

ggplot(data = y_pred,
       aes(x = m, y = y)) +
  geom_point() +
  scale_x_continuous(transform = "log1p") +
  scale_y_continuous(transform = "log1p") +
  theme_bw()

my_color <- "#0465cf"
ggplot(data = y_pred,
         aes(x = year,
             y = m)) +
  geom_ribbon(aes(ymin = ll, ymax = uu),
              fill = my_color,
              alpha = .4) +
  geom_ribbon(aes(ymin = l, ymax = u),
              fill = my_color,
              alpha = .4) +
  geom_line() +
  geom_point() +
  geom_point(aes(x = year, y = y),
             inherit.aes = FALSE,
             color = 2) +
  facet_wrap(id ~ .,
             scales = "free_y") +
  theme_bw() +
  labs(y = expression(mu),
       x = "Year")
