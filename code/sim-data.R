library(dplyr)
library(ggplot2)
library(sf)
library(geoarrow)
library(arrow)
library(drmr)
library(bayesplot)

set.seed(2025)

source("utils/sim-pop.R")

##--- map where the data will be simulated ----

my_map <- readRDS(file = "data/sim/map-temp.rds")

my_dt <- st_drop_geometry(my_map) |>
  mutate(year_id = year - min(year) + 1,
         .before = year) |>
  as_tibble() |>
  arrange(id, year)

my_map <- my_map |>
  select(id) |>
  distinct()

areas <- st_area(my_map) |>
  units::set_units("km^2") |>
  units::set_units(NULL)

areas <- areas / 1000

##--- other helper quantities ----

fmat <-
  system.file("fmat.rds", package = "drmr") |>
  readRDS()

adj_mat <- my_map |>
  st_geometry() |>
  st_buffer(dist = 2500) |>
  drmr::gen_adj()

## row-standardized matrix
adj_mat <-
  t(apply(adj_mat, 1, \(x) x / (sum(x))))

plot(st_geometry(my_map))

##--- processing SST data from Copernicus ----

satellite_sst <-
  open_dataset("data/copernicus/processed.parquet") |>
  st_as_sf()

my_dt <-
  my_map |>
  st_join(satellite_sst,
          st_intersects)

my_dt <- my_dt |>
  st_drop_geometry() |>
  rename("sal" = "so",
         "sst" = "thetao") |>
  mutate(across(sal:vo, \(x) units::set_units(x, NULL))) |>
  tidyr::pivot_longer(cols = sal:vo,
                      names_to = "qty",
                      values_to = "obs")

my_dt <- my_dt |>
  group_by(id, time, qty) |>
  summarise(min = min(obs, na.rm = TRUE),
            max = max(obs, na.rm = TRUE),
            avg = mean(obs, na.rm = TRUE),
            med = median(obs, na.rm = TRUE)) |>
  ungroup() |>
  tidyr::pivot_wider(names_from = qty,
                     values_from = min:med,
                     names_vary = "slowest")

n_time <- my_dt$time |>
  unique() |>
  length()

ggplot(data = my_dt,
       aes(x = time)) +
  geom_ribbon(aes(ymin = min_sst,
                  ymax = max_sst),
              fill = 2, alpha = .4) +
  geom_line(aes(y = avg_sst)) +
  facet_wrap(~ id,
             scales = "free_y") +
  theme_bw()

my_dt |>
  filter(lubridate::year(time) %in% c(1993, 1996, 2000,
                                      2005, 2010, 2021)) |>
  left_join(my_map, y = _, by = "id") |>
  ggplot(data = _) +
  geom_sf(aes(fill = avg_sst)) +
  geom_sf_text(aes(label = id)) +
  facet_wrap(~ time) +
  scale_fill_viridis_c(option = "H") +
  theme_bw()
  
##--- sim recruitment ----

sim_tempt_2 <- my_dt$avg_sst - mean(my_dt$avg_sst)

avg_lrec <- -2

t_max <- - 2
b2 <- -0.25
b1 <- - 2 * t_max * b2
betas <- c(avg_lrec, b1, b2)

tibble(
    temps =
      seq(from = -5, to = 6.5,
          length.out = 200)
) |>
  mutate(rec = exp(1 + b1 * temps + b2 * temps * temps)) |>
  ggplot(data = _,
         aes(x = temps,
             y = rec)) +
  geom_line() +
  theme_bw()

lrec <- cbind(1, c(sim_tempt_2), c(sim_tempt_2) * c(sim_tempt_2)) %*%
  betas

lrec <- as.numeric(lrec)

lrec <- matrix(lrec, ncol = max(my_dt$id), nrow = n_time)

alpha_r <- .4
sigma_r <- .25

ar_rec <- vector(mode = "numeric", length = n_time)

ar_rec[1] <- rnorm(1, sd = sigma_r)

for (i in 2:n_time)
  ar_rec[i] <- alpha_r * ar_rec[i - 1] +
    rnorm(1, sd = sigma_r)

## log areas is to make it "area" dependent!
lrec <- apply(lrec, 2, \(x) x + ar_rec) 

tibble(temp = c(sim_tempt_2),
       rec = exp(c(lrec))) |>
  ggplot(data = _,
         aes(x = temp,
             y = rec)) +
  geom_point(alpha = .5) +
  stat_function(fun = \(x, b0, b1, b2) {
    exp(b0 + b1 * x + b2 * x * x)
  }, args = list(b0 = avg_lrec,
                 b1 = b1, b2 = b2),
  color = 2,
  lwd = 1.2) +
  theme_bw()

##--- simulating lambdas ----

m <- .25

init <- rnorm(n = NROW(fmat) - 1,
              mean = avg_lrec,
              sd = .3) |>
  sort(decreasing = TRUE) |>
  exp()

plot(init, type = "h")
points(init, pch = 19)

fmat2 <- ifelse(fmat > 0, 0, 1)

lambdas <-
  pop_dyn(n_patches = max(my_dt$id),
          n_time = n_time,
          n_ages = NROW(fmat),
          f_a_t = fmat,
          neg_mort = matrix(- m,
                            nrow = n_time,
                            ncol = max(my_dt$id)),
          init = init,
          init_type = 3,
          recruitment = exp(lrec))

dim_lambdas <-
  dim(lambdas)

dimnames(lambdas) <- list("age" = seq_len(dim_lambdas[1]),
                          "year" = sort(unique(lubridate::year(my_dt$time))),
                          "id" = seq_len(dim_lambdas[3]))

##--- * viz lambdas ----

lambdas_df <- array2DF(lambdas,
                       responseName = "density") |>
  mutate(age = as.integer(age),
         year = as.integer(year),
         id = as.integer(id))

patches <- sample(x = seq_len(max(my_dt$id)),
                  size = 5)

times <- c(1, 5, 10, 20, 29)
times <- min(lambdas_df$year) + times - 1

lambdas_df |>
  filter(id %in% patches,
         year %in% times) |>
  ggplot(data = _,
         aes(x = age,
             y = density)) +
  geom_segment(aes(x = age,
                   xend = age,
                   y = 0,
                   yend = density),
               lty = 3) +
  geom_line() +
  geom_point() +
  facet_grid(year ~ id,
             scales = "free_y") +
  theme_bw() +
  labs(y = expression(lambda),
       x = "Age-group")

##--- trying to viz cohorts ----

lambdas_df <-
  lambdas_df |>
  mutate(time = year - min(year) + 1) |>
  group_by(id) |>
  mutate(cohort_id = time - age,
         .before = density) |>
  ungroup()

lambdas_df |>
  filter(id %in% patches,
         cohort_id %in% 11:15) |>
  ggplot(data = _,
         aes(x = age,
             y = density)) +
  geom_segment(aes(x = age,
                   xend = age,
                   y = 0,
                   yend = density),
               lty = 3) +
  geom_line() +
  geom_point() +
  facet_grid(cohort_id ~ id,
             scales = "free_y") +
  theme_bw() +
  labs(y = expression(lambda),
       x = "Age-group")

## do this make sense at all?!

##--- computing expected densities ----

selectivity <-
  dpois(x = 0:15, lambda = 6.5)

selectivity <- selectivity / max(selectivity)

selectivity[1] <- 0

sel_df <- tibble(age_group = seq_along(selectivity),
                 selectivity = selectivity)

areas <- tibble(id = seq_along(areas),
                area = areas)

##--- generating "overall" densities ----

mus_df <- left_join(lambdas_df,
                    sel_df, by = c("age" = "age_group")) |>
  left_join(areas, by = "id")

mus_df <- mus_df |>
  select(- cohort_id) |>
  group_by(id, year, time) |>
  summarise(avg_dens = sum(density),
            avg_dens_sel = as.numeric(crossprod(density, selectivity))) |>
  ungroup()

##--- incorporating to real data ----

my_dt <-
  my_dt |>
  mutate(year = lubridate::year(time),
         .before = time) |>
  select(- time) |>
  left_join(select(mus_df, - time),
            by = c("year", "id")) |>
  relocate(starts_with("avg_dens"),
           .after = year)

ggplot(data = my_dt,
       aes(x = avg_sst,
           y = avg_dens)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  theme_bw()

##--- generating data ----

phi_gamma <- 1.3

my_dt <-
  my_dt |>
  mutate(true_dens =
           sapply(avg_dens,
                  \(x) {
                    rgamma(1, phi_gamma, phi_gamma / x)
                  }),
         true_dens_sel =
           sapply(avg_dens_sel,
                  \(x) {
                    rgamma(1, phi_gamma, phi_gamma / x)
                  }),
         .after = avg_dens_sel)

rho <- .34

my_dt <- my_dt |>
  mutate(absence = rbinom(n = NROW(my_dt),
                          size = 1,
                          prob = rho)) |>
  mutate(y = (1 - absence) * true_dens) |>
  mutate(y_sel = (1 - absence) * true_dens_sel)

##--- plotting y vs temp ----

ggplot(data = my_dt,
       aes(x = avg_sst,
           y = y)) +
  geom_point(alpha = .4) +
  geom_rug(alpha = .4) +
  geom_smooth(se = FALSE) +
  scale_y_continuous(trans = "log1p") +
  theme_bw()

my_dt |>
  filter(year %in% times) |>
  left_join(my_map, y = _, by = "id") |>
  ggplot(data = _) +
  geom_sf(aes(fill = y)) +
  facet_wrap(~ year) +
  scale_fill_viridis_c(option = "H",
                       trans = "log1p") +
  theme_bw()

##--- what if we fit our model? ----

##--- * No init ----

my_dt <- my_dt |>
  mutate(time = year - min(year) + 1,
         std_sst = avg_sst - mean(avg_sst),
         .after = year)

model_fit_noinit <-
  fit_drm(.data = my_dt,
          y_col = "y", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "id", ## vector of patches
          f_mort = fmat[, 1:29],
          n_ages = NROW(fmat), ## number of age groups
          family = "gamma",
          seed = 2025,
          formula_zero = ~ 1,
          formula_rec = ~ std_sst + I(std_sst * std_sst),
          formula_surv = ~ 1,
          iter_sampling = 1000, ## number of samples after warmup
          iter_warmup = 1000, ## number of warmup samples
          parallel_chains = 4,
          init = "pathfinder",
          chains = 4,
          .toggles = list(time_ar = 1,
                          movement = 0,
                          est_mort = 0),
          .priors = list(pr_alpha_a = 5, pr_alpha_b = 5))

viz_pars <-
  c("tau[1]",
    "alpha[1]",
    "phi[1]")

color_scheme_set("mix-red-blue")

mcmc_combo(model_fit_noinit$draws$draws(variables = viz_pars),
           combo = c("trace", "dens_overlay"),
           facet_args = list(labeller = label_parsed))

## true parameters:
true_pars <- c("alpha[1]" = alpha_r,
               "tau[1]"   = sigma_r,
               "phi[1]"   = phi_gamma,
               "coef_r[1]"  = betas[1],
               "coef_r[2]"  = betas[2],
               "coef_r[3]"  = betas[3])

model_fit_noinit$draws$summary(variables = names(true_pars))

##--- * repeating initialization ----

drm_compiled <-
  instantiate::stan_package_model(name = "drm-inits",
                                  package = "drmr")

chains <- 4
cores <- 4

drm_data <-
  model_fit_noinit$data

drm_data$init_type <- 1
drm_data$init_type <- 3
drm_data$est_init <- 0
drm_data$init_data <- init

pf_drm <-
  drm_compiled$pathfinder(data = drm_data,
                          seed = 2025,
                          num_paths = chains,
                          save_single_paths = TRUE,
                          psis_resample = FALSE)

model_fit <-
  drm_compiled$sample(data = drm_data,
                      iter_sampling = 1000,
                      iter_warmup = 1000,
                      seed = 2025,
                      chains = chains,
                      parallel_chains = cores,
                      ## threads_per_chain = 1,
                      init = pf_drm,
                      refresh = 100)

mcmc_combo(model_fit$draws(variables = viz_pars),
           combo = c("trace", "dens_overlay"),
           facet_args = list(labeller = label_parsed))

model_fit$summary(variables = names(true_pars)) |>
  bind_cols(truth = true_pars,
            ..2 = _) |>
  relocate(truth, .before = mean)

##--- comparing estimated to simulated recruitment-env relationship ----

##--- est ----

temps_eval <-
  seq(from = -5, to = 6.5,
      length.out = 200)

rec_coefs <- model_fit$draws(variables = "coef_r",
                             format = "draws_matrix")

y_max <- max_quad_x(beta1 = rec_coefs[, 2],
                    beta2 = rec_coefs[, 3])

## way off
quantile(y_max, probs = c(.05, .5, .95))
t_max

x_aux <- model.matrix(~ 0 + temps_eval + I(temps_eval * temps_eval))

dim(x_aux)
dim(rec_coefs)

rec_estimates <-
  tcrossprod(rec_coefs[, -1], x_aux)

rec_estimates <- exp(c(rec_estimates))

rec_estimates <-
  tibble(temperature = c(sapply(temps_eval, \(x) rep(x, nrow(rec_coefs)))),
         recruitment = rec_estimates)

to_plot <-
  rec_estimates |>
  group_by(temperature) |>
  summarise(ll = quantile(recruitment, probs = .05),
            l = quantile(recruitment, probs = .1),
            m = median(recruitment),
            u = quantile(recruitment, probs = .9),
            uu = quantile(recruitment, probs = .95)) |>
  ungroup() |>
  mutate(truth = exp(b1 * temperature +
                     b2 * temperature * temperature))

to_plot2 <- to_plot |>
  filter(temperature <= quantile(y_max, probs = .75),
         temperature >= quantile(y_max, probs = .25))

to_plot |>
  mutate(temperature = temperature + mean(my_dt$avg_sst)) |>
  ggplot(data = _,
         aes(x = temperature,
             y = m)) +
  geom_ribbon(aes(ymin = l, ymax = u),
              fill = 2,
              alpha = .2) +
  geom_ribbon(aes(ymin = ll, ymax = uu),
              fill = 2,
              alpha = .2) +
  geom_line(color = 2) +
  geom_line(aes(y = truth),
            color = 1,
            lty = 2) +
  labs(y = "Recruitment",
       x = "Temperature") +
  scale_x_continuous() +
  theme_bw()
