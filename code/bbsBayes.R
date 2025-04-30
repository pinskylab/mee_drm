library(bbsBayes2)
library(sf)
library(dplyr)
library(ggplot2)
library(patchwork)
library(drmr)

##--- Mallard data in Maine ----

fetch_bbs_data(level = "stop")

s <- stratify(by = "latlong", species = "Baltimore oriole")

p <- prepare_data(s, min_n_routes = 1)

p$raw_data |>
  count(year, strata_name) |>
  arrange(- n)

my_dt <-
  p$raw_data |>
  filter(country_num == "840") |>
  filter(state == "MAINE")


my_dt |>
  st_as_sf(coords = c("longitude", "latitude"),
           crs = st_crs("epsg:4326")) |>
  select(count, starts_with("n_")) |>
  plot()

##--- Maine Map -----

states <-
  rnaturalearth::ne_states(
                     country = c("United States of America"),
                     returnclass = "sf"
                 ) |>
  st_transform("epsg:4326")

states <- states |>
  filter(name == "Maine") |>
  st_geometry() |>
  st_simplify(dTolerance = 2000) |>
  st_cast("POLYGON")

states <- states[4]

states <- st_make_valid(states) |>
  nngeo::st_remove_holes()

##--- creating grid ----

my_grid <- st_make_grid(x = states, n = c(10, 10))

filter_grid <-
  st_intersects(my_grid, st_intersection(st_buffer(states, 10)),
                sparse = FALSE)

my_grid <- my_grid[filter_grid[, 1], ]

##--- counting events per cell ----

my_grid <- my_grid |>
  st_as_sf() |>
  mutate(id = row_number())

st_geometry(my_grid) <- "geometry"

my_map <- my_grid

year_id <- expand.grid(id = my_grid$id,
                       year = sort(unique(my_dt$year)))

my_grid <- st_join(my_grid,
                   st_as_sf(my_dt,
                            coords = c("longitude",
                                       "latitude"),
                            crs = st_crs("epsg:4326")),
                   st_contains) |>
  filter(!is.na(year)) |>
  st_drop_geometry()

my_grid <- my_grid |>
  group_by(id, year,
           .drop = FALSE) |>
  summarise(count = sum(count, na.rm = TRUE),
            n_routes = sum(n_routes, na.rm = TRUE),
            n_obs = sum(n_obs, na.rm = TRUE),
            n_obs_sites = sum(n_obs_sites, na.rm = TRUE),
            non_zero_weight = sum(non_zero_weight, na.rm = TRUE)) |>
  ungroup() |>
  mutate()

my_grid <- year_id |>
  filter(year != 2020) |>
  left_join(my_grid, by = c("id", "year")) |>
  as_tibble()

my_grid <- my_grid |>
  st_as_sf()

my_grid <- my_grid |>
  mutate(across(count:non_zero_weight,
                .fns = \(x) tidyr::replace_na(x, 0)))

my_grid <- left_join(my_map, my_grid,
                     by = c("id"))

my_grid <- my_grid |>
  mutate(area = units::drop_units(units::set_units(st_area(geometry), "km^2")))

my_grid <- my_grid |>
  mutate(y = count / area)

my_grid |>
  filter(year == "2021") |>
  select(y) |>
  plot()

my_grid <- my_grid |>
  filter(year > 1990)

##--- drmr model ----

