library(bbsBayes2)
library(sf)
library(dplyr)
library(arrow)
library(geoarrow)

##--- selecting species ----

sps <- readr::read_csv("data/birds/martins2024significant.csv")

fetch_bbs_data(level = "stop")

species <- "Red-bellied Woodpecker"

s <- stratify(by = "latlong", species = species)

p <- prepare_data(s, min_n_routes = 1)

p$raw_data |>
  count(year, strata_name) |>
  arrange(- n)

my_dt <-
  p$raw_data |>
  filter(country_num == "840")

my_dt |>
  st_as_sf(coords = c("longitude", "latitude"),
           crs = st_crs("epsg:4326")) |>
  select(count, starts_with("n_")) |>
  st_geometry() |>
  plot()

## convex hull of the observations (to define study region)
c_hull_pts <-
  my_dt |>
  st_as_sf(coords = c("longitude", "latitude"),
           crs = st_crs("epsg:4326")) |>
  st_geometry() |>
  st_union() |>
  st_convex_hull()

## this region will be used to crop the environmental data
st_write(c_hull_pts, dsn = "data/birds/shape/window.shp")

##--- Maine Map -----

states <-
  rnaturalearth::ne_states(
                     country = c("United States of America"),
                     returnclass = "sf"
                 ) |>
  st_transform("epsg:4326")

states <-
  states[st_intersects(states, c_hull_pts,
                       sparse = FALSE)[, 1], ] |>
  ## st_geometry() |>
  st_simplify(dTolerance = 2000) |>
  select(iso_a2, code_local, region, name)

states_merged <-
  states |>
  st_geometry() |>
  st_cast("POLYGON") |>
  st_union() |>
  st_make_valid() |>
  nngeo::st_remove_holes() |>
  st_simplify(dTolerance = 5000)

## climate: https://chelsa-climate.org/

##--- creating grid ----

my_grid <- st_make_grid(x = c_hull_pts, n = c(30, 20))

filter_grid <-
  st_intersects(my_grid,
                c_hull_pts,
                sparse = FALSE)

my_grid <- my_grid[filter_grid[, 1], ]

plot(states_merged)
plot(my_grid, col = scales::alpha(2, .2),
     add = TRUE)

filter_grid <-
  st_intersects(my_grid,
                states_merged,
                sparse = FALSE)


plot(my_grid[filter_grid[, 1], ],
     col = scales::alpha(4, .2),
     add = TRUE)

my_grid <- my_grid[filter_grid[, 1]]

##--- counting events per cell ----

my_grid <- my_grid |>
  st_as_sf() |>
  mutate(id = row_number(),
         .before = x)

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
            n_obs_sites = sum(n_obs_sites, na.rm = TRUE),
            mean_obs = mean(mean_obs, na.rm = TRUE),
            non_zero_weight = sum(non_zero_weight, na.rm = TRUE)) |>
  ungroup() |>
  mutate()

my_grid <- year_id |>
  filter(between(year, left = 1980, right = 2019)) |>
  left_join(my_grid, by = c("id", "year")) |>
  as_tibble()

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
  filter(year == "2019") |>
  select(y) |>
  plot()

write_dataset(my_grid, "data/birds/no_env.parquet")
