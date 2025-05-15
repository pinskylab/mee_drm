library(bbsBayes2)
library(dplyr)
library(sf)
library(stars)
library(ggplot2)

year <- 1980

xx <- read_stars("data/chelsa/tas_2019.tif")

##--- bird data ----

fetch_bbs_data(level = "stop")

species <- "Red-bellied woodpecker"

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

st_write(c_hull_pts, dsn = "data/birds/shape/window.shp")

my_dt <-
  read_stars("data/chelsa/tas_01_1980.tif")

my_dt <- my_dt[c_hull_pts]

plot(my_dt)
