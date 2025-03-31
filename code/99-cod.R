library(dplyr)
library(sf)
library(ggplot2)

aux_latlon <- function(geom, lon = TRUE) {
  out <-
    geom |>
    sf::st_coordinates() |>
    apply(2, min)
  ifelse(lon, out[1], out[2])
}

##--- reading data ----

train_path <-
  "https://raw.githubusercontent.com/pinskylab/Approx-Bayes-Comp-Applications/refs/heads/main/data/jude_cod_train.csv"
test_path  <-
  "https://raw.githubusercontent.com/pinskylab/Approx-Bayes-Comp-Applications/refs/heads/main/data/jude_cod_test.csv"

train_dt <- readr::read_csv(train_path)
test_dt <- readr::read_csv(test_path)

train_years <-
  train_dt |>
  pull(year) |>
  range()

test_years <-
  test_dt |>
  pull(year) |>
  range()

##--- formating ----

full_dt <- bind_rows(train_dt, test_dt) |>
  group_by(haulid, region, year, month,
           stratum, lat, lon, sppocean) |>
  summarise(num = sum(NUMLEN, na.rm = TRUE),
            bt = mean(bottemp, na.rm = TRUE),
            depth = mean(depth, na.rm = TRUE)) |>
  ungroup()

full_dt |>
  pull(bt) |>
  summary()

full_dt |>
  mutate(bt = is.na(bt)) |>
  pull(bt) |>
  mean()

## 2015 useles
full_dt |>
  group_by(year) |>
  summarise(na_bt = mean(is.na(bt))) |>
  print(n = Inf)

full_dt <- full_dt |>
  filter(!is.na(bt))

##--- making it spatial ----

full_dt <- full_dt |>
  st_as_sf(coords = c("lon", "lat"),
           crs = st_crs(4326))

##--- bands ----

poly_border <- st_bbox(full_dt) |>
  st_as_sfc()

dat_gunion <- full_dt |>
  st_geometry() |>
  st_union()

dat_chull <- dat_gunion |>
  st_convex_hull()

##--- getting US map to refine polygons ----

states <-
  rnaturalearth::ne_states(
                     country = c("United States of America",
                                 "Canada"),
                     returnclass = "sf"
                 ) |>
  st_transform(4326)

bbox_buffer <- st_bbox(full_dt)
bbox_shrink <- bbox_buffer
## expanding bbox size
bbox_buffer[1:2] <- bbox_buffer[1:2] - 3
bbox_buffer[3:4] <- bbox_buffer[3:4] + 3
## shrinking bbox size
bbox_shrink[4] <- bbox_shrink[4] - .1

states <-
  states |>
  filter(st_intersects(geometry,
                       poly_border,
                       sparse = FALSE)[, 1]) |>
  st_geometry()

shore <-
  st_intersection(states,
                  st_set_crs(st_as_sfc(bbox_shrink),
                                       4326)) |>
  st_union() |>
  st_make_valid()

coastline <-
  shore |>
  st_cast("POLYGON") |>
  st_cast("LINESTRING") |>
  st_intersection(st_set_crs(poly_border,
                             4326)) |>
  st_union()

water <- 
  st_sym_difference(shore,
                    st_set_crs(st_as_sfc(bbox_shrink),
                               4326)) |>
  st_make_valid()

water <- water[2] |>
  st_cast("POLYGON")

water <- water[1] |>
  nngeo::st_remove_holes() |>
  st_make_valid()

dat_chull <-
  dat_chull |>
  st_intersection(water)

plot(water, col = scales::alpha(4, .2))
plot(dat_chull, col = scales::alpha(4, .2),
     border = 2, lwd = 2, add = TRUE)

## new study regions
study_region <-
  st_make_grid(x = dat_chull,
               n = c(1, 11))[dat_chull] |>
  st_make_valid() |>
  st_intersection(dat_chull)

study_region <- 
  study_region |>
  st_make_valid() |>
  st_as_sf()

st_geometry(study_region) <- "geometry"

study_region <- study_region |>
  mutate(area = st_area(geometry)) |>
  mutate(area = units::set_units(area, "km^2")) |>
  filter(as.numeric(area) > 0)

plot(water, col = scales::alpha(4, .2))
plot(st_geometry(study_region),
     col = scales::alpha(4, .2),
     border = 2, lwd = 2, add = TRUE)

plot(study_region)

##--- aggregating the data to bands ----

study_region <-
  study_region |>
  mutate(id = row_number(),
         .before = 1) |>
  relocate(area, .before = geometry)

band_dt <-
  study_region |>
  st_join(full_dt, st_contains) |>
  st_drop_geometry()

band_dt <-
  band_dt |>
  group_by(id, year) |>
  summarise(area = mean(area),
            region = first(region),
            num = sum(num),
            bt = mean(bt),
            depth = mean(depth)) |>
  ungroup()

plot(study_region["id"])
