library(dplyr)
library(sf)
library(terra)
library(arrow)
library(geoarrow)
library(ggplot2)

my_dt <- open_dataset("data/birds/no_env.parquet") |>
  st_as_sf()

my_map <- my_dt |>
  filter(year == 1980) |>
  select(id)

my_vec <- vect(my_map)

##--- process env data ----

avg_files <-
  list.files("data/chelsa",
             pattern = "avg",
             full.names = TRUE)

avg_maps <-
  lapply(avg_files, \(x) {
    year <- gsub("\\D", "", x) |>
      as.integer()
    out <- rast(x)
    extract(out, my_vec, mean,
            na.rm = TRUE) |>
      mutate(year = year)
  })

avg_maps <- bind_rows(avg_maps)

avg_maps <-
  avg_maps |>
  tidyr::pivot_longer(c(2, 4:NCOL(avg_maps)),
                      values_to = "tavg") |>
  select(- name) |>
  filter(!is.na(tavg))

my_dt <- left_join(my_dt, avg_maps,
                   by = c("id" = "ID",
                          "year"))

my_dt <-
  my_dt |>
  filter(!is.na(tavg))
