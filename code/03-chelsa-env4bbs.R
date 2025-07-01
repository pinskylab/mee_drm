library(dplyr)
library(sf)
library(terra)
library(ggplot2)

my_dt <- readRDS("data/birds/no_env.rds")

my_map <- st_read("data/birds/shape/grid.shp")

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

avg_maps <- lapply(avg_maps, \(x) {
  names(x)[2] <- "tavg"
  x
})

avg_maps <- bind_rows(avg_maps)

my_dt <- left_join(my_dt, avg_maps,
                   by = c("id" = "ID",
                          "year"))

## if the data is in tens of kelvin:
## my_dt$tavg / 10 - 273.15

my_dt <-
  my_dt |>
  filter(!is.na(tavg))

my_dt <-
  my_dt |>
  mutate(tavg = tavg / 10 - 273.15) |>
  mutate(fahrenreit = tavg * 1.8 + 32)

my_dt <- my_dt |>
  mutate(id = as.integer(factor(id))) |>
  arrange(id, year)

saveRDS(my_dt, "data/birds/processed.rds")
