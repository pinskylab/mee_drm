library(dplyr)
library(sf)
library(terra)
library(ggplot2)

temps <- rast("data/chelsa/tas_2019_avg.tif")

plot(temps)

temps2 <-
  temps * .1 - 273.15

temps3 <-
  temps2 * 1.8 + 32

plot(temps3)

## somewhat similar to:
## https://www.climate.gov/news-features/featured-images/new-maps-annual-average-temperature-and-precipitation-us-climate
