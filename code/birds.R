library(dplyr)
library(sf)
library(ggplot2)

weather <- readr::read_csv("data/birds/weather.csv")
vehic   <- readr::read_csv("data/birds/VehicleData/VehicleData.csv")
fifty   <- readr::read_csv("data/birds/50-StopData/1997ToPresent_SurveyWide/fifty1.csv")

