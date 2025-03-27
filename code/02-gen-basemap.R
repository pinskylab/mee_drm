library(sf)
library(dplyr)
library(arrow)
library(geoarrow)
library(ggplot2)
library(patchwork)

##--- load satellite data ----

satellite_sst <-
  open_dataset("data/copernicus/processed.parquet") |>
  st_as_sf()

## Cordinate Reference System
base_crs <- st_crs(satellite_sst)

##--- study region bounding box ----

region <- c("xmin" = -75.8, "ymin" = 35.16,
            "xmax" = -65.64, "ymax" = 44.43) |>
  st_bbox() |>
  st_as_sfc() |>
  st_set_crs(base_crs)

##--- generating 100 patches ----

patches <-
  st_make_grid(region, n = c(10, 10),
               what = "polygons")

##--- filtering patches that intersect with satellite data ----

satellite_region <-
  satellite_sst |>
  select(geometry) |>
  distinct() |>
  st_geometry() |>
  st_union() |>
  st_convex_hull()

##--- base map ---

base_map <-
  rnaturalearth::ne_states(
                     country = c("United States of America",
                                 "Canada"),
                     returnclass = "sf"
                 ) |>
  st_transform(base_crs)

base_map <-
  base_map |>
  select(postal) |>
  st_crop(region)

##--- averaging satellite data ----

patches <-
  patches |>
  st_as_sf() |>
  mutate(id = row_number(),
         .before = x)

st_geometry(patches) <- "geometry"

patches_map <- patches

patches <- patches |>
  st_join(satellite_sst,
          st_contains) |>
  rename("sal" = "so",
         "sst" = "thetao") |>
  mutate(across(sal:vo, \(x) units::set_units(x, NULL))) |>
  st_drop_geometry() |>
  tidyr::pivot_longer(cols = sal:vo,
                      names_to = "qty",
                      values_to = "obs")

patches <-
  patches |>
  group_by(id, time, qty) |>
  summarise(min = min(obs, na.rm = FALSE),
            max = max(obs, na.rm = FALSE),
            avg = mean(obs, na.rm = FALSE),
            med = median(obs, na.rm = FALSE)) |>
  ungroup() |>
  tidyr::pivot_wider(names_from = qty,
                     values_from = min:med,
                     names_vary = "slowest")

patches <- patches |>
  filter(!is.na(min_sal)) |>
  mutate(across(min_sal:med_vo,
                .fns = \(x) as.numeric(scale(x)),
                .names = "cs_{.col}"))

patches <- patches |>
  mutate(year = lubridate::year(time),
         .before = time)

patches <- patches_map |>
  left_join(patches, by = "id") |>
  filter(!is.na(min_sal))

##--- producing maps ----

range_sst <-
  satellite_sst |>
  filter(time == lubridate::as_date("1993-04-01")) |>
  mutate(sst = units::set_units(thetao, NULL)) |>
  pull(sst) |>
  range()

plot_a <-
  satellite_sst |>
  filter(time == lubridate::as_date("1993-04-01")) |>
  mutate(sst = units::set_units(thetao, NULL)) |>
  ggplot(data = _) +
  geom_sf(aes(color = sst),
          pch = 15) +
  geom_sf(data = base_map) +
  geom_sf_text(data = base_map,
               aes(label = postal),
               size = 2) +
  scale_color_viridis_c(option = "H") +
  labs(x = NULL, y = NULL,
       title = "(a)",
       color = "SST") +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.01, 0.8),
        legend.key.size = unit(.3, units = "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 5),
        legend.text.position = "left",
        legend.background = element_rect(color = 1,
                                         linewidth = .1))


plot_b <-
  patches |>
  filter(year == 1993) |>
  mutate(sst = med_sst) |>
  ggplot(data = _) +
  geom_sf(aes(fill = sst),
          pch = 15) +
  geom_sf(data = base_map) +
  geom_sf_text(data = base_map,
               aes(label = postal),
               size = 2) +
  scale_fill_viridis_c(option = "H",
                       limits = range_sst) +
  labs(x = NULL, y = NULL,
       title = "(b)") +
  guides(fill = "none") +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.title = element_blank())

plot_a + plot_b

ggsave(filename = "overleaf/img/sst_sim.pdf",
       width = 6, height = 3.08,
       dpi = 300)

##--- save data ----

write_parquet(patches,
              "data/basemap.parquet")
