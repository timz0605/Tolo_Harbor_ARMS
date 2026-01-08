
library(tidyverse)
library(sf)
library(gstat)
library(raster)
library(stars)
library(viridis)
library(ggspatial)
library(geodata)

setwd("C:/Baker_Lab/ARMS_data/COI_TOLO_Harbour/analyses/spatial_interpolation/")

# 1. Download HK Map (Level 0 is the country outline)
# path = tempdir() saves it to a temporary folder so it doesn't clutter your PC
hk_map <- gadm(country = "HKG", 
               level = 0, 
               path = "C:/Baker_Lab/ARMS_data/COI_TOLO_Harbour/analyses/",
               version = "3.6",
               resolution = 1) %>% 
  st_as_sf()

sz_map <- gadm(country = "CHN", 
               level = 0, 
               path = "C:/Baker_Lab/ARMS_data/COI_TOLO_Harbour/analyses/",
               version = "3.6",
               resolution = 1) %>% 
  st_as_sf()

# 2. Define the "Box" for Tolo Harbor and Shenzhen
tolo_bbox <- st_bbox(c(xmin = 114.15, xmax = 114.50, 
                       ymin = 22.38, ymax = 22.70), 
                     crs = st_crs(hk_map))

sz_bbox <- st_bbox(c(xmin = 114.10, xmax = 114.50, 
                     ymin = 22.50, ymax = 22.65), 
                   crs = st_crs(sz_map))

# 3. Crop the HK and Shenzhen map to just the Tolo region
tolo_region <- st_crop(hk_map, tolo_bbox)
sz_region <- st_crop(sz_map, sz_bbox)

land_layer <- rbind(
  tolo_region %>% dplyr::select(geometry), 
  sz_region %>% dplyr::select(geometry)
)

# Check what we have so far
ggplot() + 
  geom_sf(data = tolo_region) + 
  ggtitle("Basemap of Tolo Harbor Region")


# 4. Load point data and do filtering and summarizing
point_data <- read_csv("C:/Baker_Lab/ARMS_data/COI_TOLO_Harbour/analyses/spatial_interpolation/mirs_tolo_combined.csv")


# --- Summarise stations ---

point_data_summary <- point_data %>% 
  dplyr::select(-`Sample No`, -`Dates`) %>%
  filter(Depth == "Bottom Water") %>% 
  group_by(Station) %>%
  summarise(
    Latitude  = first(Latitude),
    Longitude = first(Longitude),
    across(`5dOD`:`VSS`, mean, na.rm = TRUE),
    .groups = "drop"
  )

# write_csv(point_data_summary, "tolo_point_data_summary.csv")

names(point_data_summary)[names(point_data_summary) == "5dOD"] <- "OD"

station_sf <- st_as_sf(
  point_data_summary,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)

sites <- data.frame(
  x = c(114.221275, 114.290992, 114.35645, 114.438564),
  y = c(22.43725, 22.463206, 22.501042, 22.542639),
  label = c("Center Island", "Che Lei Pai", "Port Island", "Tung Ping Chau"),
  color = c("#df8d71", "#d8b847", "#75884b", "#5b859e")
)


# --- Grid ---
grid <- st_bbox(
  c(xmin=114.16, xmax=114.47, ymin=22.40, ymax=22.58),
  crs = st_crs(4326)
) %>%
  st_as_stars(dx = 0.0045, dy = 0.0045)
st_crs(grid) <- 4326

# --- IDW ---
interpolated <- gstat::idw(
  formula = Sal ~ 1, #change the parameter we are interpolating
  locations = station_sf,
  newdata = grid,
  idp = 2
)

# --- Proper 6 km buffer (projected) ---
station_sf_m <- st_transform(station_sf, 32650)
valid_buffer <- st_buffer(station_sf_m, dist = 6000) %>% st_union()
valid_buffer <- st_transform(valid_buffer, 4326)

# --- Convert raster → points ---
interpolated_sf <- st_as_sf(interpolated, as_points = TRUE, na.rm = FALSE)

# --- Clip ---
interpolated_clipped_sf <- st_filter(interpolated_sf, valid_buffer)

# --- Prepare for geom_tile ---
coords <- st_coordinates(interpolated_clipped_sf)

plot_data <- interpolated_clipped_sf %>%
  st_drop_geometry() %>%
  mutate(
    x = coords[, 1],
    y = coords[, 2]
  ) %>%
  filter(!is.na(var1.pred))

Sal <- ggplot() + #change the name of the new figure
  geom_tile(
    data = plot_data,
    aes(x = x, y = y, fill = var1.pred),
    color = "white", # border color, NA or white
    size = 0.1
  ) +
  geom_sf(data = land_layer, fill = "grey90", color = "white") +
  geom_sf(data = station_sf, size = 1.5, shape = 19) +
  geom_point(data = sites, 
             aes(x = x, y = y), 
             size = 1.5,
             shape = 19, 
             color = "red") +
  scale_fill_gradient2(
    low = "#2c7bb6", mid = "#ffffbf", high = "#d7191c",
    midpoint = mean(point_data_summary$Sal, na.rm = TRUE), # remember to changge the parameter here
    name = "Sal (PSU)", # change the legend title
    na.value = "transparent"
  ) +
  annotation_scale(
    location = "br",      # 'br' stands for bottom-right
    width_hint = 0.18,     # Adjusts how much of the plot width it occupies
    line_width = 1,
    text_cex = 0.8,
    style = "bar"       # Can be 'ticks' or 'bar'
  ) +
  coord_sf(
    xlim = c(114.16, 114.47),
    ylim = c(22.40, 22.58),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(color = "transparent"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = c(0.03, 0.97),
    legend.justification = c("left", "top"), 
    legend.background = element_rect(fill = alpha("white", 0.8)),
    legend.title = element_text(size = 10, face = "bold")
  )

Sal

ggsave("idw_Sal.png", plot = Sal, width = 8, height = 6, dpi = 600)


library(patchwork)

env_plot_patch <- OD + AmN + Chlo + DO + FC + OrthoP + Phaeo + Sal + 
  plot_layout(
    ncol = 2, 
    axes = "collect"
  ) +
  plot_annotation(tag_levels = "a")
env_plot_patch

ggsave("idw_env.png", 
       plot = env_plot_patch, 
       width = 12, 
       height = 16,
       dpi = 600)


################################
# Regression and tests to show gradients 
################################

gradient_site <- c("TM2", "TM4", "TM6", "TM7", "TM8", "MM17", "MM5")

gradient_station <- ggplot() + #change the name of the new figure
  geom_sf(data = land_layer, fill = "grey90", color = "white") +
  geom_sf(
    data = station_sf %>% filter(Station %in% gradient_site), 
    size = 2, 
    shape = 19
  ) +
  geom_sf_text(
    data = station_sf %>% filter(Station %in% gradient_site),
    aes(label = Station),
    size = 3,
    nudge_y = -0.007,      # Adjusts vertical position slightly
    fontface = "bold"
  ) +
  geom_point(data = sites, 
             aes(x = x, y = y), 
             size = 2,
             shape = 19, 
             color = "red") +
  annotation_scale(
    location = "br",      # 'br' stands for bottom-right
    width_hint = 0.18,     # Adjusts how much of the plot width it occupies
    line_width = 1,
    text_cex = 0.8,
    style = "bar"       # Can be 'ticks' or 'bar'
  ) +
  coord_sf(
    xlim = c(114.16, 114.47),
    ylim = c(22.40, 22.58),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    axis.title = element_blank()
  )
  
gradient_station  

ggsave("gradient_stations.png", 
       plot = gradient_station, 
       width = 7, 
       height = 5,
       dpi = 600)

###################
# Regression tests
###################

library(ggpmisc)
library(scales)
library(MetBrewer)

met_pal <- met.brewer("Juarez", n = 7)
site_color <- setNames(met_pal, gradient_site)


gradient_data <- point_data_summary %>%
  filter(Station %in% gradient_site) %>%
  arrange(match(Station, gradient_site)) %>%
  dplyr::select(Station, OD, AmN, Chlo, DO, FC, OrthoP, Phaeo, Sal) %>%
  mutate(distance_km = c(0, 2.27, 4.86, 7.57, 11.89, 17.95, 22.92)) %>% 
  mutate(`log(FC)` = log(FC + 1)) %>% # log transform FC to reduce skewness
  dplyr::select(-FC) # remove original FC


# 1. Reshape the data so all parameters are in one column
long_data <- gradient_data %>%
  pivot_longer(
    cols = -c(Station, distance_km), # Keep these, pivot everything else
    names_to = "parameter", 
    values_to = "value"
  )

# 2. Create the plot
gradient_regression <- ggplot(long_data, aes(x = distance_km, y = value)) +
  geom_smooth(
    method = "lm",
    se = TRUE, 
    linewidth = 1,
    color = "darkgrey",
    fill = "lightgrey",
    alpha = 0.2
    ) +
    geom_point(aes(color = Station), size = 3, alpha = 0.8, ) +
  scale_color_manual(values = site_color) + 
  labs(
    x = "Distance (km) from Inner Harbor Monitoring Station TM2",
    y = "Value",
    color = "Station"
    ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_blank()
        ) +
  scale_y_continuous(n.breaks = 5) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
  stat_poly_eq(
    aes(label = paste(..eq.label.., 
                      ..p.value.label.., 
                      ..rr.label..,
                      sep = "~")),
    parse      = TRUE,
    size        = 3,
    label.x = "right",
    label.y = "top",
    color = "black"
  )

gradient_regression

ggsave("gradient_regression.jpg", 
       plot = gradient_regression, 
       width = 9, 
       height = 12,
       dpi = 600)



