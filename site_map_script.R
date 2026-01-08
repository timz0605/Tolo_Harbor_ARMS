

library(ggmap)
library(tidyverse)
library(tmap)
library(shadowtext)
library(maps)
library(patchwork)
library(ggspatial)

## site_col <- c("#ca0e12", "#f6bd21", "#2aa7de", "#25377f")

site_col <- c("#df8d71", "#d8b847", "#75884b", "#5b859e")

# Set the bounding box for the area of Hong Kong you want to map

api_secret <- "AIzaSyAeje3nHm6FDhbObD2xhSFuQmcvMd2o294"
register_google(key = api_secret)

# Retrieve the map tile

setwd("C:/Baker_Lab/ARMS_data/COI_Tolo_Harbour/analyses/site_map")


### Approach 2

point_data <- data.frame(
  x = c(114.221275, 114.290992, 114.35645, 114.438564),
  y = c(22.43725, 22.463206, 22.501042, 22.542639),
  label = c("Center Island", "Che Lei Pai", "Port Island", "Tung Ping Chau"),
  color = c("#df8d71", "#d8b847", "#75884b", "#5b859e")
  )

get_googlemap(center = c(lon = 114.28, lat = 22.42), 
              zoom = 11, 
              maptype = "satellite", 
              size = c(640, 640), 
              language = "en-EN") %>% 
  ggmap() +
  geom_point(data = point_data, 
             aes(x = x, y = y, color = label), 
             size = 3, 
             shape = 21, 
             fill = site_col) +
  scale_color_manual(values = point_data$color, 
                     name = "Location") +
  guides(color = guide_legend(override.aes = list(shape = 21, size = 3, stroke = 1, fill = site_col)))


### Approach 3

annotation_data <- data.frame(
  x = c(114.43),
  y = c(22.50),
  label = c("Mirs Bay")
)

site_map <- get_googlemap(center = c(lon = 114.34, lat = 22.5), 
              zoom = 11, 
              maptype = "satellite", 
              size = c(640, 640), 
              language = "en-EN") %>% 
  ggmap() +
  geom_point(data = point_data, 
             aes(x = x, y = y, color = label, fill = label), 
             size = 5,
             shape = 21, 
             stroke = 1.5, 
             color = "white") +
  geom_shadowtext(data = annotation_data, aes(x = x, y = y, label = label), size = 7) +
  scale_color_manual(values = point_data$color, 
                     name = "Location") +
  scale_fill_manual(values = point_data$color, 
                    name = "Location") +
  guides(color = guide_legend(override.aes = list(shape = 21, size = 4, stroke = 1, fill = point_data$color)),
         fill = guide_legend(override.aes = list(shape = 21, size = 4, stroke = 1))) +
  theme(text = element_text(size = 16),
        axis.title = element_blank(),
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 18, face = "bold"), 
        legend.position = c(0.97, 0.03),
        legend.justification = c("right", "bottom"), 
        panel.background = element_rect(fill = "lightgray"), 
        legend.background = element_rect(fill = alpha("white", 0.7)))

site_map

ggsave("ARMS_site_map_updated.png", plot = site_map, width = 8.7, height = 8, dpi = 300)


### Map for Hong Kong alone

annotation_data_1 <- data.frame(
  x = c(114.19, 114.03),
  y = c(22.3, 22.6),
  label = c("Hong Kong", "Shenzhen")
)

HKSZ <- get_googlemap(center = c(lon = 114.08, lat = 22.5), 
                      zoom = 10, 
                      maptype = "satellite", 
                      size = c(640, 640), 
                      language = "en-EN") %>% 
  ggmap() +
  geom_shadowtext(data = annotation_data_1, aes(x = x, y = y, label = label), size = 7) +
  geom_rect(aes(xmin = 114.15, xmax = 114.5, ymin = 22.38, ymax = 22.64), 
           alpha = 0, color = "red", linewidth = 2) +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude") +
  theme(text = element_text(size = 16),
        axis.title = element_text(size = 18, face = "bold"),
  )
HKSZ


HKSZ <- HKSZ + ggtitle("(a)")
site_map <- site_map + ggtitle("(b)")

site_map_composite <- HKSZ + site_map
site_map_composite

ggsave("ARMS_site_map_composite.jpeg", plot = site_map_composite, width = 16, height = 8, dpi = 300)

