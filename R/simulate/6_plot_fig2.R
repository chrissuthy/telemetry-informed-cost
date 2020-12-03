library(ggplot2)
library(dplyr)
library(raster)
library(viridis)
library(patchwork)
library(ggnewscale)
library(purrr)
library(reshape2)
library(RColorBrewer)
select = dplyr::select
extract = raster::extract
mutate = dplyr::mutate


#----SMALL ups----

landscape <- readRDS("output/small ups/model_data/landscape.RData")[[1]]
tracks <- readRDS("output/small ups/model_data/tracks_all.RData") 
traps <- readRDS("output/small ups/model_data/traps.RData")
track <- tracks %>%
  pluck(1) %>%
  filter(id %in% c(5, 20, 30, 15)) %>%
  filter(times < (90*24)*0.25)
landscape <- landscape^2
cost <- 1  
landscape <- exp(cost * landscape)
df.l <- as.data.frame(landscape, xy=T)
df_rect <- data.frame(
  xmin = 0,
  ymin = 0,
  xmax = max(df.l$x)+0.1,
  ymax = max(df.l$y)+0.1
)

# Track
p1 <- ggplot() +
  geom_rect(data = df_rect, 
            color = "black", size = 2,
            mapping = aes(
              xmin = xmin, xmax = xmax,
              ymin = ymin,  ymax = ymax)) +
  geom_tile(data = df.l, aes(x=x, y=y), fill = "white") +
  new_scale("fill") +
  geom_tile(data = df.l, aes(x=x, y=y, fill=layer)) +
  scale_fill_viridis_c(alpha = 0.8, direction=-1) +
  geom_path(data=track, aes(x=x, y=y, group = id), 
            color = alpha("white", 0.65),
            size=0.60) +
  geom_path(data=track, aes(x=x, y=y, group = id),  
            color = alpha("red", 0.65),
            size=0.35) +
  coord_equal() +
  theme_minimal() +
  labs(x=NULL, y=NULL, title = "High-resolution") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size=12),
    legend.position = "none",
    panel.grid = element_blank())

#----BIG ups----

landscape <- readRDS("output/big ups/model_data/landscape.RData")[[1]]
df.l <- as.data.frame(landscape, xy=T)
tracks <- readRDS("output/big ups/model_data/tracks_all.RData") 
traps <- readRDS("output/big ups/model_data/traps.RData")
track <- tracks %>%
  pluck(1) %>%
  filter(id %in% c(5, 20, 30, 15)) %>%
  filter(times < (90*24)*0.25)
landscape <- landscape^2
cost <- 1  
landscape <- exp(cost * landscape)
df.l <- as.data.frame(landscape, xy=T)

# Track
p2 <- ggplot() +
  geom_rect(data = df_rect, 
            color = "black", size = 2,
            mapping = aes(
              xmin = xmin, xmax = xmax,
              ymin = ymin,  ymax = ymax)) +
  geom_tile(data = df.l, aes(x=x, y=y), fill = "white") +
  new_scale("fill") +
  geom_tile(data = df.l, aes(x=x, y=y, fill=layer)) +
  scale_fill_viridis_c(alpha = 0.8, direction=-1) +
  geom_path(data=track, aes(x=x, y=y, group = id), 
            color = alpha("white", 0.65),
            size=0.60) +
  geom_path(data=track, aes(x=x, y=y, group = id), 
            color = alpha("red", 0.65),
            size=0.35) +
  coord_equal() +
  theme_minimal() +
  labs(x=NULL, y=NULL, title = "Low-resolution") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size=12),
    legend.position = "none",
    panel.grid = element_blank())


#----FINAL plot----
p1 + p2
