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
df.l <- as.data.frame(landscape, xy=T)
tracks <- readRDS("output/small ups/model_data/tracks_all.RData") 
traps <- readRDS("output/small ups/model_data/traps.RData")
track <- tracks %>%
  pluck(1) %>%
  filter(id %in% c(5, 20, 30, 15)) %>%
  filter(times < (90*24)*0.25)

# Track
p1 <- ggplot(data = df.l, aes(x=x, y=y)) +
  geom_tile(fill = "gray70") +
  new_scale("fill") +
  geom_tile(aes(fill=layer)) +
  scale_fill_viridis_c(alpha = 0.5) +
  geom_path(data=track, aes(x=x, y=y, group = id), 
            color = alpha("black", 0.65),
            size=0.35) +
  coord_equal() +
  theme_minimal() +
  labs(x=NULL, y=NULL, title = "High-res") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text = element_text(size=16),
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

# Track
p2 <- ggplot(data = df.l, aes(x=x, y=y)) +
  geom_tile(fill = "gray70") +
  new_scale("fill") +
  geom_tile(aes(fill=layer)) +
  scale_fill_viridis_c(alpha = 0.5) +
  geom_path(data=track, aes(x=x, y=y, group = id), 
            color = alpha("black", 0.65),
            size=0.35) +
  coord_equal() +
  theme_minimal() +
  labs(x=NULL, y=NULL, title = "Low-res") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text = element_text(size=16),
    legend.position = "none",
    panel.grid = element_blank())


#----FINAL plot----
p1 + p2
