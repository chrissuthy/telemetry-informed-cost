library(ggplot2)
library(dplyr)
library(raster)
library(viridis)
library(patchwork)
library(ggnewscale)
library(purrr)
library(reshape2)
select = dplyr::select
extract = raster::extract
mutate = dplyr::mutate


#----Load data----

landscape <- readRDS("output/small ups/model_data/landscape.RData")[[1]]
tracks <- readRDS("output/small ups/model_data/tracks_all.RData") 
traps <- readRDS("output/small ups/model_data/traps.RData")
track <- tracks %>%
  pluck(1) %>%
  filter(id == 20)
prob_thin <- 0.1
landscape <- landscape^2
cost <- 1  
landscape <- exp(cost * landscape)
df.l <- as.data.frame(landscape, xy=T)

#----Pixels for traps----

# Get pixels for each trap, assign trap number
tt_p <- extract(x = landscape, y = traps, cellnumber=T)[,1]
traps_pxs <- traps %>% 
  mutate(pixel = tt_p) %>%
  mutate(trap = 1:n())


#----Pixels for fixes----

# Get pixels for each fix, assign
fx_p <- extract(x = landscape, y = track[,c("x","y")], cellnumber=T)[,1]
fixes_pxs <- track %>%
  mutate(pixel = fx_p) %>%
  as.data.frame()


#----Tallying captures by days (across 10 sims)----

tracks_w_traps <- fixes_pxs %>%
  # Join traps to fixes via pixels
  left_join(
    x = ., y = traps_pxs, 
    by = c("x", "y", "pixel")) %>% 
  # Remove fixes not at traps
  na.omit() %>% 
  # Convert times to K
  mutate(K = as.numeric(
    cut(x = times, 
        breaks = seq(from = 0, to = (90*24), by = 24), 
        right = F, include.lowest = T))) %>%
  # Select only relevant data
  select(id, trap, K) %>%
  # Group it
  group_by(id, trap, K) %>%
  # Summarize obs into n traps per K
  summarise(caps = n()) %>%
  # Keep only unique rows 
  distinct() %>%
  as.data.frame()

# Make a big df of all possible options (to fill in missing data)
tofill_df <- expand.grid(
  id = 20,      # individuals
  trap = 1:100, # 100 traps
  K = 1:90,     # 90 days
  caps = 0      # 0 caps (filler)
)

# From all options, only keep those not found in the real data
filler_df <- setdiff(
  tofill_df[,c("id", "trap", "K")], 
  tracks_w_traps[,c("id", "trap", "K")]) %>%
  as.data.frame()
filler_df$caps = 0

# Combine filler data so all options are included
cap_hist <- rbind(tracks_w_traps, filler_df)


#----Converting to 4D array----

set.seed(23)

# Sim 1:10, id 1:50, trap 1:100, K 1:90
y <- acast(cap_hist, id~trap~K, value.var = "caps", fill=0)
dim(y)
y[y>0] <- 1

# THINNING
to_thin <- y[y>0]
y[y>0] <- rbinom(length(to_thin), 1, prob_thin)

# Final scr data
scr <- apply(y, 1:2, sum)

traps$dets <- as.numeric(scr)
traps$pres <- 0
traps$pres[traps$dets > 0] <- 1
traps$pres <- factor(traps$pres, levels = c(0, 1))
traps$dets2 <- traps$dets
traps$dets2[traps$dets2 == 0] <- NA

# # Simulate detection data
# traps$pres <- rbinom(nrow(traps), 1, 0.1)
# traps$dets <- 0
# traps$dets[traps$pres > 0] <- rpois(sum(traps$pres > 0), 3)
# traps$pres <- factor(traps$pres, levels = c(0, 1))

#----Main figure----

# Landscape outline
df_rect <- data.frame(
  xmin = 0,
  ymin = 0,
  xmax = max(df.l$x)+0.1,
  ymax = max(df.l$y)+0.1
)

#traps$dets <- traps$dets+1

traps <- traps %>%
  mutate(y = if_else(y > 21 & y < 22, 21.625, y))

# Landscape
p1 <- ggplot() +
  # geom_rect(data = df_rect, 
  #           color = "black", size = 2,
  #           mapping = aes(
  #             xmin = xmin, xmax = xmax,
  #             ymin = ymin,  ymax = ymax)) +
  geom_tile(data = df.l, aes(x=x, y=y, fill=layer)) +
  scale_fill_viridis_c(NULL, direction = -1) +
  geom_rect(data = df_rect, fill=NA,
            color = "black", size = 1,
            mapping = aes(
              xmin = 11.9, xmax = 25.9,
              ymin = 8.3,  ymax = 22.3)) +
  geom_point(data=traps, pch = 3, color = "black",
             aes(x=x, y=y), size = 0.9, stroke = 0.8) +
  coord_equal() +
  theme_minimal() +
  labs(x=NULL, y=NULL, title = "(a) Generate cost surface") +
  theme(
    plot.title = element_text(hjust=0.5, face = "bold", size=16),
    axis.text = element_text(size=16),
    legend.position = "none",
    panel.grid = element_blank())

# Track
p2 <- ggplot() +
  lims(x = c(11.9, 25.9), 
       y = c(8.3, 22.3)) +
  geom_rect(data = df_rect, fill="yellow",
            color = "black", size = 2,
            mapping = aes(
              xmin = 11.9, xmax = 25.9,
              ymin = 8.3, ymax = 22.3)) +
  geom_tile(data = df.l, aes(x=x, y=y), fill = "white") +
  new_scale("fill") +
  geom_tile(data = df.l, aes(x=x, y=y, fill=layer)) +
  scale_fill_viridis_c("Covariate", alpha=0.55, direction = -1) +
  geom_path(data=track, aes(x=x, y=y), 
            color = alpha("white", 0), size=0.8) +
  geom_path(data=track, aes(x=x, y=y), 
            color = alpha(viridis(5)[2], 0.4), size=0.7) +
  geom_rect(data = df_rect, fill=NA,
            color = "black", size = 2,
            mapping = aes(
              xmin = 11.9, xmax = 25.9,
              ymin = 8.3, ymax = 22.3)) +
  coord_equal() +
  theme_minimal() +
  labs(x=NULL, y=NULL, title = "(b) Simulate movement data") +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(hjust=0.5, face = "bold", size=16),
    axis.text = element_blank(),
    legend.position = "none",
    panel.grid = element_blank()); p2

# Traps
p3 <- ggplot() +
  lims(x = c(11.9, 25.9), 
       y = c(8.3, 22.3)) +
  # White background
  geom_rect(data = df_rect, fill=NA,
            color = "black", size = 2,
            mapping = aes(
              xmin = 11.9, xmax = 25.9,
              ymin = 8.3,  ymax = 22.3)) +
  # Cost surface
  geom_tile(data = df.l, aes(x=x, y=y), fill = "white") +
  # Movement track
  geom_path(data=track, aes(x=x, y=y), 
            color = alpha(viridis(5)[2], 0.2), size=0.9) +
  # SCR traps
  new_scale("color") +
  # # Border
  # geom_point(data=traps, pch = 21, color = "black",
  #            aes(x=x, y=y, size = dets+0.75)) +
  # Fill
  geom_point(data=traps, pch = 3,
             aes(x=x, y=y), size = 1.6, stroke=0.9) +
  geom_point(data=traps, pch = 1,
             aes(x=x, y=y, color = pres, size = dets+3), stroke=3) +
  geom_point(data=traps, pch = 16,
             aes(x=x, y=y, color = pres, size = dets+3)) +
  # Color
  scale_color_manual(values=c(NA, "black")) +
  # Label SCR traps
  geom_text(data=traps, aes(x=x, y=y, label = dets2, size = dets2+5), 
            color = "white", fontface="bold") +
  # Formattting
  #coord_equal() +
  theme_minimal() +
  labs(x=NULL, y=NULL, title = "(c) Collect as SCR data") +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(hjust=0.5, face = "bold", size=16),
    axis.text = element_blank(),
    legend.position = "none",
    panel.grid = element_blank()); p3

# Together 1 x 3
ps_small <- p1 + p2 + p3
ps_small
