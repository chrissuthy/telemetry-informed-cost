library(dplyr)
library(ggplot2)
library(tidyr)

#----Raw data----

# Get the tracks object
tracks <- readRDS("output/10sims_tracks_forTraps.RData")

# Make thihs into a full df
df <- tracks %>%
  do.call(rbind, .) %>%
  mutate(sim = rep(1:10, each = 90*24*50))
  
# Plot all of these data
ggplot() +
  geom_path(data = df, aes(x=x, y=y, color = as.factor(id), group = as.factor(id))) +
  geom_point(data = traps, aes(x=x, y=y), size = 0.2) +
  facet_wrap(~sim, nrow = 2) +
  theme_minimal() +
  theme(aspect.ratio = 1,
        legend.position = "none")


#----Pixels for traps----
  
# Get pixel values for each trap, assign trap number
traps 
trap_pxs <- extract(x = landscape_r, y = traps, cellnumber=T)[,1]
traps_pxs <- traps %>% 
  mutate(pixel = trap_pxs) %>%
  mutate(trap = 1:n())


#----Pixels for tracks----

# Summarize spatial locations of fixes
fixes_per_xy <- df %>% 
  filter(sim == 1) %>%
  group_by(id, x, y) %>%
  summarise(n = n()) %>%
  ungroup()

# Get pixels for thhose locations
fix_pxs <- extract(x = landscape_r, y = fixes_per_xy[,c("x","y")], cellnumber=T)[,1]

# Assign those pixels to te original data
fixes_pxs <- fixes_per_xy %>%
  mutate(pixel = fix_pxs)


#----Tallying captures----

# Join trap locs to track locs
# keep only fixes on traps
# tally the number of captures per trap
tracks_w_traps <- left_join(
  x = fixes_pxs, y = traps_pxs, 
  by = c("x", "y", "pixel")) %>% 
  na.omit() %>% 
  group_by(id, trap) %>% 
  summarise(caps = n()*n)

# Number of individuals captured
n_inds <- length(unique(tracks_w_traps$id))

# Number of captures per individual
n_caps_per_ind <- tracks_w_traps %>%
  group_by(id) %>%
  summarise(caps_total = sum(caps)) %>%
  pull(caps_total)
hist(n_caps_per_ind, breaks = 20)

# Number of individuals captured on 2 traps
n_traps_per_ind <- tracks_w_traps %>%
  group_by(id) %>%
  summarise(traps_total = n()) %>%
  pull(traps_total)
hist(n_traps_per_ind, breaks = 20)
as.numeric(table(n_traps_per_ind)[2])
