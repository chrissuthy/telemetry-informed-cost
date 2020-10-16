library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)

#----Raw data----

# Get the tracks object
tracks <- readRDS("output/oct14_N50_alpha2of2_10sims.RData")

# Make thihs into a full df
df <- tracks %>%
  do.call(rbind, .) %>%
  mutate(sim = rep(1:10, each = 90*24*50))


#----Pixels for traps----
  
# Get pixels for each trap, assign trap number
tt_p <- extract(x = landscape_r, y = traps, cellnumber=T)[,1]
traps_pxs <- traps %>% 
  mutate(pixel = tt_p) %>%
  mutate(trap = 1:n())


#----Pixels for fixes----

# Get pixels for each fix, assign
fx_p <- extract(x = landscape_r, y = df[,c("x","y")], cellnumber=T)[,1]
fixes_pxs <- df %>%
  mutate(pixel = fx_p)


#----Tallying captures by days----

# NEED TO FILL IN TIMES 1:(90*24) OR K 1:90
tracks_w_traps <- fixes_pxs %>%
  left_join(
    x = ., y = traps_pxs, 
    by = c("x", "y", "pixel")) %>% 
  filter(sim == 1) %>%
  na.omit() %>% 
  mutate(K = as.numeric(
    cut(x = times, 
        breaks = seq(1, 90*24, by = 24), 
        right = F))) %>%
  select(id, trap, K) %>%
  group_by(id, trap, K) %>%
  summarise(caps = n()) %>%
  distinct()

# Convert this into a 3D array
y <- acast(tracks_w_traps, id~trap~K, value.var = "caps", fill=0)
dim(y)



