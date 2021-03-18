# # # # # # # # # # # # # # #
#                           #
#  ____HEY!!_DO THIS!!____  #
#     Parameterize sims     #
#                           #
# # # # # # # # # # # # # # #

# Right here:
select_ups  <- c("small ups", "big ups")[1]

# # # # # # # # # # # # # # #
#                           #
#       GREAT, THANKS!      #
#                           #
# # # # # # # # # # # # # # #


library(dplyr)
library(raster)
library(ggplot2)
library(tidyr)
library(reshape2)
library(NLMR)
library(sf)
select = dplyr::select
extract = raster::extract
mutate = dplyr::mutate

#### THINNING ####
prob_thin <- 1


#----Load Data---

# Landscape
file <- paste0("output/", select_ups, "/model_data/landscape.RData")
landscape_r <- readRDS(file)[[1]] # colnames?

# State-space
file <- paste0("output/", select_ups, "/model_data/ss.RData")
ss <- readRDS(file)

# tracks
file <- paste0("output/", select_ups, "/model_data/tracks_all.RData")
tracks_all <- readRDS(file)

#----Make the traps----

# Movement model
sigma <- 1  #4.3

# SCR
N <- 100

# Trap array (space)
trap_array <- ss %>%
  filter(x >= (min(x)+3*sigma) & x < (max(x)-3*sigma)) %>%
  filter(y >= (min(y)+3*sigma) & y < (max(y)-3*sigma))

# Make traps
trap_array_sf <- st_as_sf(trap_array, coords = c('y', 'x'))
traps <- st_make_grid(trap_array_sf, n = c(15,15)) %>% # alternatively could use cellsize
  as_Spatial() %>%
  coordinates() %>%
  as.data.frame() %>%
  select(x = V2, y = V1) %>%
  extract(landscape_r, ., cellnumbers=T, df=T) %>%
  pull(cells) %>%
  xyFromCell(landscape_r, .) %>%
  as.data.frame()
rownames(traps) <- NULL
table(diff(traps$y))

nrow(traps)
plot(landscape_r)
lines(extent(ss))
points(traps)
lines(extent(traps))

#----Raw movement data----

# Get the tracks object
tracks <- tracks_all

# Make thihs into a full df
df <- tracks %>%
  do.call(rbind, .) %>%
  mutate(sim = rep(1:100, each = 90*24*N))

# # Plot
# par(mfrow = c(2,5), pty = "s")
# for(i in 1:10){
#   plot(landscape_r)
#   lines(y~x, data = df[df$sim==i,], col = alpha("red", 0.3))
#   points(traps, pch = 20) 
# }
# par(mfrow=c(1,1))

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
  select(id, trap, K, sim) %>%
  # Group it
  group_by(id, trap, K, sim) %>%
  # Summarize obs into n traps per K
  summarise(caps = n()) %>%
  # Keep only unique rows 
  distinct() %>%
  as.data.frame()

# Make a big df of all possible options (to fill in missing data)
tofill_df <- expand.grid(
  id = 1:N,     # 100 individuals
  trap = 1:nrow(traps), # 100 traps
  K = 1:90,     # 90 days
  sim = 1,      # There will always be 1:10 sims anyway
  caps = 0      # 0 caps (filler)
)

# From all options, only keep those not found in the real data
filler_df <- setdiff(
  tofill_df[,c("id", "trap", "K", "sim")], 
  tracks_w_traps[,c("id", "trap", "K", "sim")]) %>%
  as.data.frame()
filler_df$caps = 0

# Combine filler data so all options are included
cap_hist <- rbind(tracks_w_traps, filler_df)


#----Converting to 4D array----

# Sim 1:10, id 1:100, trap 1:100, K 1:90
y <- acast(cap_hist, sim~id~trap~K, value.var = "caps", fill=0)
dim(y)
y[y>0] <- 1

# THINNING
to_thin <- y[y>0]
y[y>0] <- rbinom(length(to_thin), 1, prob_thin)


#---Summary stats----

n_inds      <- c()
mean_caps   <- c()
caps_total  <- c()
caps_2traps <- c()
caps_3traps <- c()
caps_4traps <- c()
for(i in 1:100){
  
  tmp_y0 <- y[i,,,]
  
  ncap <- apply(tmp_y0,1, sum) # sum of captures for each individual
  tmp_y <- tmp_y0[ncap>0,,] # reduce the y array to include only captured individuals
  
  # Number of individuals captured
  n_inds[i] <- nrow(tmp_y)
  
  # Total captures
  caps_total[i] <- sum(tmp_y0)
  
  # Average number of captures per individual
  mean_caps[i] <- caps_total[i]/n_inds[i]
  
  # Number of individuals captured on 2 traps
  tmp_y_collapsed <- apply(tmp_y, 1:2, sum)
  tmp_y_collapsed[tmp_y_collapsed>0] <- 1
  caps_2traps[i] <- sum(apply(tmp_y_collapsed, 1, sum) == 2)
  
  # Number of individuals captured on 3 traps
  caps_3traps[i] <- sum(apply(tmp_y_collapsed, 1, sum) == 3)
  
  # Number of individuals captured on 4 traps
  caps_4traps[i] <- sum(apply(tmp_y_collapsed, 1, sum) == 4)
  
}


#---Plots----

par(mfrow=c(2,3))

# Number of individuals captured
hist(n_inds, breaks = 20, main = "n inds. captured", xlab = NA, col = "gray80")

# Average number of captures per individual
hist(mean_caps, breaks = 20, main = "Avg. captures per ind.", xlab = NA, col = "gray80")

# Total captures
hist(caps_total, breaks = 20, main = "Total captures", xlab = NA, col = "gray80")

# 2 traps
hist(caps_2traps, breaks = 20, main = "Inds. captured on 2 traps", xlab = NA, col = "gray80")

# 3 traps
hist(caps_3traps, breaks = 20, main = "Inds. captured on 3 traps", xlab = NA, col = "gray80")

# 3 traps
hist(caps_4traps, breaks = 20, main = "Inds. captured on 4 traps", xlab = NA, col = "gray80")

par(mfrow=c(1,1))

#----Save SCR data----

# Collapse to captured
y_ALL <- list()
for(i in 1:100){
  
  # Reduce to captured
  tmp_y <- y[i,,,]
  captured <- apply(tmp_y,1,sum)>0
  y_ALL[[i]] <- tmp_y[captured,,]
  
}

file <- paste0("output/", select_ups, "/model_data/y_NoThin_moreTraps.RData")
saveRDS(y_ALL, file)

file <- paste0("output/", select_ups, "/model_data/traps_moreTraps.RData")
saveRDS(traps, file)
