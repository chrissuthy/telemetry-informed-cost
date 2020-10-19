library(dplyr)
library(raster)
library(ggplot2)
library(tidyr)
library(reshape2)
library(NLMR)
select = dplyr::select
extract = raster::extract
mutate = dplyr::mutate

#----Make the traps----

# Settings: Parameters
upsilon <- 250
sigma <- 1000

# Settings: Landscape
ncol <- nrow <- 131
rr <- 250
hr95_lim <- (3*sigma) + (3*upsilon)

# Landscape
landscape0 <- nlm_gaussianfield(
  ncol = ncol, nrow = nrow, 
  resolution = rr, autocorr_range = 7)
landscape_r <- landscape0^2 # I STILL DON'T LIKE THIS! RESTRICTS LOW COST
landscape <- as.data.frame(landscape_r, xy=T)

# State-space
ss <- landscape %>%
  filter(x >= (min(x)+hr95_lim) & x < (max(x)-hr95_lim)) %>%
  filter(y >= (min(y)+hr95_lim) & y < (max(y)-hr95_lim))
nrow(ss)

# Trap array (space)
trap_array <- ss %>%
  filter(x >= (min(x)+3*sigma) & x < (max(x)-3*sigma)) %>%
  filter(y >= (min(y)+3*sigma) & y < (max(y)-3*sigma))

# Traps
t_coords_all <- unique(trap_array$x)
t_coords_label <- c(0, rep(c(1, rep(0, 7)),floor(length(t_coords_all)/8)),1,0)
t_coords_select <- data.frame(
  coord = t_coords_all, 
  indx = t_coords_label) %>%
  filter(indx == 1) %>%
  pull(coord)
traps <- expand.grid(x = t_coords_select, y = t_coords_select)


#----Raw movement data----

# Get the tracks object
tracks <- readRDS("output/oct18_N50_alpha2of2_psi05_10sims.RData")

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
  distinct()

# Make a big df of all possible options (to fill in missing data)
tofill_df <- expand.grid(
  id = 1:50,    # 50 individuals
  trap = 1:100, # 100 traps
  K = 1:90,     # 90 days
  sim = 1,      # There will always be 1:10 sims anyway
  caps = 0      # 0 caps (filler)
)

# From all options, only keep those not found in the real data
filler_df <- setdiff(
  tofill_df[,c("id", "trap", "K", "sim")], 
  tracks_w_traps[,c("id", "trap", "K", "sim")])
filler_df$caps = 0

# Combine filler data so all options are included
cap_hist <- rbind(tracks_w_traps, filler_df)

#----Converting to 4D array----

# Sim 1:10, id 1:50, trap 1:100, K 1:90
y <- acast(cap_hist, sim~id~trap~K, value.var = "caps", fill=0)
dim(y)
y[y>0] <- 1


#---Summary stats----

n_inds      <- c()
mean_caps   <- c()
caps_total  <- c()
caps_2traps <- c()
caps_3traps <- c()
caps_4traps <- c()
for(i in 1:10){
  
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
hist(n_inds, breaks = 20, main = "n inds. captured", xlab = NA)

# Average number of captures per individual
hist(mean_caps, breaks = 20, main = "Avg. captures per ind.", xlab = NA)

# Total captures
hist(caps_total, breaks = 20, main = "Total captures", xlab = NA)

# 2 traps
hist(caps_2traps, breaks = 20, main = "Inds. captured on 2 traps", xlab = NA)

# 3 traps
hist(caps_3traps, breaks = 20, main = "Inds. captured on 3 traps", xlab = NA)

# 3 traps
hist(caps_4traps, breaks = 20, main = "Inds. captured on 4 traps", xlab = NA)

par(mfrow=c(1,1))


