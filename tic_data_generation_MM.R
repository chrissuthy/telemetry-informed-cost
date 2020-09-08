# Gates Dupont
# Adapted from code by Chris Sutherland
# August - September 2020

library(oSCR)
library(gdistance)
library(NLMR)
library(landscapetools)
library(raster)
library(sf)
library(trajr)
library(circular)
library(ggplot2)
library(viridis)
library(dplyr)
library(patchwork)
area <- raster::area


set.seed(202006)


"LANDSCAPE"

#----Create the landscape----

# Set autocorrelation range
autocorr <- 3

# Simulate landscape
gauss_nlm <- nlm_gaussianfield(ncol = 100, nrow = 100, resolution = 1,
                               autocorr_range = autocorr)

# Create state-space raster from landscape (axes are different)
use_this <- gauss_nlm
landscape0 <- rasterFromXYZ(cbind(coordinates(use_this)/10,values(use_this)))
landscape <- landscape0^2
statespace <- as.data.frame(landscape, xy=T)[,c("x", "y")] 
statespace <- as.matrix(statespace)


"SCR"

#----SCR parameters----

# Abundance to density
abundance <- 300
area <- prod(round(apply(bbox(statespace),1,diff),2))

abs.density <- abundance / area
pix.density <- abundance / nrow(statespace)

# Activity centers
idx <- sample(x = 1:nrow(statespace), size = abundance)
ac <- statespace[idx,]
rownames(ac) <- NULL

# Detection
p0 <- 0.1
K <- 5


"MOVEMENT MODEL PARAMETERS"

#----Movement model parameters----

# MM parameters
sigma <- 0.35   # space use - step length
psi <- 0.7         # This is pr(moved)
theta <- 1.1       # Space use - home range
alpha2 <- 1        # Cost value
n_fixes <- 90*24   # Number of pings
moved <- rbinom(n_fixes, 1, psi) # Vector of decisions to move


"SCR traps"

#----Create traps----

# Traps
traplocs <- as.matrix(
  expand.grid(X = seq(5-3*2*sigma,5+3*2*sigma,length=7),
              Y = seq(5-3*2*sigma,5+3*2*sigma,length=7)))

# Number of traps
n_traps <- nrow(traplocs)


"COST SURFACE"

#----Cost surface----

# Cost surface
logcost <- alpha2*landscape
cost <- exp(logcost)

# Create tansition surface
tr1 <- transition(
  cost,
  transitionFunction=function(x) 1/mean(x),
  directions=16)

# Corrected transition surface
tr1Corr <-geoCorrection(
  tr1,
  multpl=FALSE,
  scl=FALSE)


"DATA-GENERATION (telemetry)"

#---Start tracking individuals----

# Tag some of the individuals
telemetry_n <- 16
which_telemetered_acs <- sample(1:nrow(ac), telemetry_n)
telemetered_acs <- ac[which_telemetered_acs,]

tracks = list()
for(j in 1:telemetry_n){
  
  # Get starting position for the individual
  sbar <- telemetered_acs[j,]
  
  # Create vector of steps
  s.grid <- rep(NA, n_fixes)
  s.grid[1] <- extract(stack(landscape),matrix(sbar,nrow=1), cellnumber=T)[1]
  
  # Set progress bar
  print(paste("Individual", j, "of", telemetry_n))
  pb = txtProgressBar(min = 2, max = n_fixes, initial = 1, style=3) 
  
  # Loop to generate fix locations
  for(i in 2:n_fixes){
    
    setTxtProgressBar(pb,i)
    
    # If the animal doesn't move, assign same loc and skip the rest
    if(moved[i]==0){
      s.grid[i] <- s.grid[i-1]
      next # Next mean skip the rest of the for loop
    }
    
    # If the animal moved...
    
    # Calculate ecological distance from the last position to each pixel (for sigma)
    D <- costDistance(tr1Corr,  statespace[s.grid[i-1],], statespace)
    
    # Calculate euclidean distance from the activity center to each pixel (for theta)
    Dac <- sqrt( (statespace[s.grid[1],1] - statespace[,1] )^2  + (statespace[s.grid[1],2] - statespace[,2])^2  )
    
    # Create a movement kernel using ecoD & sigma and eucDac & theta
    kern <- exp((-D*D/(2*sigma*sigma))  - (Dac*Dac/(2*theta*theta)) )
    
    # Assign the current position a probability of zero (can't stay in the same spot, b/c Psi=1)
    kern[s.grid[i-1]] <- 0
    
    # Normalize the kernel cell values into probabilities
    kern <- kern/rowSums(kern)
    
    # Sample a pixel for the next step using kernel
    s.grid[i]<- sample( 1:length(kern), 1, prob=kern)
  }
  
  
  #----Organize track----
  
  # Compile track from state-space using sampled cellnumbers
  obs <- statespace[s.grid,]
  obs <- as.data.frame(obs)
  
  # RETURN RESULT
  tracks[[j]] <- cbind(obs, id = which_telemetered_acs[j])
  
}

# Combine tracks into a df
telemetered_df = do.call(rbind, tracks) %>%
  mutate(times = 1:nrow(.))


#----Plot the resulting tracks----s

# Plot of state-space, surface, & track
ggplot() +
  geom_tile(data=as.data.frame(cost, xy=T), aes(x=x, y=y, fill = layer)) + 
  scale_fill_viridis("Cost", option = "D") + ggtitle("Individual track") +
  geom_path(data = telemetered_df, aes(x=x, y=y, group=id, color=times), size=0.2) +
  scale_color_viridis("Step number", option="C") +
  geom_point(data=as.data.frame(telemetered_acs), aes(x=x, y=y), 
             fill = "black", color="white", pch = 21) +
  coord_equal() + theme_minimal()


#----Thin the data----

# Degree of thinning
prop_thin <- 1  # percent for thinning

# Do the thinning
telemetered_thin_df <- telemetered_df %>%
  group_by(id) %>% # for each individual
  slice(., seq(1, nrow(.), length=prop_thin*nrow(.))) %>% # thin the data (take indicies)
  ungroup() # ungroup at the end


#----Plot the tracks and the thinned data per-pixel frequencies----

# Cost surface to df
df_cost_surface <- as.data.frame(landscape, xy=T)

# Raw tracks
p1 <- ggplot() +
  geom_tile(data=df_cost_surface, aes(x=x, y=y, fill = layer)) + 
  scale_fill_viridis("Cost", option = "D") + ggtitle("Individual track") +
  geom_path(data = telemetered_df, aes(x=x, y=y, color=times, group = id), size=0.2) +
  scale_color_viridis("Step number", option="C") +
  geom_point(data=as.data.frame(telemetered_acs), aes(x=x, y=y), 
             fill = "black", color="white", pch = 21) +
  coord_equal() + theme_minimal()

# Per-pixel frequency
p2 <- telemetered_thin_df %>%
  group_by(x, y) %>%
  summarise(count = n()) %>%
  ggplot(data = ., aes(x, y)) +
  geom_tile(aes(fill = count)) +
  ggtitle("Frequency of fixes per-pixel (n=90x24)") +
  scale_fill_viridis("Count") +
  coord_equal() +
  theme_minimal()

# Together
p1 + p2


#----Thinned tracks to matrix of pings per pixel----

# ppfreq template
ppFreq_template <- df_cost_surface %>%
  select(x, y) %>%
  mutate(count = 0)

# data collection matrix
n_by_cell_freq <- matrix(NA,nrow=telemetry_n,ncell(cost))

# ppfreq for each individual
for(i in 1:telemetry_n){
  
  tmp_ppfreq <- telemetered_thin_df %>% # use thinned data
    filter(id == which_telemetered_acs[i]) %>% # grab individual i
    group_by(x, y) %>% # group by xy coords to count
    summarise(count = n()) %>% # count within each group
    right_join(., ppFreq_template, by = c("x", "y")) %>% # right join to entire cost surface
    mutate(count = count.x + count.y) %>% # count the number of fixes per CS pixel
    mutate(count = if_else(is.na(count), 0, count)) %>% # NA + 0 to 0
    select(x, y, count) %>% # select important vars
    arrange(x) %>% # arrange/sort x and y to maintain order
    arrange(y) %>%
    pull(count) # pull the data
  
  # Fill in the matrix
  n_by_cell_freq[i,] <- tmp_ppfreq
  
}


#----Spatial encounter histories----

# (Cost) distance matrix between ACs and traps
d <- costDistance(tr1Corr, ac, traplocs)

# Capture probability
a1 <- 1/(2*sigma^2)                                              ############ SIGMA vs THETA? ################
probcap <- plogis(alpha2) * exp(-a1 * d^2)                       ############ IS PLOGIS RIGHT? ################

# Empty encounter data frame, data-collection matrix
Y <- matrix(NA, nrow=abundance, ncol=n_traps)

# Loop through each indvidual
for(i in 1:nrow(Y)){
  Y[i,] <- rbinom(n_traps, K, probcap[i,])
}

# Reduce to only captured individuals
Y <- Y[apply(Y,1,sum)>0,]

# Summary stats
nrow(Y) # Sample size n
sum(apply(Y, 1, sum)) # Captures per individual


