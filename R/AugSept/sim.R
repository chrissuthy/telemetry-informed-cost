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
#library(patchwork)
area <- raster::area

set.seed(202006)


"SIMULATION SETUP"

#----Create the landscape----

# Set autocorrelation range
autocorr <- 3

# Simulate landscape
gauss_nlm <- nlm_gaussianfield(ncol = 200, nrow = 200, resolution = 1,
                               autocorr_range = autocorr)

# Create state-space raster from landscape (axes are different)
use_this <- gauss_nlm
landscape0 <- rasterFromXYZ(cbind(coordinates(use_this)/10,values(use_this)))
landscape <- landscape0^2
statespace <- as.data.frame(landscape, xy=T)[,c("x", "y")] 
statespace <- as.matrix(statespace)


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

# MM parameters
sigma <- 0.7/6.5   # space use - step length
sbar <- ac[1,]     # Starting position, just one individual for now
psi <- 0.7         # This is pr(moved)
theta <- 0.7       # Space use - home range
alpha2 <- 1        # Cost value
nfix <- 90*24      # Number of pings
moved <- rbinom(nfix, 1, psi) # Vector of decisions to move


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


#---Start tracking individuals----

# Create vector of steps
s.grid <- rep(NA, nfix)
s.grid[1] <- extract(stack(landscape),matrix(sbar,nrow=1), cellnumber=T)[1]


# Loop to generate fix locations
for(i in 2:nfix){
  
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

# Convert to data frame with time
obs <- obs %>%
  as.data.frame() %>%
  mutate(times = 1:nrow(.))

# Plot of state-space, surface, & track
ggplot() +
  geom_tile(data=as.data.frame(cost, xy=T), aes(x=x, y=y, fill = layer)) + 
  scale_fill_viridis("Cost", option = "D") + ggtitle("Individual track") +
  geom_path(data = obs, aes(x=x, y=y, color=times), size=0.2) +
  scale_color_viridis("Step number", option="C") +
  #geom_point(data=as.data.frame(telemetered_acs), aes(x=X, y=Y), 
  #           fill = "black", color="white", pch = 21) +
  coord_equal() + theme_minimal() +
  xlim(min(obs$x)-1, max(obs$x)+1) + ylim(min(obs$y)-1, max(obs$y)+1)



