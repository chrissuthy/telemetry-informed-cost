# Gates Dupont
# September 2020

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


"LANDSCAPE"

#----Create the landscape----

# Set autocorrelation range
autocorr <- 3

# Simulate landscape
rr <- (0.5/4)
lncol <- lnrow <- 10/rr
gauss_nlm <- nlm_gaussianfield(ncol = lncol, nrow = lnrow, resolution = rr,
                               autocorr_range = autocorr)

# Create state-space raster from landscape (axes are different)
use_this <- gauss_nlm
landscape0 <- rasterFromXYZ(cbind(coordinates(use_this),values(use_this)))
landscape <- landscape0^2
statespace <- as.data.frame(landscape, xy=T)[,c("x", "y")] 
statespace <- as.matrix(statespace)


"SCR: population"

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


"MOVEMENT MODEL PARAMETERS"

#----Movement model parameters----

# MM parameters
alpha2 <- 1          # Cost value
sigma <- 0.5         # Space use - home range
upsilon <- sigma/4   # Space use - step length
psi <- 0.7           # This is pr(moved)
nfix <- 90*24        # Number of pings


"SCR: parameters"

#---SCR traps---

# Detection
p0 <- 0.2
K <- 5
d0 <- pix.density

# Traps - 2 sigma design, 4 sigma from edges (when sigma = 0.5)
traplocs <- as.matrix(
  expand.grid(X = seq(5-3*2*sigma,5+3*2*sigma,length=7),
              Y = seq(5-3*2*sigma,5+3*2*sigma,length=7)))

# Number of traps
n_traps <- nrow(traplocs)

# Plot
plot(landscape, col = viridis(1000))
points(traplocs, col = "white", pch = 21)
points(traplocs, col = "red", pch = 20)
points(ac, pch = 20, col = "black")


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


"DATA-GENERATION: CAMERAS"

#----Spatial encounter histories----

D <- costDistance(tr1Corr, ac, traplocs) # opposite, ac -> trap, not trap -> pixel
probcap <- p0 * exp(-(D*D)/(2*sigma*sigma))

Y <- matrix(NA, nrow = abundance, ncol = n_traps)
for(i in 1:nrow(Y)){
  Y[i,] <- rbinom(n_traps, K, probcap[i,])
}
Y <- Y[apply(Y,1,sum)>0,]
nrow(Y) # number of individuals
sum(apply(Y, 1, sum)) # Number of captures
table(apply(Y>0, 1, sum)) # Individuals captured on X traps (sp recaps)


"DATA-GENERATION: TELEMETRY"

#---Start tracking individuals----

# Tag some of the individuals
telemetry_n <- 8
which_telemetered_acs <- sample(1:nrow(ac), telemetry_n)
telemetered_acs <- ac[which_telemetered_acs,]

# Parallel setup
library(doParallel) # For parallelization
ncores = 4 # Number of available cores -1 to leave for computer
cl = makeCluster(ncores) # Make the cluster with that many cores
registerDoParallel(cl)  
#clusterExport(cl, varlist = c("e2dist"), envir = environment()) # Export required function to the cores

# Data-collection list
tracks = list()

# Cluster!
t0 <- Sys.time()
tracks <- foreach(j=1:telemetry_n, .packages = c(.packages())) %dopar% {
  
  # Generating movement based on pr(move)
  moved <- rbinom(nfix, 1, psi) # This should be different for each individual
  
  # Get starting position for the individual
  sbar <- telemetered_acs[j,]
  
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
    
    # Calculate ecological distance from the last position to each pixel (for upsilon)
    D <- costDistance(tr1Corr,  statespace[s.grid[i-1],], statespace)
    
    # Calculate euclidean distance from the activity center to each pixel (for sigma)
    Dac <- sqrt( (statespace[s.grid[1],1] - statespace[,1] )^2  + (statespace[s.grid[1],2] - statespace[,2])^2  )
    
    # Create a movement kernel using ecoD & upsilon and eucDac & sigma
    kern <- exp((-D*D/(2*upsilon*upsilon))  - (Dac*Dac/(2*sigma*sigma)) )
    
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

# Post-cluster cleanup
stopCluster(cl)
tf <- Sys.time()
tf-t0


# Combine tracks into a df
telemetered_df = do.call(rbind, tracks) %>%
  as.data.frame() %>%
  group_by(id) %>%
  mutate(times = 1:n()) %>%
  ungroup()


#----Plot the resulting tracks----s

# Plot of state-space, surface, & track
ggplot() +
  geom_tile(data=as.data.frame(cost, xy=T), aes(x=x, y=y, fill = layer)) + 
  scale_fill_gradientn("Cost", colors = c("darkgreen", "lightgreen")) + ggtitle("Individual track") +
  geom_path(data = telemetered_df, aes(x=x, y=y, group=id, color = id), size=0.2) +
  scale_color_viridis("Step number", option="C") +
  geom_point(data=as.data.frame(telemetered_acs), aes(x=x, y=y), 
             fill = "black", color="white", pch = 21) +
  coord_equal() + theme_minimal()


"MODEL"

#----Fit the model----

iterations = 1
n_tagged_inds_selected <- 4

source("R/gates_new/scr_move_cost_like.R")

# Data-collection setup
out_it <- matrix(NA, nrow = iterations, ncol = 6)
colnames(out_it) <- c("alpha2", "upsilon", "psi", "sigma", "p0", "d0")

# Grab data and fit the model
for(j in 1:iterations){
  
  # Data-collection setup
  cost.data <- list()
  teldata_raw <- list()
  h <- sample(x = 1:telemetry_n, size=n_tagged_inds_selected)
  
  # Select k of n simulated individuals
  for(i in 1:n_tagged_inds_selected){
    
    # Subset track
    tmp_track <- tracks[[h[i]]][,1:2]
    
    # trimS type buffer
    trimS <- 4*upsilon
    
    # Get extent from track
    tmp_move_ext <- extent(c(min(tmp_track$x)-trimS, max(tmp_track$x)+trimS,
                             min(tmp_track$y)-trimS, max(tmp_track$y)+trimS))
    
    # Crop landscape to extent
    tmp_landscape_crop <- crop(landscape, tmp_move_ext)
    
    # Assign individual-specific objects
    teldata_raw[[i]] <- tmp_track
    cost.data[[i]] <- as.matrix(as.data.frame(tmp_landscape_crop, xy=T))
    
  }
  
  # NLM likelihood evaluation
  t0 <- Sys.time()
  mm <- nlm(scr_move_cost_like, mod = "gauss",
            c(0, 0, 0, 0, 0, 0), hessian = T,
            teldata = teldata_raw, 
            spatdata = cost.data,
            landscape = landscape,
            K = K, scr_y = Y, trap_locs = traplocs,
            dist = "lcp", popcost=T, popmove=T, fixcost=F, use.sbar=T, prj=NULL)
  tf <- Sys.time()
  t_total <- tf-t0
  
  # Back-transform estimates
  est <- mm$estimate
  
  final <- c()
  final[1] <- est[1]
  final[2] <- exp(est[2])
  final[3] <- plogis(est[3])
  final[4] <- exp(est[4])
  final[5] <- plogis(est[5])
  final[6] <- exp(est[6])
  
  # Append
  out_it[j,] <- final
}

# Percent relative bias
rb <- data.frame(
  alpha2 = 100*((final[1]-alpha2)/alpha2),
  upsilon = 100*((final[2]-upsilon)/upsilon),
  psi = 100*((final[3]-psi)/psi),
  sigma = 100*((final[4]-sigma)/sigma),
  p0 = 100*((final[5]-p0)/p0),
  d0 = 100*((final[6]-d0)/d0)
)

# Print results
print(t_total)
print(rb)


