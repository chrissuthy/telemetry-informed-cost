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


#----Create the landscape----

# Set autocorrelation range
autocorr <- 3

# Simulate landscape
gauss_nlm <- nlm_gaussianfield(ncol = 100, nrow = 100, resolution = 1,
                               autocorr_range = autocorr)

# Create state-space raster from landscape (axes are different)
use_this <- gauss_nlm
landscape <- rasterFromXYZ(cbind(coordinates(use_this)/10,
                                 values(use_this)))


#----Create cost surface----

# Set cost surface parameter
cost_parameter <- 1.5

# Calculate per-pixel cost values (surface)
cost_surface <- exp(cost_parameter * landscape^2)

# Create tansition surface
transistion_surface <- geoCorrection(
  transition(
    cost_surface,
    transitionFunction = function(x) (1/(mean(x))),
    direction = 16),
  scl = F)


#---Create state-space----

# Set resolution
res <- 0.125

# Create state-space
statespace <-expand.grid(
    X = seq(round(bbox(cost_surface)[1,1])+res/2,
            round(bbox(cost_surface)[1,2])-res/2,
            res),
    Y = seq(round(bbox(cost_surface)[2,1])+res/2,
            round(bbox(cost_surface)[2,2])-res/2,
            res))


#----SCR parameters----

# Activity centers
acs <- cbind(X = round(runif(1,0,10),2),
             Y = round(runif(1,0,10),2))

# Space-use paramter (step length?)
sigma <- 0.35


#----Cost surface----

# Points for LCP
x1y1 <- matrix(acs[1,],1,2,byrow=TRUE)
x2y2 <- coordinates(cost_surface)

# Calculate LCP per-pixel from AC
telem_rsf_dmat <- costDistance(
  transistion_surface,
  fromCoords = x1y1,
  toCoords = x2y2)

# Create distance surface based on sigma
telem_rsf_lammat <- exp(1-telem_rsf_dmat^2/(2*sigma^2))

# Create raster from cost and sigma
r_dlcp <- rasterFromXYZ(cbind(coordinates(cost_surface),cost = c(telem_rsf_lammat)))

# Clip this raster
r_dlcp_ext <- extent(
  acs[1,1]-3*sigma, acs[1,1]+3*sigma,
  acs[1,2]-3*sigma, acs[1,2]+3*sigma)

# Crop the raster
r_dlcp <- crop(r_dlcp, r_dlcp_ext)

# d_lcp raster to df for plotting
df_dlcp <- as.data.frame(r_dlcp, xy=T)

# Plot dlcp surface
ggplot(data = df_dlcp, aes(x = x, y = y)) +
  geom_raster(aes(fill = cost)) +
  geom_point(data = as.data.frame(acs), 
             aes(x = X, y = Y), color = "red") +
  coord_equal() +
  scale_fill_viridis()


#----Per-pixel turning angles----

# Just make a first step to test TrajAngles
step_1 <- c(acs[,1]+0.2, acs[,2]-0.1)
steps <- as.data.frame(rbind(acs, step_1))
steps$times <- 0:(nrow(steps)-1)
colnames(steps) <- c("x", "y", "times")
rownames(steps) <- NULL

# Plot d_lcp and first step
plot(r_dlcp, col=viridis(1000))
points(acs, pch = 20, cex = 2, col = "red")
lines(steps, lwd=2)

# Loop through pixels, calculate angle from x1y1 for each
theta <- c()
theta_probd <- c()
for(i in 1:nrow(df_dlcp)){
  step_start <- nrow(steps)-1
  step_end <- nrow(steps)
  
  # Bind last step with each possible step
  tmp_coords <- rbind(
    steps[step_start:step_end,c("x","y")],
    df_dlcp[i,c("x","y")])
  rownames(tmp_coords) <- NULL
  
  tmp_coords$times <- 0:(nrow(tmp_coords)-1)
  
  # Make this a trj object
  trj <- TrajFromCoords(tmp_coords)
  
  # Calculate angle
  theta[i] <- TrajAngles(trj = trj)
  
  theta_probd[i] <- dwrappedcauchy(x = theta[i], mu=circular(0), rho=0.2)
  #theta_probd[i] <- dnorm(x = theta[i], mean = 0, sd = 0.3) # Similar to dnorm
}

# Visual check
plot(rasterFromXYZ(cbind(df_dlcp[,1:2], theta)))
points(acs, pch = 20, cex = 2, col = "red")
lines(steps, lwd=2)


#----Compare the two----

# Plot
par(mfrow=c(1,2))

# Cost Surface
plot(r_dlcp, col=viridis(1000), main = "Cost relative w.r.t AC")
points(acs, pch = 20, cex = 3, col = "white")
lines(steps, lwd=4, col="white")
points(acs, pch = 20, cex = 2, col = "black")
lines(steps, lwd=2, col="black")

# Turning angle
plot(rasterFromXYZ(cbind(df_dlcp[,1:2], theta)), 
     col = plasma(1000), main = "Turning angle w.r.t. step 1")
points(acs, pch = 20, cex = 3, col = "white")
lines(steps, lwd=4, col="white")
points(acs, pch = 20, cex = 2, col = "black")
lines(steps, lwd=2, col="black")

# End plot
par(mfrow=c(1,1))


#----Converting to probability rasters----

b0_cost <- -1
b1_cost <- 1
prob_cost <- exp(b0_cost + b1_cost*df_dlcp$cost)/sum(exp(b0_cost + b1_cost*df_dlcp$cost))

# Plot
par(mfrow=c(1,2))

# Cost Surface probabilities
plot(rasterFromXYZ(cbind(df_dlcp[,1:2], prob_cost)),
     col=viridis(1000), main = "Cost probabilities")
points(acs, pch = 20, cex = 3, col = "white")
lines(steps, lwd=4, col="white")
points(acs, pch = 20, cex = 2, col = "black")
lines(steps, lwd=2, col="black")

# Turning angle probabilities
plot(rasterFromXYZ(cbind(df_dlcp[,1:2], theta_probd)),
     col=plasma(1000), main = "Turning angle probabilities")
points(acs, pch = 20, cex = 3, col = "white")
lines(steps, lwd=4, col="white")
points(acs, pch = 20, cex = 2, col = "black")
lines(steps, lwd=2, col="black")

# Reset plot
par(mfrow=c(1,1))


#----Get next step----

# Multiplicative
r_nextStep <- rasterFromXYZ(cbind(df_dlcp[,1:2], prob_cost)) * rasterFromXYZ(cbind(df_dlcp[,1:2], theta_probd))

# Plot this to see it
plot(r_nextStep, col = viridis(1000))

# Make it a df to sample
df_nextStep <- as.data.frame(r_nextStep, xy=T)

# Sample according to multiplicative probabilities
idx <- base::sample(x = nrow(df_nextStep), size = 1, prob = df_nextStep[,"layer"], replace=T)
nextStep <- df_nextStep[idx,1:2]

# Plot the sampled point
points(nextStep, pch = 20, cex=2, col="red")

# If we do this a ton, can we recover the distribution?
idx <- base::sample(x = nrow(df_nextStep), size = 500000, prob = df_nextStep[,"layer"], replace=T)
nextStep <- df_nextStep[idx,1:2]

# Yes, we can!
nextStep %>%
  group_by(x, y) %>%
  summarise(count = n()) %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = count)) +
  scale_fill_viridis() +
  coord_equal()

# But it's adding some noise to the ladscape due to the turning angle


#----Make this iterative----

# Just make a first step to test TrajAngles
step_1 <- df_dlcp %>% arrange(desc(cost)) %>% slice(1) %>% select(x,y) %>% matrix(1,2)
steps <- as.data.frame(rbind(acs, step_1))
steps$times <- 0:(nrow(steps)-1)
colnames(steps) <- c("x", "y", "times")
rownames(steps) <- NULL

# Plot d_lcp and first step
plot(r_dlcp, col=viridis(1000))
points(acs, pch = 20, cex = 2, col = "red")
lines(steps, lwd=2, col = "red")

# Loop through pngs
pings <- 90*24*10
for(i in 1:pings){
  
  cat(i, "\n")
  
  # Loop through ss pixels, calculate angle from x1y1 for each
  theta <- c()
  theta_probd <- c()
  for(j in 1:nrow(df_dlcp)){
    
    # Last step for relative turning angle
    step_start <- nrow(steps)-1
    step_end <- nrow(steps)
    
    # Bind last step with each possible step
    tmp_coords <- rbind(
      steps[step_start:step_end,c("x","y")],
      df_dlcp[j,c("x","y")])
    rownames(tmp_coords) <- NULL
    
    tmp_coords$times <- 0:(nrow(tmp_coords)-1)
    
    # Make this a trj object
    trj <- TrajFromCoords(tmp_coords)
    
    # Calculate angle
    theta[j] <- TrajAngles(trj = trj)
    
    # Calculate probability density of the relative turning angle surface
    theta_probd[j] <- dwrappedcauchy(x = theta[j], mu=circular(0), rho=0.2)
    #theta_probd[j] <- dnorm(x = theta[j], mean = 0, sd = 0.3) # Similar to dnorm
  }
  
  ######## Step-specific least-cost path based on activity center specific cost surface. Still in progress.
  
  # # Create tansition surface
  # tmp_transistion_surface <- geoCorrection(
  #   transition(
  #     r_dlcp,
  #     transitionFunction = function(x) (1/(mean(x))),
  #     direction = 16),
  #   scl = F)
  # 
  # # Points for LCP
  # tmp_x1y1 <- steps %>% slice(nrow(.)) %>% select(x,y) %>% as.matrix() %>% unlist() %>% as.numeric()
  # tmp_x2y2 <- coordinates(r_dlcp)
  # 
  # # Calculate LCP per-pixel from AC
  # tmp_telem_rsf_dmat <- costDistance(
  #   tmp_transistion_surface,
  #   fromCoords = tmp_x1y1,
  #   toCoords = tmp_x2y2)
  # 
  # # Final: per-pixel cost distance based on current step.
  # tmp_r_dlcp <- rasterFromXYZ(cbind(df_dlcp[,c("x", "y")], t(tmp_telem_rsf_dmat)))
  # #plot(tmp_r_dlcp)
  # 
  # # Convert to data frame and convert to probabilities
  # tmp_df_dlcp <- as.data.frame(tmp_r_dlcp, xy=T)
  # colnames(tmp_df_dlcp) <- c("x", "y", "cost")
  # 
  # b0 <- -1
  # b1 <- 1
  # tmp_prob_stepCost <- 1 - exp(b0 + b1*tmp_df_dlcp$cost)/sum(exp(b0 + b1*tmp_df_dlcp$cost))
  
  ########
  
  # Calculate probability density of next step using multiplicative relationship between cost and theta
  cost_nextStep <- cbind(df_dlcp[,1:2], prob_cost) # Original, using AC-specific cost
  #cost_nextStep <- cbind(df_dlcp[,1:2], tmp_prob_stepCost) # swapping out prob_cost here with step-specific perspective of LCP
  theta_nextStep <- cbind(df_dlcp[,1:2], theta_probd)
  r_nextStep <- rasterFromXYZ(cost_nextStep) * rasterFromXYZ(theta_nextStep)
  
  # Make it a df to sample
  df_nextStep <- as.data.frame(r_nextStep, xy=T)
  colnames(df_nextStep) <- c("x", "y", "prob")
  
  # Sample according to multiplicative probabilities
  idx <- base::sample(x = nrow(df_nextStep), size = 1, prob = df_nextStep[,"prob"], replace=F)
  nextStep <- df_nextStep[idx,1:2]
  nextStep$times <- steps$times[nrow(steps)] + 1
  
  # Append selected step to next step
  steps <- rbind(steps, nextStep)
  
}

steps2 <- data.frame(
  x = unlist(steps$x),
  y = unlist(steps$y),
  times = unlist(steps$times)
)

# Plot d_lcp and first step
plot(r_dlcp, col=viridis(1000))
points(acs, pch = 20, cex = 2, col = "red")
lines(steps2, lwd=2, col = "red")

# Plot dlcp surface with tracks
ggplot(data = df_dlcp, aes(x = x, y = y)) +
  geom_raster(aes(fill = cost)) +
  geom_point(data = as.data.frame(acs), 
             aes(x = X, y = Y), color = "red") +
  geom_path(data=steps2, aes(x=x, y=y), color = "red", lwd=0.2) +
  coord_equal() +
  scale_fill_viridis()

# Checking if we can recover original distribution
p2 <- steps2[3:nrow(steps2),] %>% # Get rid of manually-selected first pings to clean up raster
  group_by(x, y) %>%
  summarise(count = n()) %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = count)) +
  ggtitle("Frequency of per-pixel (n=90x24x10)") +
  scale_fill_viridis() +
  coord_equal()

# Oirignal distribution
p1 <- ggplot(data = df_dlcp, aes(x = x, y = y)) +
  geom_raster(aes(fill = cost)) +
  ggtitle("Cost surface") +
  coord_equal() +
  scale_fill_viridis()

p1+p2
