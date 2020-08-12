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
  
  #theta_probd[i] <- dwrappedcauchy(x = theta[i], mu=circular(0), rho=0.05)
  theta_probd[i] <- dnorm(x = theta[i], mean = 0, sd = 2.25) # Similar to dnorm
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

