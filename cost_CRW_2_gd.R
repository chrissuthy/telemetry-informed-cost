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


#----Home range bias----

# Same pixels as cost surface
df_cost_surface <- as.data.frame(cost_surface, xy=T)
p0 <- 1
dmat <- e2dist(x = acs, y = df_cost_surface[,c("x", "y")] )
pmat <- p0 * exp(-dmat * dmat/(2 * sigma * sigma))
p_ss <- cbind(df_cost_surface[,c("x", "y")], p = as.numeric(pmat))


#----Start tracking!----

# Surface parameters
b_0 <- -1
b_lcp <- -5
b_ac <- 5

# Initiate a tag
steps = data.frame(
  x = acs[1,"X"],
  y = acs[1,"Y"],
  times = 0
)

# First step into a pixel centroid
steps[2,] <- p_ss %>% arrange(desc(p)) %>% slice(1) %>% mutate(p = 1)

# Number of pings
n_pings = 90*24 # Ping every hour for 90 days

# Set progress bar
pb = txtProgressBar(min = 1, max = n_pings, initial = 1, style=3) 

# Tag 'em and watch 'em run!
for(i in 1:n_pings){
  
  setTxtProgressBar(pb,i)
  
  #---Home range bias----
  
  # p_ss # Constant
  
  
  #----LCP----
  
  # Points for LCP
  x1y1 <- as.numeric(matrix(steps[nrow(steps),1:2],1,2,byrow=TRUE))
  x2y2 <- coordinates(cost_surface)
  
  # Calculate LCP per-pixel from AC
  telem_rsf_dmat <- costDistance(
    transistion_surface,
    fromCoords = x1y1,
    toCoords = x2y2)
  
  # Convert LCP data to data frame
  df_dlcp <- as.data.frame(cbind(coordinates(cost_surface), dlcp = as.numeric(telem_rsf_dmat)))
  
  
  #----Combine LCP + AC----
  
  # Linear transform
  X_lcp <- df_dlcp$dlcp
  X_ac <- p_ss$p
  
  pi <- exp(b_0 + b_lcp*X_lcp + b_ac*X_ac)/sum(exp(b_0 + b_lcp*X_lcp + b_ac*X_ac))
  
  pi_df <- as.data.frame(cbind(coordinates(cost_surface), pi = pi))
  
  
  #----Sample position t+1----
  
  # Sample according to multiplicative probabilities
  idx <- base::sample(x = nrow(pi_df), size = 1, prob = pi_df[,"pi"])
  nextStep <- pi_df[idx,1:2]
  nextStep$times <- steps$times[nrow(steps)] + 1
  
  # Append selected step to next step
  steps <- rbind(steps, nextStep)
  
  # Plot each iteration
  # xlim <- c(acs[1,1]-3*sigma, acs[1,1]+3*sigma)
  # ylim <- c(acs[1,2]-3*sigma, acs[1,2]+3*sigma)
  # plot(rasterFromXYZ(pi_df), ylim=ylim, xlim=xlim, col = viridis(1000)[350:1000])
  # points(acs, cex = 1.5, pch = 20, col = "white")
  # points(acs, cex = 1.2, pch = 20, col = "black")
  # lines(steps[,c("x","y")], lwd=0.2)
  # points(steps[nrow(steps),], col = "red", cex = 1.2, pch = 20)
  
  
}


#----Plotting results----

# Plot track
p1 <- ggplot() +
  geom_tile(data=df_cost_surface, aes(x=x, y=y, fill = layer)) + 
  scale_fill_viridis("Cost", option = "D") + ggtitle("Individual track") +
  geom_path(data = steps, aes(x=x, y=y, color=times), size=0.2) +
  scale_color_viridis("Step number", option="C") +
  geom_point(data=steps[1,], aes(x=x, y=y), fill = "black", color="white", pch = 21) +
  coord_equal() + theme_minimal()


# Count pings per cell
ppFreq <- steps[3:nrow(steps),] %>% # Get rid of manually-selected first pings to clean up raster
  group_by(x, y) %>%
  summarise(count = n())

p2 <- ggplot(data = ppFreq, aes(x, y)) +
  geom_tile(aes(fill = count)) +
  ggtitle("Frequency of per-pixel (n=90x24)") +
  scale_fill_viridis() +
  coord_equal() +
  theme_minimal()

# Together
p1 + p2


#----Thin the telemetry data----

# Get rid of first couple of manually-selected positions
steps <- steps[steps$times > 2,]

# Calculate thinning rate
n_fixes <- nrow(steps)
prop_thin <- 0.2
n_thin <- n_fixes*prop_thin

# Thin the data
final_data <- steps[seq(1, nrow(steps), length = n_thin),]

# Plot the thhinned data
final_data %>% # Get rid of manually-selected first pings to clean up raster
  group_by(x, y) %>%
  summarise(count = n()) %>% 
  ggplot(data = ., aes(x, y)) +
  geom_tile(aes(fill = count)) +
  ggtitle("Frequency of per-pixel (n=90x24)") +
  scale_fill_viridis() +
  coord_equal() +
  theme_minimal()











