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


"SIMULATION SETUP"

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


#---Create state-space---- (Though I don't really use this?)

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

# Abundance to density
abundance <- 300
area <- prod(round(apply(bbox(cost_surface),1,diff),2))

abs.density <- abundance / area
pix.density <- abundance / nrow(statespace)

# Activity centers
acs <- cbind(X = round(runif(abundance,0,10),2),
             Y = round(runif(abundance,0,10),2))

# Space-use paramter (step length?)
sigma <- 0.35

# Detection
p0 <- 0.1
K <- 5


#----Movement model parameters----

# Tag some of the individuals
telemetry_n <- 16
which_telemetered_acs <- sample(1:nrow(acs), telemetry_n)
telemetered_acs <- acs[which_telemetered_acs,]
n_fixes <- 90*24 # Every hour for three month study period


# Surface parameters  ############ ARE THESE RIGHT? ################
b_0 <- -1
b_lcp <- -5
b_ac <- 5



"MOVEMENT MODEL"


#----Start movement model----

# Cost surface
df_cost_surface <- as.data.frame(cost_surface, xy=T)
p0_telem <- 1  ############ I feel like this should be the same as above ################

# Data-collection object
telemetered_tracks <- list()

# Loop through each individual ac
for(i in 1:telemetry_n){
  
  # Grab individual i
  acs <- matrix(telemetered_acs[i,],1,2,byrow=T)
  colnames(acs) <- c("X", "Y")
  
  
  #----Home range bias----
  
  # Same pixels as cost surface
  dmat <- e2dist(x = acs, y = df_cost_surface[,c("x", "y")] )
  pmat <- p0_telem * exp(-dmat * dmat/(2 * sigma * sigma))
  p_ss <- cbind(df_cost_surface[,c("x", "y")], p = as.numeric(pmat))
  
  
  #----Start tracking!----
  
  # Initiate a tag
  steps = data.frame(
    x = acs[1,"X"],
    y = acs[1,"Y"],
    times = 0
  )
  
  # First step into a pixel centroid
  steps[2,] <- p_ss %>% arrange(desc(p)) %>% slice(1) %>% mutate(p = 1)
  
  # Set progress bar
  print(paste("Individual", i, "of", telemetry_n))
  pb = txtProgressBar(min = 1, max = n_fixes, initial = 1, style=3) 
  
  # Tag 'em and watch 'em run!
  for(j in 1:n_fixes){
    
    setTxtProgressBar(pb,j)
    
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
    
    # Write out track to data-collection object
    telemetered_tracks[[i]] <- cbind(steps, id = which_telemetered_acs[i])
    
    # Plot each iteration
    # xlim <- c(acs[1,1]-3*sigma, acs[1,1]+3*sigma)
    # ylim <- c(acs[1,2]-3*sigma, acs[1,2]+3*sigma)
    # plot(rasterFromXYZ(pi_df), ylim=ylim, xlim=xlim, col = viridis(1000)[350:1000])
    # points(acs, cex = 1.5, pch = 20, col = "white")
    # points(acs, cex = 1.2, pch = 20, col = "black")
    # lines(steps[,c("x","y")], lwd=0.2)
    # points(steps[nrow(steps),], col = "red", cex = 1.2, pch = 20)
    
  }
  
}


#----Thin the telemetry data----

# Prep the telemetry data
telemetered_df <- do.call(rbind, telemetered_tracks) %>%
  filter(times > 0) # getting rid of starting manual fixes.

# Thin the telemetry data
prop_thin <- 0.25  # percent for thinning
telemetered_thin_df <- telemetered_df %>%
  group_by(id) %>% # for each individual
  mutate(times = c(1:(max(times)))) %>% # just resetting the times
  slice(., seq(1, max(times), length = max(times)*prop_thin)) %>% # thin the data (take indicies)
  ungroup() # ungroup at the end


#----Plot the tracks and the thinned data per-pixel frequencies----

# Raw tracks
p1 <- ggplot() +
  geom_tile(data=df_cost_surface, aes(x=x, y=y, fill = layer)) + 
  scale_fill_viridis("Cost", option = "D") + ggtitle("Individual track") +
  geom_path(data = telemetered_df, aes(x=x, y=y, color=times, group = id), size=0.2) +
  scale_color_viridis("Step number", option="C") +
  geom_point(data=steps[1,], aes(x=x, y=y), fill = "black", color="white", pch = 21) +
  coord_equal() + theme_minimal()

# Per-pixel frequency
p2 <- telemetered_thin_df %>%
  group_by(x, y) %>%
  summarise(count = n()) %>%
  ggplot(data = ., aes(x, y)) +
  geom_tile(aes(fill = count)) +
  ggtitle("Frequency of per-pixel (n=90x24)") +
  scale_fill_viridis() +
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
n_by_cell_freq <- matrix(NA,nrow=telemetry_n,ncell(cost_surface))

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



