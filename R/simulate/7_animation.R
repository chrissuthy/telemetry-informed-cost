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


library(NLMR)
library(dplyr)
library(plotrix)
library(raster)
library(ggplot2)
library(viridis)
library(ggforce)
library(patchwork)
library(doParallel)
library(gdistance)
library(oSCR)
library(ggnewscale)
library(magick)

extract <- raster::extract
select <- dplyr::select

ref.Dmat_i <- function(from, to, local_ss_r, Dmat_i){
  from_cell <- as.numeric(raster::extract(x = local_ss_r, y = from, cellnumbers=T)[,1])
  to_cell <- as.numeric(raster::extract(x = local_ss_r, y = to, cellnumbers=T)[,1])
  result <- matrix(Dmat_i[from_cell, to_cell], nrow = nrow(from))
  return(result)
}


"STARTING PARAMS"

#----Simulation settings----

# Manual settings
nfix = 100
sims = 100

# Cost
alpha2 <- 2

# Movement model
psi <- 1
scale_factor_cost <- 4  #a2 = 1, sf = 2.5 OR a2 = 2, sf = 4
sigma <- 1
sigma_sf <- sigma * scale_factor_cost
small_ups <- 0.25
small_ups_sf <- small_ups*scale_factor_cost
upsilon <- ifelse(select_ups == "small ups", small_ups, sigma)
upsilon_sf <- upsilon * scale_factor_cost


# SCR
N <- 100

# Statespace
ncol <- nrow <- 137 #125 for 3sig move buffer, ups=0.25, sig=1
#rr <- upsilon # Actually, could this even by 2 ups? since ups_sf = ups*4
rr <- small_ups
autocorr <- 6

# Derived 
hr95 <- sqrt(5.99) * sigma
hr95_lim <- (3*sigma) * 1.5 # THIS MAY NEED UPDATING

#----GET SS----

# Do this here so we only need to get ss once
for_ss <- nlm_gaussianfield(
  ncol = ncol, nrow = nrow, 
  resolution = rr, autocorr_range = autocorr)

# Make ss using aggregated pixels
ss <- for_ss %>%
  #aggregate(., fact = 4) %>% # THIS MESSES UP TRAPS AND ISNT USED BECAUSE RESOLUTION IS DIFFERENT!!!
  as.data.frame(., xy=T) %>%
  select(x, y) %>%
  filter(x >= (min(x)+hr95_lim) & x < (max(x)-hr95_lim)) %>%
  filter(y >= (min(y)+hr95_lim) & y < (max(y)-hr95_lim))


"START SIM"

#----START SIM----

set.seed(1)

#----Landscape----

landscape0 <- nlm_gaussianfield(
  ncol = ncol, nrow = nrow, 
  resolution = rr, autocorr_range = autocorr)
landscape_r <- landscape0^2
landscape <- as.data.frame(landscape_r, xy=T)


#----Activity centers----

acs <- ss[sample(1:nrow(ss), size = N, replace = T),c("x","y")]
rownames(acs) <- NULL
acs_df <- as.data.frame(acs)


"COST SURFACE"

#----Cost surface----

# Cost surface
logcost <- alpha2*landscape_r
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


"SIMULATE MOVEMENT"

# Data-collection list
tracks = list()
j = 20

# Generating movement based on pr(move)
moved <- rbinom(nfix, 1, psi) # This should be different for each individual

# Get starting position for the individual
sbar <- acs[j,]
sbar_x <- as.numeric(sbar[1])
sbar_y <- as.numeric(sbar[2])
#step_max <- ((sqrt(5.99) * sigma) + 3*upsilon) * 1.5###### OLD
step_max <- (sqrt(5.99) * sigma_sf) * 1.05 # 1.5x b/c sig_sf NEW


# Subset statespace for local evaluation
local_ss <- landscape %>%
  filter(x <= (sbar_x+step_max) & x >= (sbar_x-step_max)) %>%
  filter(y <= (sbar_y+step_max) & y >= (sbar_y-step_max)) %>%
  mutate(cell_true = cellFromXY(landscape_r, xy = .)) %>%
  mutate(cell = 1:n())

# Make local ss into raster
local_ss_r <- rasterFromXYZ(local_ss[,1:3])

# Check this out
if(F){
  plot(local_ss_r)
  points(sbar)
  draw.circle(x = sbar_x, y = sbar_y, radius = sqrt(5.99)*sigma)
}

# Create Indvidual D matrix (from = rows, to = cols)
Dmat_i <- costDistance(
  tr1Corr,  
  local_ss %>% select(x,y) %>% as.matrix(),
  local_ss %>% select(x,y) %>% as.matrix())

# Create vector of steps
s.grid <- rep(NA, nfix)
s.grid[1] <- extract(x = local_ss_r, y = sbar, cellnumber=T)[1]

# Calculate ecological distance from the activity center to each pixel (for sigma)
# This only has to happen once
dac_to <- as.matrix(local_ss[,c("x","y")])

Dac <- ref.Dmat_i(
  from = sbar, 
  to = dac_to, 
  local_ss_r, Dmat_i)


# Loop to generate fix locations
for(i in 2:nfix){
  
  # If the animal doesn't move, assign same loc and skip the rest
  if(moved[i]==0){
    s.grid[i] <- s.grid[i-1]
    next # Next mean skip the rest of the for loop
  }
  
  ## If the animal moved...
  
  # Calculate ecological distance from the last position to each pixel (for upsilon)
  d_from <- as.matrix(local_ss[local_ss$cell == s.grid[i-1], c("x","y")])
  d_to <- as.matrix(local_ss[,c("x","y")])
  
  D <- ref.Dmat_i(
    from = d_from,
    to = d_to,
    local_ss_r,
    Dmat_i
  )
  
  # Distance from activity center only needs to be calculated once (above)
  
  # Create a movement kernel using ecoD & upsilon and eucDac & sigma
  kern <- exp((-D*D/(2*upsilon_sf*upsilon_sf))  - (Dac*Dac/(2*sigma_sf*sigma_sf)) )
  
  # Assign the current position a probability of zero (can't stay in the same spot, b/c Psi=1)
  kern[s.grid[i-1]] <- 0
  
  # Normalize the kernel cell values into probabilities
  kern <- kern/rowSums(kern)
  
  # Plot all of this
  if(T){
    
    # jpeg(plot_file, quality = 100, width = 540, res = 100)
    # 
    # par(mfrow=c(2,2))
    
    # # Landscape
    # plot(local_ss_r, main = "Covariate")
    
    
    p1 <- ggplot(data = local_ss, aes(x = x, y = y, fill = layer)) +
      geom_tile() +
      scale_fill_viridis(option = "C") +
      ggtitle("Surface") +
      coord_equal() +
      theme_void() +
      theme(legend.position = "none",
            plot.title = element_text(hjust=0.5))
    
    # # Step dist
    # plot(rasterFromXYZ(cbind(local_ss[,1:2], z = as.numeric(D))), main = "Upsilon eco")
    
    p2 <- ggplot(data = cbind(local_ss[,1:2], z = as.numeric(D)), 
           aes(x = x, y = y, fill = z)) +
      geom_tile() +
      scale_fill_viridis(option = "D", direction = -1) +
      ggtitle("Step") +
      coord_equal() +
      theme_void() +
      theme(legend.position = "none",
            plot.title = element_text(hjust=0.5))
    
    # # AC dist
    # plot(rasterFromXYZ(cbind(local_ss[,1:2], z = as.numeric(Dac))), main = "Sigma eco")
    
    p3 <- ggplot(data = cbind(local_ss[,1:2], z = as.numeric(Dac)), 
           aes(x = x, y = y, fill = z)) +
      geom_tile() +
      scale_fill_viridis(option = "D", direction = -1) +
      ggtitle("Activity center") +
      coord_equal() +
      theme_void() +
      theme(legend.position = "none",
            plot.title = element_text(hjust=0.5))
    
    # # Kern
    # plot(rasterFromXYZ(cbind(local_ss[,1:2], z = as.numeric(kern))), main = "Next step probabilities")
    # points(sbar, pch = 20)
    # lines(local_ss[s.grid,c("x","y")], col = alpha("black", 0.3))
    
    path <- local_ss[s.grid,c("x","y")] %>%
      na.omit() %>%
      mutate(step = 1:n())
    pos <- path %>%
      slice(nrow(.))
    
    p4 <- ggplot() +
      geom_tile(data = cbind(local_ss[,1:2], z = as.numeric(kern)), 
                aes(x = x, y = y, fill = z)) +
      scale_fill_viridis(option = "E", direction = 1) +
      geom_point(data=sbar, aes(x=x, y=y), color = "white", size=3) +
      geom_point(data=sbar, aes(x=x, y=y), color = "red", size=2.5) +
      new_scale("color") +
      geom_path(data = path, aes(x=x, y=y, color = step), size=0.5) +
      scale_color_gradientn(colors = c(
        alpha("white", 0),alpha("white", 0.2),alpha("black", 0.5))) +
      geom_point(data=pos, aes(x=x, y=y), pch = 20, color = "white", size = 2) +
      geom_point(data=pos, aes(x=x, y=y), pch = 20, color = "red", size = 1) +
      ggtitle("Next step prob.") +
      coord_equal() +
      theme_void() +
      theme(legend.position = "none",
            plot.title = element_text(hjust=0.5))
    
      # Compile final plot
      final_plot <- (p1 + p2) / (p3 + p4)
      
      # Save plot
      plot_file <- paste0("sim_", sprintf("%03d",i-1), ".jpeg")
      plot_path <- paste0("./output/animation/", select_ups)
      
      ggsave(
        filename = plot_file,
        path = plot_path,
        plot = final_plot,
        device = "jpeg",
        width = 6,
        height = 6,
        units = c("in"),
        dpi = 300)
  }
  
  # Sample a pixel for the next step using kernel
  s.grid[i]<- sample( 1:length(kern), 1, prob=kern)
}
  

if(F){
  
  plot_path <- paste0("./output/animation/", select_ups)
  
  ## list file names and read in
  imgs <- list.files(plot_path, full.names = TRUE)
  img_list <- lapply(imgs, image_read)
  
  ## join the images together
  img_joined <- image_join(img_list)
  
  ## animate at 2 frames per second
  img_animated <- image_animate(img_joined, fps = 2)

  ## image path
  img_file <- paste0(plot_path, "/track.gif")
  
  ## save to disk
  image_write(image = img_animated,
              path = img_file)
  
}
