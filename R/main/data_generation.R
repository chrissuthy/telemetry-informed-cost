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

extract <- raster::extract
select <- dplyr::select

ref.Dmat_i <- function(from, to, local_ss_r, Dmat_i){
  from_cell <- as.numeric(raster::extract(x = local_ss_r, y = from, cellnumbers=T)[,1])
  to_cell <- as.numeric(raster::extract(x = local_ss_r, y = to, cellnumbers=T)[,1])
  result <- matrix(Dmat_i[from_cell, to_cell], nrow = nrow(from))
  return(result)
}

#----Simulation settings----

# Manual settings
nfix = 90*24
sims = 10

# Cost
alpha2 <- 2

# Movement model
psi <- 0.9
scale_factor_cost <- 4 # For scaling sigma later
upsilon <- 0.25 #0.6
upsilon_sf <- upsilon * scale_factor_cost
sigma <- 1  #4.3
sigma_sf <- sigma * scale_factor_cost

# SCR
N <- 50

# Statespace
ncol <- nrow <- 137 #125 for 3sig move buffer, ups=0.25, sig=1
rr <- upsilon # Actually, could this even by 2 ups? since ups_sf = ups*4
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
  aggregate(., fact = 4) %>%
  as.data.frame(., xy=T) %>%
  select(x, y) %>%
  filter(x >= (min(x)+hr95_lim) & x < (max(x)-hr95_lim)) %>%
  filter(y >= (min(y)+hr95_lim) & y < (max(y)-hr95_lim))


#----Initialize items to save---

landscape_ALL <- list()
teldata_raw_ALL <- list()
cost.data_ALL <- list()


#----Start sims----
t1 <- Sys.time()
tracks_all <- list()

# Parallel setup
ncores = detectCores()-1 # Number of available cores -1 to leave for computer
cl = makeCluster(ncores) # Make the cluster with that many cores
registerDoParallel(cl)  
clusterExport(cl, varlist = c("e2dist"), envir = environment()) # Export required function to the cores

for(sim in 1:sims){

  set.seed(sim)
  
  #----Landscape----
  
  landscape0 <- nlm_gaussianfield(
    ncol = ncol, nrow = nrow, 
    resolution = rr, autocorr_range = autocorr)
  landscape_r <- landscape0^2
  landscape <- as.data.frame(landscape_r, xy=T)
  
  # Save landscape to output
  landscape_ALL[[sim]] <- landscape_r

  
  #----Activity centers----
  
  acs <- ss[sample(1:nrow(ss), size = N, replace = T),c("x","y")]
  rownames(acs) <- NULL
  acs_df <- as.data.frame(acs)
  
  
  #----Plotting----
  
  # p1 <- ggplot() +
  #   geom_tile(data = landscape, aes(x=x, y=y, fill=layer)) +
  #   scale_fill_viridis_c("Landscape") +
  #   geom_rect(aes(xmin=min(ss$x), xmax=max(ss$x), ymin=min(ss$y), ymax=max(ss$y)),
  #             fill=alpha("white",0.25)) +
  #   geom_circle(data = acs_df, aes(x0=x, y0=y, r = hr95),
  #               lwd = 0.2, color = alpha("white", 0.65)) +
  #   coord_equal() +
  #   theme_minimal()
  # p1
  # text1 <- paste0(
  #   "N = ", N, "\n",
  #   "Psi = ", psi, "\n",
  #   "Upsilon = ", upsilon, "\n",
  #   "Sigma = ", sigma, "\n\n",
  #   "pixels = ", nrow(ss), "\n",
  #   "resolution = ", "upsilon"
  # )
  # 
  # p2 <- ggplot() +
  #   annotate("text", x = 0, y = 0, size=8, 
  #            label = text1) + 
  #   theme_void()
  # 
  # p1+p2
  
  
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
  
  # Cluster!
  t0 <- Sys.time()
  tracks <- foreach(j=1:N, .packages = c(.packages())) %dopar% {
    
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
      if(F){
        par(mfrow=c(2,2))
        
        # Landscape
        plot(local_ss_r, main = "Covariate")
        
        # Step dist
        plot(rasterFromXYZ(cbind(local_ss[,1:2], z = as.numeric(D))), main = "Upsilon eco")
        
        # AC dist
        plot(rasterFromXYZ(cbind(local_ss[,1:2], z = as.numeric(Dac))), main = "Sigma eco")
        
        # Kern
        plot(rasterFromXYZ(cbind(local_ss[,1:2], z = as.numeric(kern))), main = "Next step probabilities")
        lines(local_ss[s.grid,c("x","y")])
        
        par(mfrow=c(1,1))
      }
      
      # Sample a pixel for the next step using kernel
      s.grid[i]<- sample( 1:length(kern), 1, prob=kern)
    }
    
    #----Organize track----
    
    # Compile track from state-space using sampled cellnumbers
    obs <- local_ss[s.grid,c("x","y")]
    obs <- as.data.frame(obs)
    rownames(obs) <- NULL
    
    # RETURN RESULT
    tracks[[j]] <- cbind(obs, id = j)
    
  }
  
  # Post-cluster cleanup
  # stopCluster(cl)
  tf <- Sys.time()
  tf-t0
  
  # Combine tracks into a df
  telemetered_df = do.call(rbind, tracks) %>%
    as.data.frame() %>%
    group_by(id) %>%
    mutate(times = 1:n()) %>%
    ungroup()
  
  write.table(
    telemetered_df, 
    file = paste0("simout/sim", sim, ".txt")
  )
  
  tracks_all[[sim]] <- telemetered_df
  
  
  "COMPILE TELEM"

  #----Compile telemetry data----

  # Data-collection setup
  cost.data <- list()
  teldata_raw <- list()

  # Select k of n simulated individuals
  for(i in 1:N){

    # Subset track
    tmp_track <- tracks[[i]][,1:2]

    # trimS type buffer
    trimS <- 3*upsilon_sf

    # Get extent from track
    tmp_move_ext <- extent(c(min(tmp_track$x)-trimS, max(tmp_track$x)+trimS,
                             min(tmp_track$y)-trimS, max(tmp_track$y)+trimS))

    # Crop landscape to extent
    tmp_landscape_crop <- raster::crop(landscape_r, tmp_move_ext)

    # Assign individual-specific objects
    teldata_raw[[i]] <- tmp_track
    cost.data[[i]] <- as.matrix(raster::as.data.frame(tmp_landscape_crop, xy=T))

  }

  teldata_raw_ALL[[sim]] <- teldata_raw
  cost.data_ALL[[sim]] <- cost.data
  
  
}
stopCluster(cl)
t2 <- Sys.time()
t2-t1


# Plot of state-space, surface, & track
ggplot() +
  geom_tile(data=as.data.frame(cost, xy=T), aes(x=x, y=y, fill = layer)) +
  scale_fill_gradientn("alpha2*cost", colors = c("darkgreen", "lightgreen")) + ggtitle("Individual track") +
  geom_path(data = telemetered_df, aes(x=x, y=y, group=id, color = id), size=0.2) +
  scale_color_viridis("Individual", option="C") +
  geom_point(data=as.data.frame(acs), aes(x=x, y=y),
             fill = "black", color="white", pch = 21) +
  coord_equal() + theme_minimal()


#----SAVE DATA FOR MODELS----
saveRDS(ss,              "output/model_data/ss.RData")
saveRDS(teldata_raw_ALL, "output/model_data/teldata_raw.RData")
saveRDS(cost.data_ALL,   "output/model_data/cost_data.RData")
saveRDS(landscape_ALL,   "output/model_data/landscape.RData")
saveRDS(tracks_all,      "output/model_data/tracks_all.RData")  




