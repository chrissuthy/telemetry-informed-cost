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
psi <- 0.5
upsilon <- 250 #600
sigma <- 1000  #4300

# SCR
N <- 50

# Statespace
ncol <- nrow <- 131 #175
rr <- upsilon
autocorr <- 7

# Derived 
hr95 <- sqrt(5.99) * sigma
hr95_lim <- (3*sigma) + (3*upsilon) # this is 3 sigma


#----Start sims----
t1 <- Sys.time()
tracks_all <- list()

for(sim in 1:sims){
  
  set.seed(sim)

  #----Landscape----
  
  landscape0 <- nlm_gaussianfield(
    ncol = ncol, nrow = nrow, 
    resolution = rr, autocorr_range = autocorr)
  landscape_r <- landscape0^2 # I STILL DON'T LIKE THIS! RESTRICTS LOW COST
  landscape <- as.data.frame(landscape_r, xy=T)
  
  ss <- landscape %>%
    filter(x >= (min(x)+hr95_lim) & x < (max(x)-hr95_lim)) %>%
    filter(y >= (min(y)+hr95_lim) & y < (max(y)-hr95_lim))
  nrow(ss)
  
  
  #----Activity centers----
  
  acs <- ss[sample(1:nrow(ss), size = N, replace = T),c("x","y")]
  rownames(acs) <- NULL
  acs_df <- as.data.frame(acs)
  
  
  #----Plotting----
  
  # p1 <- ggplot() +
  #   geom_tile(data = landscape, aes(x=x, y=y, fill=layer)) +
  #   scale_fill_viridis_c("Landscape") +
  #   geom_rect(aes(xmin=min(ss$x), xmax=max(ss$x), ymin=min(ss$y), ymax=max(ss$y)), 
  #             fill=alpha("white",0.1)) +
  #   geom_circle(data = acs_df, aes(x0=x, y0=y, r = hr95), 
  #               lwd = 0.2, color = alpha("white", 0.65)) +
  #   coord_equal() +
  #   theme_minimal()
  # 
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
  
  # Parallel setup
  ncores = detectCores()-1 # Number of available cores -1 to leave for computer
  cl = makeCluster(ncores) # Make the cluster with that many cores
  registerDoParallel(cl)  
  clusterExport(cl, varlist = c("e2dist"), envir = environment()) # Export required function to the cores
  
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
    step_max <- (sqrt(5.99) * sigma) + 3*upsilon
    
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
    Dac <- ref.Dmat_i(
      from = sbar, 
      to = local_ss %>% select(x,y) %>% as.matrix(), 
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
      D <- ref.Dmat_i(
        from = local_ss %>% filter(cell == s.grid[i-1]) %>% select(x,y) %>% as.matrix(),
        to = local_ss %>% select(x,y) %>% as.matrix(),
        local_ss_r,
        Dmat_i
      )
      
      # Distance from activity center only needs to be calculated once (above)
      
      # Create a movement kernel using ecoD & upsilon and eucDac & sigma
      kern <- exp((-D*D/(2*upsilon*upsilon))  - (Dac*Dac/(2*sigma*sigma)) )
      
      # Assign the current position a probability of zero (can't stay in the same spot, b/c Psi=1)
      kern[s.grid[i-1]] <- 0
      
      # Normalize the kernel cell values into probabilities
      kern <- kern/rowSums(kern)
      
      # Plot all of this
      if(T){
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
  stopCluster(cl)
  tf <- Sys.time()
  tf-t0
  
  # Combine tracks into a df
  telemetered_df = do.call(rbind, tracks) %>%
    as.data.frame() %>%
    group_by(id) %>%
    mutate(times = 1:n()) %>%
    ungroup()
  
  tracks_all[[sim]] <- telemetered_df
    
  
}
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
 
#saveRDS(telemetered_df, file = "output/oct13_N50_alpha2of2.RData")  