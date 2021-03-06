# THIS FILE IS TO TEST DAN's FIX TO THE LIKELIHOOD TO ALLOW IRREGULAR SCR STATE_SPACES (scr_ss)
# THIS WAS A SUCCESS!POINT ESTIMATES REMAINED THE SAME.

# # # # # # # # # # # # # # #
#                           #
#  ____HEY!!_DO THIS!!____  #
#     Parameterize sims     #
#                           #
# # # # # # # # # # # # # # #

# Right here:
select_ups  <- c("small ups", "big ups")[1]
select_ntel <- c(1, 3, 5)[1]
share_sig   <- c(TRUE, FALSE)[2]

# # # # # # # # # # # # # # #
#                           #
#       GREAT, THANKS!      #
#                           #
# # # # # # # # # # # # # # #


#----Get started----

library(gdistance)
library(dplyr)
# library(doParallel)


#----Directory setup----

# Create that directory
file_id <- paste0("est_ntel=", select_ntel, "_", "share=", share_sig)
new_dir <- paste0("./output/", select_ups, "/", file_id)
# if(!dir.exists(new_dir)){dir.create(new_dir)}


#----Load data----

# Telemetry data
file      <- paste0("./output/", select_ups, "/model_data/teldata_raw.RData")
teldata   <- readRDS(file)

# Landscape
file      <- paste0("./output/", select_ups, "/model_data/landscape.RData")
landscape <- readRDS(file)

# Y array
file      <- paste0("./output/", select_ups, "/model_data/y.RData")
y         <- readRDS(file)

# Traps
file      <- paste0("./output/", select_ups, "/model_data/traps.RData")
traps     <- readRDS(file) %>% as.matrix()
colnames(traps) <- c("X", "Y")

# SS
file      <- paste0("./output/", select_ups, "/model_data/ss.RData")
ss        <- readRDS(file) # colnames?

# Sampling occasions
K         <- 90

# Make SCR state-space
# This is just the 10,000 cells from 100x100, so it's NOT aggregated at fact = 4,
# and the likelihood does not split up the state-space and the cost landscape
scr_ss <- list()
for(i in 1:length(landscape)){
  
  scr_ss[[i]] <- crop(landscape[[i]], extent(ss)) %>%
    raster::aggregate(fact = 4)
}


#----Number of tagged individuals----

inds <- 1:select_ntel


#----Items for fitting movememt model----

# Number of sims
sims <- length(y)

# Additionals
if(select_ups == "small ups"){
  
  if(share_sig == TRUE){
    
    p <- c(1,                         # alpha2
           log(2.5*0.25),             # ups
           qlogis(0.9),               # psi
           log(2.5*1),                # sig
           qlogis(0.1),               # p0
           log(100/ncell(scr_ss[[1]])) # d0
    )
    
    out <- matrix(NA, nrow = sims, ncol = 6)
    colnames(out) <- c("alpha2", "upsilon", "psi", "sig", "p0", "d0")
    
  }else if(share_sig == FALSE){
    
    p <- c(1,                           # alpha2
           log(2.5*0.25),               # ups
           qlogis(0.9),                 # psi
           log(2.5*1),                  # sig
           qlogis(0.1),                 # p0
           log(100/ncell(scr_ss[[1]])), # d0
           log(2.5*1)                   # sig_mm
    )
    
    out <- matrix(NA, nrow = sims, ncol = 7)
    colnames(out) <- c("alpha2", "upsilon", "psi", "sig", "p0", "d0", "sig_mm")
  }
  
  
  
}else if(select_ups == "big ups"){
  
  if(share_sig == TRUE){
    
    p <- c(1,                         # alpha2
           log(2.5*1),                # ups
           qlogis(0.9),               # psi
           log(2.5*1),                # sig
           qlogis(0.1),               # p0
           log(100/ncell(scr_ss[[1]])) # d0
    )
    
    out <- matrix(NA, nrow = sims, ncol = 6)
    colnames(out) <- c("alpha2", "upsilon", "psi", "sig", "p0", "d0")
    
  }else if(share_sig == FALSE){
    
    p <- c(1,                           # alpha2
           log(2.5*1),                  # ups
           qlogis(0.9),                 # psi
           log(2.5*1),                  # sig
           qlogis(0.1),                 # p0
           log(100/ncell(scr_ss[[1]])), # d0
           log(2.5*1)                   # sig_mm
    )
    
    out <- matrix(NA, nrow = sims, ncol = 7)
    colnames(out) <- c("alpha2", "upsilon", "psi", "sig", "p0", "d0", "sig_mm")
  }
  
}




#----Fit movement model----

source("./R/likelihoods/telem-cost-nll-becky2.R")

# Parallel setup
# ncores = detectCores() # Number of available cores -1 to leave for computer
# cl = makeCluster(ncores) # Make the cluster with that many cores
# registerDoParallel(cl)  

sim <- 1
file <- paste0(getwd(), "/output/", select_ups, "/model_data/cost_data_light/cost_data_", sim, ".RData")
spatdata <- readRDS(file)

teldata1   = teldata[[sim]][inds]
spatdata1  = spatdata[inds]
landscape1 = landscape[[sim]]
scr_ss1 = scr_ss[[sim]]
K1 = K 
scr_y1 = y[[sim]]
rm(teldata, spatdata, landscape, scr_ss, K, y)



#----Add an empty buffer around the state-space----
# This is removed by likelihood via Dan's edit
# and should allow for irregular rasters
scr_ss2 <- raster::extend(x = scr_ss1, y = c(1,1))

sqrt(ncell(scr_ss1))
sqrt(ncell(scr_ss2))

plot(scr_ss1)
plot(scr_ss2)
lines(extent(scr_ss2))

nrow(coordinates(scr_ss1))
nrow(coordinates(scr_ss2)[which(!is.na(values(scr_ss2))),])



# NLM likelihood evaluation
mmscreco <- nlm(
  scr_move_cost_like_rw,
  p,
  mod = "gauss", share_sigma = share_sig,
  hessian = F, print.level = 2,
  teldata   = teldata1,
  spatdata  = spatdata1,
  landscape = landscape1,
  scr_ss = scr_ss2,
  K = K1, scr_y = scr_y1, trap_locs = traps,
  dist = "lcp", popcost=T, popmove=T, fixcost=F, use.sbar=T, prj=NULL)

est <- mmscreco$estimate

# file <- paste0("./output/", select_ups, "/", file_id, "/mmscreco_", sim, ".txt")
# write.table(est, file)

# Back-transform point estimates
final <- c()
final[1] <- est[1]
final[2] <- exp(est[2])
final[3] <- plogis(est[3])
final[4] <- exp(est[4])
final[5] <- plogis(est[5])
final[6] <- exp(est[6])

if(share_sig == FALSE){final[7] <- exp(est[7])}

# Output
results[[sim]] <- final


# results <- foreach(sim=1:sims, .packages = c(.packages())) %dopar% {
#   
#   file <- paste0(getwd(), "/output/", select_ups, "/model_data/cost_data_light/cost_data_", sim, ".RData")
#   spatdata <- readRDS(file)
#   
#   # NLM likelihood evaluation
#   mmscreco <- nlm(
#     scr_move_cost_like,
#     p,
#     mod = "gauss", share_sigma = share_sig,
#     hessian = F, #print.level = 2,
#     teldata   = teldata[[sim]][inds],
#     spatdata  = spatdata[inds],
#     landscape = landscape[[sim]],
#     scr_ss = scr_ss[[sim]],
#     K = K, scr_y = y[[sim]], trap_locs = traps,
#     dist = "lcp", popcost=T, popmove=T, fixcost=F, use.sbar=T, prj=NULL)
#   
#   est <- mmscreco$estimate
#   
#   file <- paste0("./output/", select_ups, "/", file_id, "/mmscreco_", sim, ".txt")
#   write.table(est, file)
#   
#   # Back-transform point estimates
#   final <- c()
#   final[1] <- est[1]
#   final[2] <- exp(est[2])
#   final[3] <- plogis(est[3])
#   final[4] <- exp(est[4])
#   final[5] <- plogis(est[5])
#   final[6] <- exp(est[6])
#   
#   if(share_sig == FALSE){final[7] <- exp(est[7])}
#   
#   # Output
#   final
#   
# }
# stopCluster(cl)
# tf <- Sys.time()
# t_total <- tf-t0

# FINAL OUTPUT AS MATRIX
out[1:ncell(out)] <- matrix(unlist(results), nrow=sims, byrow=T)




