# # # # # # # # # # # # # # #
#                           #
#  ____HEY!!_DO THIS!!____  #
#     Parameterize sims     #
#                           #
# # # # # # # # # # # # # # #

# Right here:
select_ups  <- c("small ups", "big ups")[NULL]

# # # # # # # # # # # # # # #
#                           #
#       GREAT, THANKS!      #
#                           #
# # # # # # # # # # # # # # #


library(gdistance)
library(dplyr)
select = dplyr::select

#----Directory setup----

# Create that directory
file_id <- paste0("est_noMove")
new_dir <- paste0("output/", select_ups, "/", file_id)
if(!dir.exists(new_dir)){dir.create(new_dir)}


#----Load data----

# Landscape
file      <- paste0("output/", select_ups, "/model_data/landscape.RData")
landscape <- readRDS(file)

# Y array
file      <- paste0("output/", select_ups, "/model_data/y.RData")
y         <- readRDS(file)

# Traps
file      <- paste0("output/", select_ups, "/model_data/traps.RData")
traps     <- readRDS(file) %>% as.matrix()
colnames(traps) <- c("X", "Y")

# SS
file      <- paste0("output/", select_ups, "/model_data/ss.RData")
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


#----Fit movement model----

source("R/likelihoods/scr_cost_like.R")

out <- matrix(NA, nrow = sim, ncol = 4)
colnames(out) <- c("alpha2", "sig", "p0", "d0")

# Number of sims
sims <- length(y)

# Parallel setup
ncores = detectCores()-1 # Number of available cores -1 to leave for computer
cl = makeCluster(ncores) # Make the cluster with that many cores
registerDoParallel(cl)  

# Cluster!
t0 <- Sys.time()
results <- foreach(sim=1:sims, .packages = c(.packages())) %dopar% {
  
  # NLM likelihood evaluation
  screco <- nlm(
    scr_cost_like, mod = "gauss",
    c(1,                         # alpha2
      log(2.5),                 # sigma
      qlogis(0.1),               # p0
      log(100/ncell(scr_ss[[1]])) # d0
    ),
    hessian = T,
    landscape = landscape[[sim]],
    scr_ss = scr_ss[[sim]],
    K = K, scr_y = y[[sim]], trap_locs = traps,
    dist = "lcp")

  est <- screco$estimate
  
  file <- paste0("output/", select_ups, "/", file_id, "/screco_", sim, ".txt")
  write.table(est, file)
  
  # Back-transform point estimates
  final <- c()
  final[1] <- est[1]
  final[2] <- exp(est[2])
  final[3] <- plogis(est[3])
  final[4] <- exp(est[4])

  # Output
  final
  
}
stopCluster(cl)
tf <- Sys.time()
t_total <- tf-t0

# FINAL OUTPUT AS MATRIX
out[1:ncell(out)] <- matrix(unlist(results), nrow=sims, byrow=T)

# SAVE IT ALL
file <- paste0("output/", select_ups, "/", file_id, "/results.RData")
saveRDS(results, file)
file <- paste0("output/", select_ups, "/", file_id, "/results.txt")
write.table(out, file)