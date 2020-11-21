library(gdistance)
library(dplyr)
select = dplyr::select

#----Load data----

teldata   <- readRDS("output/model_data/teldata_raw.RData") # colnames? # fine
#spatdata_old  <- readRDS("output/model_data/cost_data.RData") # colnames? # fine
landscape <- readRDS("output/model_data/landscape.RData") # colnames?
y         <- readRDS("output/model_data/y.RData")
traps     <- readRDS("output/model_data/traps.RData") %>% as.matrix()
colnames(traps) <- c("X", "Y")
ss        <- readRDS("output/model_data/ss.RData") # colnames?
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

#inds <- sample(1:length(teldata[[1]]), size = 8)


#----Fit movement model----

source("R/models/scr_cost_like.R")

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
  
  file <- paste0("output/noMove_out/screco_", sim, ".txt")
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
saveRDS(results, "output/noMove_out/results.RData")
file <- paste0("output/noMove_out/screco_results.txt")
write.table(out, file)