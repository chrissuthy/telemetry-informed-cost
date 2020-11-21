library(gdistance)
library(dplyr)


#----Load data----

teldata   <- readRDS("output/model_data/teldata_raw.RData") # colnames? # fine
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

inds <- 1:5


#----Fit movement model----

source("R/main/models/scr_move_cost_like.R")

# Number of sims
#sims <- length(y)
sims_start <- 1
sims_end <- 10
sims_all <- seq(sims_start, sims_end, by = 1)
nsims <- length(sims_all)

out <- matrix(NA, nrow = nsims, ncol = 6)
colnames(out) <- c("alpha2", "upsilon", "psi", "sig", "p0", "d0")
results <- list()

# # Parallel setup
# ncores = detectCores() # Number of available cores -1 to leave for computer
# cl = makeCluster(ncores) # Make the cluster with that many cores
# registerDoParallel(cl)  

# Cluster!
t0 <- Sys.time()
#results <- foreach(sim=1:sims, .packages = c(.packages())) %dopar% {
for(i in 1:length(sims_all)){
  
  sim <- sims_all[i]
  
  file <- paste0(getwd(), "/output/model_data/cost_data_light/cost_data_", sim, ".RData")
  spatdata <- readRDS(file)
  
  set.seed(sim)
  
  
  # NLM likelihood evaluation
  mmscreco <- nlm(
    scr_move_cost_like,
    c(1,                         # alpha2
      log(0.625),                # ups
      qlogis(0.9),               # psi
      log(2.5),                  # sig
      qlogis(0.1),               # p0
      log(100/ncell(scr_ss[[1]])) # d0
    ),
    mod = "gauss",
    hessian = F, #print.level = 2,
    teldata   = teldata[[sim]][inds],
    spatdata  = spatdata[inds],
    landscape = landscape[[sim]],
    scr_ss = scr_ss[[sim]],
    K = K, scr_y = y[[sim]], trap_locs = traps,
    dist = "lcp", popcost=T, popmove=T, fixcost=F, use.sbar=T, prj=NULL)

  est <- mmscreco$estimate
  
  file <- paste0("output/mm_out/mmscreco_", sim, ".txt")
  write.table(est, file)
  
  # Back-transform point estimates
  final <- c()
  final[1] <- est[1]
  final[2] <- exp(est[2])
  final[3] <- plogis(est[3])
  final[4] <- exp(est[4])
  final[5] <- plogis(est[5])
  final[6] <- exp(est[6])
  
  results[[i]] <- final
  
  # Output
  # final
  
}
# stopCluster(cl)
tf <- Sys.time()
t_total <- tf-t0

# FINAL OUTPUT AS MATRIX
out[1:ncell(out)] <- matrix(unlist(results), nrow=nsims, byrow=T)

# SAVE IT ALL
file <- paste0("output/mm_out/results_", sims_start, "_", sims_end, ".RData")
saveRDS(results, file)
file <- paste0("output/mm_out/mmscreco_results_", sims_start, "_", sims_end, ".txt")
write.table(out, file)


