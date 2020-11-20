library(gdistance)
library(dplyr)


#----Load data----

#spatdata_old  <- readRDS("output/model_data/cost_data.RData") # colnames? # fine
spatdata_old <- readRDS(file.choose())

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


# Re-construct spatdata
spatdata <- list()
for(sim in 1:length(spatdata_old)){
  
  spatdata[[sim]] <- list()
  for(ind in 1:length(spatdata_old[[sim]])){
    
    tmp_df <- spatdata_old[[sim]][[ind]] %>%
      as.data.frame()
    
    sbar <- tmp_df %>%
      select(x, y) %>%
      colMeans() %>%
      as.numeric() %>%
      matrix(ncol = 2)
    
    tmp_r <- raster::rasterFromXYZ(tmp_df)
    
    sbar_indx <- raster::extract(x = tmp_r, y = sbar, cellnumbers=T)[,1]
    sbar_on_r <- tmp_df[sbar_indx,c("x", "y")]
    
    tmp_result <- tmp_df %>%
      select(x,y) %>%
      mutate(sbar = ifelse(
        (x == sbar_on_r[,1]) & (y == sbar_on_r[,2]),
        1,0)) %>%
      as.matrix()
    
    spatdata[[sim]][[ind]] <- tmp_result
    
  }
  
}


#----Number of tagged individuals----

inds <- c(26, 41, 15, 16, 4)


#----Fit movement model----

source("R/main/models/scr_move_cost_like.R")

# Number of sims
sims <- length(y)

out <- matrix(NA, nrow = sims, ncol = 6)
colnames(out) <- c("alpha2", "upsilon", "psi", "sig", "p0", "d0")

# Parallel setup
ncores = detectCores() # Number of available cores -1 to leave for computer
cl = makeCluster(ncores) # Make the cluster with that many cores
registerDoParallel(cl)  

# Cluster!
t0 <- Sys.time()
results <- foreach(sim=1:sims, .packages = c(.packages())) %dopar% {

  # NLM likelihood evaluation
  mmscreco <- nlm(
    scr_move_cost_like,
    c(2,                         # alpha2
      log(1),                    # ups
      qlogis(0.9),               # psi
      log(4),                    # sig
      qlogis(0.1),               # p0
      log(50/ncell(scr_ss[[1]])) # d0
    ),
    mod = "gauss",
    hessian = T, #print.level = 2,
    teldata   = teldata[[sim]][inds],
    spatdata  = spatdata[[sim]][inds],
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
  
  # Output
  final
  
}
stopCluster(cl)
tf <- Sys.time()
t_total <- tf-t0

# FINAL OUTPUT AS MATRIX
out[1:ncell(out)] <- matrix(unlist(results), nrow=sims, byrow=T)

# SAVE IT ALL
saveRDS(results, "output/mm_out/results.RData")
file <- paste0("output/mm_out/mmscreco_results.txt")
write.table(out, file)


