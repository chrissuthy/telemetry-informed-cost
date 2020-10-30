
"PART 2: FIT MODEL TO SIMULATIONS"

#----LOAD DATA----

teldata   <- readRDS("output/model_data/teldata_raw.RData") # colnames?
spatdata  <- readRDS("output/model_data/cost_data.RData") # colnames?
landscape <- readRDS("output/model_data/landscape.RData") # colnames?
y         <- readRDS("output/model_data/y.RData")
traps     <- readRDS("output/model_data/traps.RData") # X Y
ss        <- readRDS("output/model_data/ss.RData") # colnames?
K         <- 90


#----Model fitting----

# Number of model fits
nfits <- 10

# Parallel setup
ncores = detectCores()-1 # Number of available cores -1 to leave for computer
cl = makeCluster(ncores) # Make the cluster with that many cores
registerDoParallel(cl)  

# Data-collection list
parallel_out <- list()

# Data-collection matrix
out <- matrix(NA, nrow = nfits, ncol = 6)
colnames(out) <- c("alpha2", "sigma", "p0", "d0", "upsilon", "psi")

# Cluster!
t0 <- Sys.time()
results <- foreach(sim=1:nfits, .packages = c(.packages())) %dopar% {
  
  "SCR+TELEM"
  
  # NLM likelihood evaluation
  mm <- nlm(scr_move_cost_like, mod = "gauss",
            c(0, 0, 0, 0, 0, 0), hessian = T,
            teldata   = teldata[[sim]], 
            spatdata  = spatdata[[sim]],
            landscape = landscape[[sim]],
            K = K, scr_y = y[[sim]], trap_locs = traps,
            dist = "lcp", popcost=T, popmove=T, fixcost=F, use.sbar=T, prj=NULL)
  
  # Back-transform estimates
  est <- mm$estimate
  
  # Back-transform estimates
  final <- rep(NA, ncol(out))
  final[1] <- est[1]          # alpha2
  final[2] <- exp(est[4])     # sigma
  final[3] <- plogis(est[5])  # p0
  final[4] <- exp(est[6])     # d0
  final[5] <- exp(est[2])     # upsilon
  final[6] <- plogis(est[3])  # psi
  
  file <- paste0("simout/mmscreco_", sim, ".txt")
  if(TRUE) write.table(final, file)

  # Time est
  est_t_total <- ((Sys.time()-t0)/sim)*nfits
  end_time <- t0 + est_t_total
  file <- paste0("simout/stopwatch.txt")
  if(TRUE) write.table(end_time, file)
  
  # Output
  final
}
stopCluster(cl)
tf <- Sys.time()
t_total <- tf-t0


"RESULTS"

# Combine
out[1:ncell(out)] <- matrix(unlist(results), nrow=nfits, byrow=T)

