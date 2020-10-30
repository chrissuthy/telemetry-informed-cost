"PART 2: FIT MODEL TO SIMULATIONS"

# Number of mdoel fits
nfits <- nsims*2 # fit all with scr, then with movement
model = rep(c("scr-telem", "scr+telem"), each = nsims)

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
  
  # Account for doublee-loop through
  sim_c = ifelse(sim>nsims, sim-nsims, sim)
  
  "MODELS"
  
  if(model[sim] == "scr-telem"){
    
    "SCR-TELEM"
    
    # NLM likelihood evaluation
    mm <- nlm(scr_cost_like, mod = "gauss",
              c(0, 0, 0, 0), hessian = T,
              landscape = landscape_ALL[[sim_c]],
              K = K, scr_y = Y_ALL[[sim_c]], trap_locs = traplocs,
              dist = "lcp", prj=NULL)
    
    # Back-transform estimates
    est <- mm$estimate
    
    # Back transform estimates
    final <- rep(NA, ncol(out))
    final[1] <- est[1]          # alpha2
    final[2] <- exp(est[2])     # sigma
    final[3] <- plogis(est[3])  # p0
    final[4] <- exp(est[4])     # d0
    
    file <- paste0("C:/Users/dupon/Desktop/sims_10052020/scr_sim", sim, ".txt")
    if(TRUE) write.table(final, file)
    
  }else if(model[sim] == "scr+telem"){
    
    "SCR+TELEM"
    
    # NLM likelihood evaluation
    mm <- nlm(scr_move_cost_like, mod = "gauss",
              c(0, 0, 0, 0, 0, 0), hessian = T,
              teldata = teldata_raw_ALL[[sim_c]], 
              spatdata = cost.data_ALL[[sim_c]],
              landscape = landscape_ALL[[sim_c]],
              K = K, scr_y = Y_ALL[[sim_c]], trap_locs = traplocs,
              dist = "lcp", popcost=T, popmove=T, fixcost=F, use.sbar=T, prj=NULL)
    
    # Back-transform estimates
    est <- mm$estimate
    
    # Back transform estimates
    final <- rep(NA, ncol(out))
    final[1] <- est[1]          # alpha2
    final[2] <- exp(est[4])     # sigma
    final[3] <- plogis(est[5])  # p0
    final[4] <- exp(est[6])     # d0
    final[5] <- exp(est[2])     # upsilon
    final[6] <- plogis(est[3])  # psi
    
    file <- paste0("C:/Users/dupon/Desktop/sims_10052020/scr_move_sim", sim, ".txt")
    if(TRUE) write.table(final, file)
    
  }
  
  # Time est
  est_t_total <- ((Sys.time()-t0)/sim)*nfits
  end_time <- t0 + est_t_total
  file <- paste0("C:/Users/dupon/Desktop/sims_10052020/stopwatch.txt")
  if(TRUE) write.table(end_time, file)
  
  # Output
  final
}
stopCluster(cl)
tf <- Sys.time()
t_total <- tf-t0


"RESULTS"

pal <- wesanderson::wes_palette("Zissou1", 5, type = "continuous")

# Combine
out[1:ncell(out)] <- matrix(unlist(results), nrow=nfits, byrow=T)

