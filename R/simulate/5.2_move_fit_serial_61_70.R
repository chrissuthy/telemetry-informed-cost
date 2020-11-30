# # # # # # # # # # # # # # #
#                           #
#  ____HEY!!_DO THIS!!____  #
#     Parameterize sims     #
#                           #
# # # # # # # # # # # # # # #

# Right here:
select_ups  <- c("small ups", "big ups")[2]
select_ntel <- c(1, 3, 5)[3]
share_sig   <- c(TRUE, FALSE)[2]

# Break it up:
sims_start <- 61
sims_end <- 70
sims_all <- seq(sims_start, sims_end, by = 1)
nsims <- length(sims_all)

# WRITING TO DROPBOX? SHOULD FIT IN THE LAST LINES OF THIS SCRIPT.

# # # # # # # # # # # # # # #
#                           #
#       GREAT, THANKS!      #
#                           #
# # # # # # # # # # # # # # #


#----Get started----

library(gdistance)
library(dplyr)


#----Directory setup----

# Create that directory
file_id <- paste0("est_ntel=", select_ntel, "_", "share=", share_sig)
new_dir <- paste0("./output/", select_ups, "/", file_id)
if(!dir.exists(new_dir)){dir.create(new_dir)}


#----Load data----

# Telemetry data
file1     <- paste0("./output/", select_ups, "/model_data/teldata_raw.RData")
teldata   <- readRDS(file1)

# Landscape
file2     <- paste0("./output/", select_ups, "/model_data/landscape.RData")
landscape <- readRDS(file2)

# Y array
file3     <- paste0("./output/", select_ups, "/model_data/y.RData")
y         <- readRDS(file3)

# Traps
file4     <- paste0("./output/", select_ups, "/model_data/traps.RData")
traps     <- readRDS(file4) %>% as.matrix()
colnames(traps) <- c("X", "Y")

# SS
file5     <- paste0("./output/", select_ups, "/model_data/ss.RData")
ss        <- readRDS(file5) # colnames?

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
    
    out <- matrix(NA, nrow = nsims, ncol = 6)
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
    
    out <- matrix(NA, nrow = nsims, ncol = 7)
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
    
    out <- matrix(NA, nrow = nsims, ncol = 6)
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
    
    out <- matrix(NA, nrow = nsims, ncol = 7)
    colnames(out) <- c("alpha2", "upsilon", "psi", "sig", "p0", "d0", "sig_mm")
  }
  
}




#----Fit movement model----

source("./R/likelihoods/scr_move_cost_like_SigmaFlag.R")

results <- list()

t0 <- Sys.time()
for(i in 1:length(sims_all)){
  
  sim <- sims_all[i]
  
  set.seed(sim)
  
  file <- paste0(getwd(), "/output/", select_ups, "/model_data/cost_data_light/cost_data_", sim, ".RData")
  spatdata <- readRDS(file)
  
  # NLM likelihood evaluation
  mmscreco <- nlm(
    scr_move_cost_like,
    p,
    mod = "gauss", share_sigma = share_sig,
    hessian = F, #print.level = 2,
    teldata   = teldata[[sim]][inds],
    spatdata  = spatdata[inds],
    landscape = landscape[[sim]],
    scr_ss = scr_ss[[sim]],
    K = K, scr_y = y[[sim]], trap_locs = traps,
    dist = "lcp", popcost=T, popmove=T, fixcost=F, use.sbar=T, prj=NULL)
  
  est <- mmscreco$estimate
  
  file <- paste0("./output/", select_ups, "/", file_id, "/mmscreco_", sim, ".txt")
  write.table(est, file)
  
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
  results[[i]] <- final
  
}
tf <- Sys.time()
t_total <- tf-t0

# FINAL OUTPUT AS MATRIX
out[1:ncell(out)] <- matrix(unlist(results), nrow=nsims, byrow=T)

# SAVE IT ALL
file <- paste0("./output/", select_ups, "/", file_id, "/results_", sims_start, "_", sims_end, ".RData")
saveRDS(results, file)
file <- paste0("./output/", select_ups, "/", file_id, "/results_", sims_start, "_", sims_end, ".txt")
write.table(out, file)

