library(gdistance)
library(dplyr)


#----Load data----

teldata   <- readRDS("output/model_data/teldata_raw.RData") # colnames? # fine
spatdata_old  <- readRDS("output/model_data/cost_data.RData") # colnames? # fine
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


# par(mfrow=c(2,4))
# for(i in 1:8){
#   plot((spatdata[[1]][inds][[i]])[,1:2], 
#        pch = 16, col = "gray80", cex = 0.5, asp = 1)
#   lines((teldata[[1]][inds][[i]]))
# }
# par(mfrow=c(1,1))



#----1 sim, first position----

sim <- 1
inds <- sample(1:length(teldata[[1]]), size = 8)


#----Fit movement model----

source("R/main/models/scr_move_cost_like.R")

t1 <- Sys.time()

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
  hessian = T, print.level = 2,
  teldata   = teldata[[sim]][inds], 
  spatdata  = spatdata[[sim]][inds],
  landscape = landscape[[sim]],
  scr_ss = scr_ss[[sim]],
  K = K, scr_y = y[[sim]], trap_locs = traps,
  dist = "lcp", popcost=T, popmove=T, fixcost=F, use.sbar=T, prj=NULL)

t2 <- Sys.time()
t_mmscreco <- t2-t1


beepr::beep()

#----Fit model w/o movement----

source("R/main/models/scr_cost_like.R")

t3 <- Sys.time()

# NLM likelihood evaluation
screco <- nlm(
  scr_cost_like, mod = "gauss",
  c(2,                         # alpha2 
    log(4),                    # sigma 
    qlogis(0.1),               # p0 
    log(50/ncell(scr_ss[[1]])) # d0 
    ), 
  hessian = T,
  landscape = landscape[[sim]],
  scr_ss = scr_ss[[sim]],
  K = K, scr_y = y[[sim]], trap_locs = traps,
  dist = "lcp")

t4 <- Sys.time()
t_screco <- t4-t3

