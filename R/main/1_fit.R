library(gdistance)


#----Load data----

teldata   <- readRDS("output/model_data/teldata_raw.RData") # colnames? # fine
spatdata  <- readRDS("output/model_data/cost_data.RData") # colnames? # fine
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
  scr_ss[[i]] <- crop(landscape[[i]], extent(ss))
}


#----1 sim, first position----

sim = 1


#----Fit movement model----

source("R/main/models/scr_move_cost_like.R")

t1 <- Sys.time()

# NLM likelihood evaluation
mmscreco <- nlm(
  scr_move_cost_like, mod = "gauss",
  c(3,            # alpha2 
    log(4),       # sigma 
    qlogis(0.2),    # p0 
    log(70/10000),  # d0 
    log(1),       # upsilon   
    qlogis(0.8)),  # psi   
  hessian = T,
  teldata   = teldata[[sim]], 
  spatdata  = spatdata[[sim]],
  landscape = landscape[[sim]],
  scr_ss = scr_ss[[sim]],
  K = K, scr_y = y[[sim]], trap_locs = traps,
  dist = "lcp", popcost=T, popmove=T, fixcost=F, use.sbar=T, prj=NULL)

t2 <- Sys.time()
t_mmscreco <- t2-t1


#----Fit model w/o movement----

source("R/main/models/scr_cost_like.R")

t3 <- Sys.time()

# NLM likelihood evaluation
screco <- nlm(
  scr_cost_like, mod = "gauss",
  c(2.1,            # alpha2 
    log(4.1),       # sigma 
    qlogis(0.2),    # p0 
    log(51/10000)), # d0 
  hessian = T,
  landscape = landscape[[sim]],
  scr_ss = scr_ss[[sim]],
  K = K, scr_y = y[[sim]], trap_locs = traps,
  dist = "lcp")

t4 <- Sys.time()
t_screco <- t4-t3

