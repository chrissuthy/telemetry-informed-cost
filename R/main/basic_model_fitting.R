library(oSCR)

# Traps
traps <- readRDS("output/oct29_traps.RData")
tdf <- data.frame(
  Trap_id = 1:nrow(traps),
  X = traps[,1],
  Y = traps[,2])

# State-space
S <- readRDS("output/oct29_ss.RData")

# Prediction data frames
pred.df.det <- data.frame(Session = factor(1))
pred.df.sig <- data.frame(Session = factor(1))
pred.df.dens <- data.frame(Session = factor(1))


for(i in 1:10){
  
  # Reduce to captured
  y <- readRDS("output/oct29_edf_10sim.RData")[i,,,]
  captured <- apply(y,1,sum)>0
  y <- y[captured,,]
  
  # Make edf
  edf <- data.frame(which(y>0, arr.ind = T))
  edf <- data.frame(
    session = rep(1,nrow(edf)),
    ind = edf$dim1,
    trap = edf$dim2,
    occ = edf$dim3)
  
  # Compile data
  data <- data2oscr(
    edf = edf,
    tdf = list(tdf),
    sess.col = 1,
    id.col = 2,
    occ.col = 4,
    trap.col = 3,
    K = 90,
    ntraps = nrow(tdf))
  
  # several objects returned, but we just want the scrFrame
  sf <- data$scrFrame
  
  # State-space
  df <- data.frame(X=S$x, Y=S$y, Z=S$layer)
  ss <- list(df)
  class(ss) <- "ssDF"
  
  # PLOT
  # plot(S[,1:2], col = "gray80", asp = 1, cex = 0.1, axes=F)
  # spiderplot(sf, add=T)
  
  # Fit model
  model_formula <- list(D~1, p0~1, sig~1)
  scr0 <- oSCR.fit(
    model_formula, # model formulation
    sf,            # the scrFrame object
    ss,            # the ssDF object
    trimS = sf$mmdm*3.05)            
  
  # Model output
  Den <- get.real(scr0, type = "dens", newdata = pred.df.dens, d.factor = 1)
  Abu <- get.real(scr0, type = "dens", newdata = pred.df.dens, d.factor = nrow(scr0$ssDF[[1]]))
  sig <- get.real(scr0, type = "sig",  newdata = pred.df.sig)
  det <- get.real(scr0, type = "det",  newdata = pred.df.det)
  
}






