ss <- expand.grid(x=2:10, y = 2:10)
ss$z <- rnorm(nrow(ss), mean = 5, sd = 1)
landscape00 <- rasterFromXYZ(ss)
landscape01 <- extend(landscape00, y = c(1,1))
ss0 <- as.data.frame(landscape01, xy=T)[,1:2] %>% as.matrix()

cost <- exp(1.5*landscape01) # Cost surface w/ proposed parameter
tr1 <- transition(cost,transitionFunction=function(x) 1/mean(x),directions=16)
tr1CorrC <- geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
D <- costDistance(tr1CorrC,ss0,ss0) # Cost distance

D
View(D)
