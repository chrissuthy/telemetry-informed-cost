#remotes::install_github("ropensci/FedData")
library(FedData)
library(magrittr)
library(raster)
library(sf)
library(oSCR)
data("nybears")



#----Get cost data----

# Extent polygon
pgon <- polygon_from_extent(extent(nybears$ss)*1.05,
                            proj4string = "+proj=utm +datum=NAD83 +zone=18")

# NLCD data
NLCD <- get_nlcd(pgon, label = "nybear_nlcd", year = 2016, force.redo = T)
NLCD <- projectRaster(NLCD, crs=crs(pgon), res=30)

# NED data
NED <- get_ned(pgon, label = "nybear_ned", force.redo = T)
NED <- projectRaster(NED, crs=crs(pgon),res=30)

# Separate covariates
forest <- NLCD %in% 41:43
elev <- NED


# Separate covariates
forest <- aggregate(crop(forest,pgon),fac=10,fun="mean")
elev <- aggregate(crop(elev,pgon),fac=10,fun="mean")
elev <- elev/1000


par(mfrow=c(1,2))
hist(values(forest), breaks = 50, main = "Forest")
hist(values(elev), breaks = 50, main = "Elevation")
par(mfrow=c(1,1))


# Plot the telemetry data
par(mfrow=c(1,2))
plot(forest); points(nybears$teldata[,c("X_UTM","Y_UTM")], pch=16,cex=0.2)
plot(elev); points(nybears$teldata[,c("X_UTM","Y_UTM")], pch=16,cex=0.2)
par(mfrow=c(1,1))

#----SCR data----

traps <- nybears$traplocs
colnames(traps) <- c("X","Y")

ss <- aggregate(crop(forest,nybears$ss), fact = 7, fun="mean")

plot(forest)
lines(pgon)
points(as.data.frame(ss, xy=T)[,1:2], col = "black", pch = 1)
points(nybears$traplocs, pch=20, col = "red", cex=2)
lines(nybears$teldata[,c("X_UTM","Y_UTM")], col = "blue")


#----SCR ecological distance----

source("R/likelihoods/scr_cost_like.R")


# FOREST

t0 <- Sys.time()
nomove.forest <- nlm(
  scr_cost_like,
  c(1,                        # alpha2
    log(2000),                # sigma
    qlogis(0.05),             # p0
    log(100/ncell(ss))        # d0
  ),
  print.level = 2,
  hessian = T,
  mod = "gauss",
  landscape = forest,
  scr_ss = ss,
  K = nybears$K, 
  scr_y = nybears$y2d, 
  trap_locs = traps,
  dist = "lcp")
tf <- Sys.time()
tf-t0

est <- nomove.forest$estimate
final.forest <- c()
final.forest[1] <- est[1]
final.forest[2] <- exp(est[2])
final.forest[3] <- plogis(est[3])
final.forest[4] <- exp(est[4])
final.forest



# ELEVATION

t0 <- Sys.time()
nomove.elev <- nlm(
  scr_cost_like,
  c(1,                        # alpha2
    log(2000),                # sigma
    qlogis(0.05),             # p0
    log(100/ncell(ss))        # d0
  ),
  print.level = 2,
  hessian = T,
  mod = "gauss",
  landscape = elev,
  scr_ss = ss,
  K = nybears$K,
  scr_y = nybears$y2d,
  trap_locs = traps,
  dist = "lcp")
tf <- Sys.time()
tf-t0

est <- nomove.elev$estimate
final.elev <- c()
final.elev[1] <- est[1]
final.elev[2] <- exp(est[2])
final.elev[3] <- plogis(est[3])
final.elev[4] <- exp(est[4])
final.elev


#----Compile results----

final <- rbind(final.forest, final.elev) %>%
  as.data.frame() %>%
  mutate(covariate = c("forest", "elevation")) %>%
  select(covariate, everything()) %>%
  mutate(V4 = V4 * ncell(ss))

colnames(final) <- c("covariate", "cost",  "sigma", "p0", "N")

final

# save.image("./output/nybears/nomove_ForestElevation.RData")
