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

# STATE_SPACE
ss <- aggregate(crop(forest,nybears$ss), fact = 7, fun="mean")

plot(forest)
lines(pgon)
points(as.data.frame(ss, xy=T)[,1:2], col = "black", pch = 1)
points(nybears$traplocs, pch=20, col = "red", cex=2)
lines(nybears$teldata[,c("X_UTM","Y_UTM")], col = "blue")


# TELDATA

# Organize raw teldata
df0 <- nybears$teldata %>%
  select(id = animalid, fix = fixnum, X_UTM, Y_UTM)

# Convert continous space fixes to pixel centroids
fix_cells <- extract(forest, df0[,c("X_UTM","Y_UTM")], cellnumbers=T)[,1]
fix_cells_xy <- xyFromCell(forest, fix_cells)
colnames(fix_cells_xy) <- c("x","y")
df <- cbind(df0, fix_cells_xy)

CU803 <- df %>%
  filter(id == "CU803") %>%
  arrange(fix) %>%
  select(x, y)

CU818 <- df %>%
  filter(id == "CU818") %>%
  arrange(fix) %>%
  select(x, y)

CU905 <- df %>%
  filter(id == "CU905") %>%
  arrange(fix) %>%
  select(x, y)

bears_teldata <- list(CU803, CU818, CU905)


# SPATDATA

spatdata_base <- as.data.frame(forest, xy = T)[,1:2]


# BEAR 1
CU803_sbar <- colMeans(CU803) %>%
  as.numeric() %>%
  matrix(ncol = 2)

sbar_indx <- raster::extract(
  x = forest, 
  y = CU803_sbar, 
  cellnumbers=T)[,1]
sbar_on_r <- spatdata_base[sbar_indx,c("x", "y")]

tmp_result <- spatdata_base %>%
  select(x,y) %>%
  mutate(sbar = ifelse(
    (x == sbar_on_r[,1]) & (y == sbar_on_r[,2]),
    1,0)) %>%
  as.matrix()

spatdata_CU803 <- tmp_result


# BEAR 2
CU818_sbar <- colMeans(CU818) %>%
  as.numeric() %>%
  matrix(ncol = 2)

sbar_indx <- raster::extract(
  x = forest, 
  y = CU818_sbar, 
  cellnumbers=T)[,1]
sbar_on_r <- spatdata_base[sbar_indx,c("x", "y")]

tmp_result <- spatdata_base %>%
  select(x,y) %>%
  mutate(sbar = ifelse(
    (x == sbar_on_r[,1]) & (y == sbar_on_r[,2]),
    1,0)) %>%
  as.matrix()

spatdata_CU818 <- tmp_result


# BEAR 3
CU905_sbar <- colMeans(CU905) %>%
  as.numeric() %>%
  matrix(ncol = 2)

sbar_indx <- raster::extract(
  x = forest, 
  y = CU905_sbar, 
  cellnumbers=T)[,1]
sbar_on_r <- spatdata_base[sbar_indx,c("x", "y")]

tmp_result <- spatdata_base %>%
  select(x,y) %>%
  mutate(sbar = ifelse(
    (x == sbar_on_r[,1]) & (y == sbar_on_r[,2]),
    1,0)) %>%
  as.matrix()

spatdata_CU905 <- tmp_result

# FINAL SPATDATA
spatdata <- list(spatdata_CU803, spatdata_CU818, CU905_sbar)




#----SCR ecological distance----

source("./R/likelihoods/scr_move_cost_like_SigmaFlag.R")


# NLM likelihood evaluation
mm_forest <- nlm(
  scr_move_cost_like,
  p,
  mod = "gauss", share_sigma = F,
  hessian = T, print.level = 2,
  teldata   = bears_teldata,
  spatdata  = spatdata,
  landscape = forest,
  scr_ss = ss,
  K = nybears$K, 
  scr_y = nybears$y2d, 
  trap_locs = traps,
  dist = "lcp", 
  popcost=T, 
  popmove=T, 
  fixcost=F, 
  use.sbar=T, 
  prj=NULL)

est <- mm_forest$estimate

# Back-transform point estimates
final <- c()
final[1] <- est[1]
final[2] <- exp(est[2])
final[3] <- plogis(est[3])
final[4] <- exp(est[4])
final[5] <- plogis(est[5])
final[6] <- exp(est[6])
final[7] <- exp(est[7])
