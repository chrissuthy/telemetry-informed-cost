#remotes::install_github("ropensci/FedData")
library(FedData)
library(magrittr)
library(raster)
library(sf)
library(oSCR)
library(dplyr)
library(gdistance)
data("nybears")



#----Get cost data----

colnames(nybears$ss) <- c("x","y")

# Extent polygon
pgon <- polygon_from_extent(extent(nybears$ss)*1.5,
                            proj4string = "+proj=utm +datum=NAD83 +zone=18")

# NLCD data
NLCD <- get_nlcd(pgon, label = "nybear_nlcd", year = 2016, force.redo = T)
NLCD <- projectRaster(NLCD, crs=crs(pgon), res=30)

# NED data
# NED <- get_ned(pgon, label = "nybear_ned", force.redo = T)
# NED <- projectRaster(NED, crs=crs(pgon),res=30)

# Separate covariates
forest0 <- NLCD %in% 41:43
# elev <- NED

# Shrink axes
forest_df <- as.data.frame(forest0, xy=T) %>%
  mutate(x = x/1000, y = y/1000)
forest1 <- rasterFromXYZ(forest_df)

# Setup
bears.ss_ext <- polygon_from_extent(
  extent(nybears$ss/1000), 
  proj4string = "+proj=utm +datum=NAD83 +zone=18")

# Make ss at 1 km resolution
cost.xy <- bears.ss_ext %>%
  spsample(., cellsize=0.5, offset = c(0.5, 0.5),
           type = "regular") %>%
  as.data.frame() %>%
  select(x = x1, y = x2)
cost.vals <- raster::extract(x = forest1, y = cost.xy, buffer = 0.25, fun = mean)
forest <- rasterFromXYZ(cbind(cost.xy, cost.vals))


# Make ss at 1 km resolution
ss <- bears.ss_ext %>%
  spsample(., cellsize=1, offset = c(0.5, 0.5),
           type = "regular") %>%
  as.data.frame() %>%
  select(x = x1, y = x2) %>%
  cbind(., layer = 1) %>%
  rasterFromXYZ()

plot(forest)
points(ss)


# par(mfrow=c(1,2))
hist(values(forest), breaks = 50, main = "Forest")
# hist(values(elev), breaks = 50, main = "Elevation")
# par(mfrow=c(1,1))

# Plot the telemetry data
# par(mfrow=c(1,2))
plot(forest); points(nybears$teldata[,c("X_UTM","Y_UTM")]/1000, pch=16,cex=0.2)
# plot(elev); points(nybears$teldata[,c("X_UTM","Y_UTM")], pch=16,cex=0.2)
# par(mfrow=c(1,1))

#----SCR data----

traps <- nybears$traplocs/1000
colnames(traps) <- c("X","Y")


plot(forest)
lines(bears.ss_ext)
points(as.data.frame(ss, xy=T)[,1:2], col = "black", pch = 1)
points(traps, pch=20, col = "red", cex=2)
lines(nybears$teldata[,c("X_UTM","Y_UTM")]/1000, col = "blue")

# Organize raw teldata
df0 <- nybears$teldata %>%
  mutate(X_UTM = X_UTM/1000,
         Y_UTM = Y_UTM/1000) %>%
  select(id = animalid, fix = fixnum, 
         X_UTM, Y_UTM)

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
bears_teldata <- list(CU803)

plot(forest)
lines(CU803, col = 6, lwd=2)
lines(CU818, col = 2, lwd=2)
lines(CU905, col = 4, lwd=2)

# SPATDATA

spatdata_base <- as.data.frame(forest, xy = T)[,1:2]


# BEAR 1
CU803_sbar <- colMeans(CU803) %>%
  as.numeric() %>%
  matrix(ncol = 2)

# Make the extent
CU803_ext <- extent(
  CU803_sbar[,1] - 4.3*4,
  CU803_sbar[,1] + 4.3*4,
  CU803_sbar[,2] - 4.3*4,
  CU803_sbar[,2] + 4.3*4
)

# Get raster from that
CU803_r <- crop(forest, CU803_ext)
CU803_r_df <- as.data.frame(CU803_r, xy=T)[,1:2]

# Find sbar on raster
sbar_indx <- raster::extract(
  x = CU803_r, 
  y = CU803_sbar, 
  cellnumbers=T)[,1]
sbar_on_r <- CU803_r_df[sbar_indx,c("x", "y")]

# Record sbar pixel as 1 in sbar =0 col
tmp_result <- CU803_r_df %>%
  select(x,y) %>%
  mutate(sbar = ifelse(
    (x == sbar_on_r[,1]) & (y == sbar_on_r[,2]),
    1,0)) %>%
  as.matrix()

# Assign
spatdata_CU803 <- tmp_result

par(mfrow=c(1,3))
plot(forest)
lines(CU803_ext)
points(spatdata_CU803[,1:2])
lines(CU803, col = 6, lwd = 2)
points(CU803_sbar, pch = 20, col = "yellow", cex = 4)


# BEAR 2
CU818_sbar <- colMeans(CU818) %>%
  as.numeric() %>%
  matrix(ncol = 2)

# Make the extent
CU818_ext <- extent(
  CU818_sbar[,1] - 4.300*4,
  CU818_sbar[,1] + 4.300*4,
  CU818_sbar[,2] - 4.300*4,
  CU818_sbar[,2] + 4.300*4
)

# Get raster from that
CU818_r <- crop(forest, CU818_ext)
CU818_r_df <- as.data.frame(CU818_r, xy=T)[,1:2]

# Find sbar on raster
sbar_indx <- raster::extract(
  x = CU818_r, 
  y = CU818_sbar, 
  cellnumbers=T)[,1]
sbar_on_r <- CU818_r_df[sbar_indx,c("x", "y")]

# Record sbar pixel as 1 in sbar =0 col
tmp_result <- CU818_r_df %>%
  select(x,y) %>%
  mutate(sbar = ifelse(
    (x == sbar_on_r[,1]) & (y == sbar_on_r[,2]),
    1,0)) %>%
  as.matrix()

# Assign
spatdata_CU818 <- tmp_result

plot(forest)
lines(CU818_ext)
points(spatdata_CU818[,1:2])
lines(CU818, col = 2, lwd = 2)
points(CU818_sbar, pch = 20, col = "yellow", cex = 4)


# BEAR 3
CU905_sbar <- colMeans(CU905) %>%
  as.numeric() %>%
  matrix(ncol = 2)


# Make the extent
CU905_ext <- extent(
  CU905_sbar[,1] - 4.300*4,
  CU905_sbar[,1] + 4.300*4,
  CU905_sbar[,2] - 4.300*4,
  CU905_sbar[,2] + 4.300*4
)

# Get raster from that
CU905_r <- crop(forest, CU905_ext)
CU905_r_df <- as.data.frame(CU905_r, xy=T)[,1:2]

# Find sbar on raster
sbar_indx <- raster::extract(
  x = CU905_r, 
  y = CU905_sbar, 
  cellnumbers=T)[,1]
sbar_on_r <- CU905_r_df[sbar_indx,c("x", "y")]

# Record sbar pixel as 1 in sbar =0 col
tmp_result <- CU905_r_df %>%
  select(x,y) %>%
  mutate(sbar = ifelse(
    (x == sbar_on_r[,1]) & (y == sbar_on_r[,2]),
    1,0)) %>%
  as.matrix()

# Assign
spatdata_CU905 <- tmp_result

plot(forest)
lines(CU905_ext)
points(spatdata_CU905[,1:2])
lines(CU905, col = 4, lwd = 2)
points(CU905_sbar, pch = 20, col = "yellow", cex = 4)

par(mfrow=c(1,1))


# FINAL SPATDATA
spatdata <- list(spatdata_CU803, spatdata_CU818, spatdata_CU905)
spatdata <- list(spatdata_CU803)

rm(NLCD, pgon, df0, CU803, CU818, CU905, spatdata_base, sbar_indx,
   tmp_result, spatdata_CU803, spatdata_CU818, spatdata_CU905,
   CU803_r, CU803_r_df, CU818_r, CU818_r_df, CU905_r, CU905_r_df)


#----SCR ecological distance----

source("./R/likelihoods/scr_move_cost_like_SigmaFlag.R")

# Starting values
p <- c(1,                          # alpha2
       log(0.6),                   # ups
       qlogis(0.5),                # psi
       log(4.3),                   # sig
       qlogis(0.05),               # p0
       log(100/ncell(ss)),         # d0
       log(4.3)                    # sig_mm
)


t0 <- Sys.time()
# NLM likelihood evaluation
mm_forest <- nlm(
  scr_move_cost_like,
  p = p,
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
tf <- Sys.time()
tf-t0

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

ASE <- mm_forest$hessian %>% solve %>% diag %>% sqrt
out <- rbind(est, ASE) %>% as.data.frame() %>% mutate(model = "~forest") %>% select(model, everything())
out <- out %>% mutate_if(is.numeric, round, digits = 4)
rownames(out) <- c("MLE", "SE")
colnames(out) <- c("model", "cost", "sigma", "psi", "sig_det", "p0", "d0", "sig_ind")
out

write.table(out, "./output/nybears/ntel=1w1_share=F_MSE_SE_forest.txt")
save.image("./output/nybears/ntel=1w1_share=F_forest.RData")
