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
# NLCD <- get_nlcd(pgon, label = "nybear_nlcd", year = 2016, force.redo = T)
# NLCD <- projectRaster(NLCD, crs=crs(pgon), res=30)

# NED data
NED <- get_ned(pgon, label = "nybear_ned", force.redo = T)
NED <- projectRaster(NED, crs=crs(pgon),res=30)

# Separate covariates
# forest0 <- NLCD %in% 41:43
elev0 <- NED/1000

# Shrink axes
elev_df <- as.data.frame(elev0, xy=T) %>%
  mutate(x = x/1000, y = y/1000)
elev1 <- rasterFromXYZ(elev_df)

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
cost.vals <- raster::extract(x = elev1, y = cost.xy, buffer = 0.25, fun = mean)
elev <- rasterFromXYZ(cbind(cost.xy, cost.vals))


# Make ss at 1 km resolution
ss <- bears.ss_ext %>%
  spsample(., cellsize=1, offset = c(0.5, 0.5),
           type = "regular") %>%
  as.data.frame() %>%
  select(x = x1, y = x2) %>%
  cbind(., layer = 1) %>%
  rasterFromXYZ()

plot(elev)
points(as.data.frame(ss, xy = T))


# par(mfrow=c(1,2))
hist(values(elev), breaks = 50, main = "elev")
# hist(values(elev), breaks = 50, main = "Elevation")
# par(mfrow=c(1,1))

# Plot the telemetry data
# par(mfrow=c(1,2))
plot(elev); points(nybears$teldata[,c("X_UTM","Y_UTM")]/1000, pch=16,cex=0.2)
# plot(elev); points(nybears$teldata[,c("X_UTM","Y_UTM")], pch=16,cex=0.2)
# par(mfrow=c(1,1))

#----SCR data----

traps <- nybears$traplocs/1000
colnames(traps) <- c("X","Y")


plot(elev)
lines(bears.ss_ext)
points(as.data.frame(ss, xy=T)[,1:2], col = "black", pch = 1)
points(traps, pch=20, col = "red", cex=2)
lines(nybears$teldata[,c("X_UTM","Y_UTM")]/1000, col = "blue", lwd = 2)

rm(NED, cost.xy, elev_df, elev0, elev1, pgon)


#----SCR ecological distance----

source("R/likelihoods/scr_cost_like.R")

# Starting values
p <- c(1,                          # alpha2
       log(4.3),                   # sig
       qlogis(0.05),               # p0
       log(100/ncell(ss))          # d0
)


t0 <- Sys.time()
# NLM likelihood evaluation
mm_elev <- nlm(
  scr_cost_like,
  p = p,
  mod = "gauss",
  hessian = T, print.level = 2,
  landscape = elev,
  scr_ss = ss,
  K = nybears$K, 
  scr_y = nybears$y2d, 
  trap_locs = traps,
  dist = "lcp")
tf <- Sys.time()
tf-t0

est <- mm_elev$estimate

# Back-transform point estimates
final <- c()
final[1] <- est[1]
final[2] <- exp(est[2])
final[3] <- plogis(est[3])
final[4] <- exp(est[4])
final

ASE <- mm_elev$hessian %>% solve %>% diag %>% sqrt
out <- rbind(est, ASE) %>% as.data.frame() %>% mutate(model = "~elev") %>% select(model, everything())
out <- out %>% mutate_if(is.numeric, round, digits = 4)
rownames(out) <- c("MLE", "SE")
colnames(out) <- c("model", "cost", "sig_det", "p0", "d0")
out

write.table(out, "./output/nybears/nomove_MSE_SE_elev.txt")

save.image("./output/nybears/nomove_elev.RData")
