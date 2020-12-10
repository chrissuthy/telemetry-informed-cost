#remotes::install_github("ropensci/FedData")
library(FedData)
library(magrittr)
library(raster)
library(sf)
library(oSCR)
library(dplyr)
library(gdistance)
library(spData)
data("nybears")


#----Get state data----

ny <- map_data("state", region="new york")[,c("long","lat")]
colnames(ny) <- c("x","y")
ny <- SpatialPoints(ny, proj4string =  CRS("+proj=longlat +datum=WGS84 +no_defs "))
ny <- spTransform(ny, "+proj=utm +datum=NAD83 +zone=18")


#----Get cost data----

colnames(nybears$ss) <- c("x","y")

# Extent polygon
pgon <- polygon_from_extent(ny, proj4string = "+proj=utm +datum=NAD83 +zone=18")

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

# Make ss at 1 km resolution
cost.xy <- polygon_from_extent(
  extent(forest1), 
  proj4string = "+proj=utm +datum=NAD83 +zone=18") %>%
  spsample(., cellsize=1, offset = c(0.5, 0.5),
           type = "regular") %>%
  as.data.frame() %>%
  select(x = x1, y = x2)
cost.vals <- raster::extract(x = forest1, y = cost.xy, buffer = 0.5, fun = mean)
forest <- rasterFromXYZ(cbind(cost.xy, cost.vals))


# Make ss at 1 km resolution
ss <- pgon %>%
  spsample(., cellsize=1, offset = c(0.5, 0.5),
           type = "regular") %>%
  as.data.frame() %>%
  select(x = x1, y = x2) %>%
  cbind(., layer = 1) %>%
  rasterFromXYZ()

plot(forest)
points(as.data.frame(ss,xy=T)[,1:2])