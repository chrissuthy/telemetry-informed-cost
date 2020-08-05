library(oSCR)
library(gdistance)
library(NLMR)
library(landscapetools)
library(raster)
library(sf)
library(ggplot2)
library(viridis)

set.seed(202006)

telemetry_informed_cost.fn <- function(){}

# Create the landscape ----------------------------------------------------
autocorr <- 3
gauss_nlm <- nlm_gaussianfield(ncol = 100, nrow = 100, resolution = 1,
                               autocorr_range = autocorr)
use_this <- gauss_nlm
landscape <- rasterFromXYZ(cbind(coordinates(use_this)/10,
                                 values(use_this)))
plot(landscape)

# Create cost surface -----------------------------------------------------
cost_parameter <- 1.5
cost_surface <- exp(cost_parameter * landscape^2)
transistion_surface <- geoCorrection(
                        transition(cost_surface,
                                   transitionFunction = function(x) (1/(mean(x))),
                                   direction = 16),
                       scl = F)


# Create state space ------------------------------------------------------
res <- 0.125
statespace <- data.frame(
               expand.grid(
                X = seq(round(bbox(cost_surface)[1,1])+res/2,
                        round(bbox(cost_surface)[1,2])-res/2,
                        res),
                Y = seq(round(bbox(cost_surface)[2,1])+res/2,
                        round(bbox(cost_surface)[2,2])-res/2,
                        res)))



# SCR parameters ----------------------------------------------------------
abundance <- 300
area <- prod(round(apply(bbox(cost_surface),1,diff),2))

abs.density <- abundance / area
pix.density <- abundance / nrow(statespace)

sigma <- 0.35
p0 <- 0.1

K <- 5

acs <- cbind(X = round(runif(abundance,0,10),2),
             Y = round(runif(abundance,0,10),2))


# Creat traps -------------------------------------------------------------

traplocs <- as.matrix(
             expand.grid(X = seq(5-3*2*sigma,5+3*2*sigma,length=7),
                         Y = seq(5-3*2*sigma,5+3*2*sigma,length=7)))



# A vizual test -----------------------------------------------------------
plot(cost_surface, col=adjustcolor(viridis(100,option="A",direction=-1),0.6),axes=FALSE)
points(statespace, col="black", pch=16, cex=0.1)
points(traplocs, pch=15, col="blue")
points(acs, pch=16, col="red", cex = 0.2)

test_n <- 5
which_test_acs <- sample(1:nrow(acs),test_n)
test_ac <- st_as_sf(x = data.frame(id = 1:test_n,
                                   X = acs[which_test_acs,1],
                                   Y = acs[which_test_acs,2]),
                    coords = c("X","Y"))  
test_95hr <- st_buffer(test_ac,sqrt(5.99)*sigma, joinStyle = )
plot(st_geometry(test_95hr), add=TRUE, col = adjustcolor("black",0.2))
plot(st_geometry(test_ac), add=TRUE, bg = "red", pch=21, cex=1.5)


# Telemetry points --------------------------------------------------------

telemetry_n <- 16
which_telemetered_acs <- sample(1:nrow(acs), telemetry_n)
telemetered_acs <- acs[which_telemetered_acs,]
n_fixes <- 90*24
n_by_cell_freq <- matrix(NA,nrow=telemetry_n,ncell(cost_surface))

view_multiplier <- 5
plotit <- TRUE

par(mfrow=c(4,4), oma=c(0,0,0,0),mar=c(1,1,1,1))
for(i in 1:telemetry_n){
  tmp_from <- matrix(telemetered_acs[i,],1,2,byrow=TRUE)
  tmp_to <- coordinates(cost_surface)
  telem_rsf_dmat <- costDistance(transistion_surface,
                                 fromCoords = tmp_from,
                                 toCoords = tmp_to)
  telem_rsf_lammat <- exp(1-telem_rsf_dmat^2/(2*sigma^2))
  n_by_cell_freq[i,] <- c(rmultinom(1,n_fixes,telem_rsf_lammat/sum(telem_rsf_lammat)))
  
  if(plotit){
    tmp_rsf <- rasterFromXYZ(cbind(coordinates(cost_surface),
                                   c(telem_rsf_lammat)))
    ext <- as(st_polygon(list(rbind(c(tmp_from[1]-view_multiplier*sigma,tmp_from[2]-view_multiplier*sigma), 
                                    c(tmp_from[1]-view_multiplier*sigma,tmp_from[2]+view_multiplier*sigma), 
                                    c(tmp_from[1]+view_multiplier*sigma,tmp_from[2]+view_multiplier*sigma), 
                                    c(tmp_from[1]+view_multiplier*sigma,tmp_from[2]-view_multiplier*sigma),
                                    c(tmp_from[1]-view_multiplier*sigma,tmp_from[2]-view_multiplier*sigma)))),Class = "Spatial")
  
    tmp_rsf_clip <- crop(tmp_rsf,ext)
    tmp_ac <- st_as_sf(x = data.frame(id = 1, 
                                      X = tmp_from[1],
                                      Y = tmp_from[2]),
                       coords = c("X","Y"))
  
    plot(tmp_rsf_clip, axes=FALSE, legend=FALSE, col = viridis(100,option="A", direction=-1, alpha = 0.5))
    plot(st_geometry(st_buffer(tmp_ac,sqrt(5.99)*sigma)), add=TRUE)
    points(coordinates(cost_surface), cex=sqrt(n_by_cell_freq[i,])/2, pch=16, col=adjustcolor(1,0.2))
    plot(tmp_ac,add=TRUE, pch=16, col=2, cex=2)
  }
}


# SCR data ----------------------------------------------------------------



