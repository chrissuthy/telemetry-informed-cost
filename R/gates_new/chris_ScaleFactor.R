library(dplyr)
library(raster)
library(ggplot2)
library(tidyr)
library(reshape2)
library(NLMR)
library(sf)
library(gdistance)
library(oSCR)

select = dplyr::select
extract = raster::extract
mutate = dplyr::mutate

# Manual settings
nfix = 90*24
sims = 10

# Cost
alpha2 <- 2

# Movement model
psi <- 0.9
upsilon <- 0.5 #bears: 0.6
sigma <- 2  #bears: 4.3

# SCR
N <- 50

# Statespace
ncol <- nrow <- 149
rr <- upsilon/2
autocorr <- 6

# Derived 
hr95 <- sqrt(5.99) * sigma
hr95_lim <- (3*sigma)

#----Landscape----

# Make the landscape
landscape0 <- nlm_gaussianfield(
  ncol = ncol, nrow = nrow, 
  resolution = rr, autocorr_range = autocorr)
landscape_r <- landscape0^2 # I STILL DON'T LIKE THIS! RESTRICTS LOW COST
landscape <- as.data.frame(landscape_r, xy=T)

# Subset the state-space
ss <- landscape %>%
  filter(x >= (min(x)+hr95_lim) & x < (max(x)-hr95_lim)) %>%
  filter(y >= (min(y)+hr95_lim) & y < (max(y)-hr95_lim))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## lets explore this sigma nonsense

# first, forget sigma altogether, and lets think of 
# 95% home range sze, and lets use number of pixels 
# as the area metric:

# distance:
euc_d <- e2dist(ss[,1:2],landscape[,1:2])

#encounter prob:
euc_p <- exp(-euc_d^2 / (2*sigma^2))

#is a pixel [y=1]>=0.05
euc_p <- euc_p >= 0.05

#a S x L matrix, so the number of TRUE's in a row
# is the number of landscape pixels in that ac's '95% HR'

euc_hr <- apply(euc_p,1,sum)

boxplot(euc_hr)

## not what happens if we do the same, but
## when there is a resistance surface (a2=2)
## the below is the same as above

alpha2 <- 2

cost <- exp(alpha2 * landscape0)
tr <- transition(cost,function(x) 1/mean(x),16)

# distance:
ecol_d <- costDistance(tr,
                       fromCoords = as.matrix(ss[,1:2]),
                       toCoords = as.matrix(landscape[,1:2]))

#encounter prob:
ecol_p <- exp(-ecol_d^2 / (2*sigma^2))

#is a pixel [y=1]>=0.05
ecol_p <- ecol_p >= 0.05

#a S x L matrix, so the number of TRUE's in a row
# is the number of landscape pixels in that ac's '95% HR'

ecol_hr <- apply(ecol_p,1,sum)

boxplot(ecol_hr)

boxplot(euc_hr,ecol_hr) #massive reduction in sace use!

## so finally, lets use a 'sigma scaling factor

sf <- 4 #arbitrary choice
sigma_cost <- sigma*sf

#encounter prob:
ecol_p2 <- exp(-ecol_d^2 / (2*sigma_cost^2))

#is a pixel [y=1]>=0.05
ecol_p2 <- ecol_p2 >= 0.05

#a S x L matrix, so the number of TRUE's in a row
# is the number of landscape pixels in that ac's '95% HR'

ecol_hr2 <- apply(ecol_p2,1,sum)

boxplot(data.frame(Euc = euc_hr,
                   Ecol_1 = ecol_hr,
                   Ecol_2 = ecol_hr2),
        ylab="Pixels used",
        xlab="Scenario")



par(mfrow=c(1,3))
plot(rasterFromXYZ(cbind(landscape[,1:2], euc_p[4000,])))
plot(rasterFromXYZ(cbind(landscape[,1:2], ecol_p[4000,])))
plot(rasterFromXYZ(cbind(landscape[,1:2], ecol_p2[4000,])))
par(mfrow=c(1,1))

