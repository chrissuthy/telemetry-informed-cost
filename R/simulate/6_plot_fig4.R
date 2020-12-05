library(dplyr)
library(raster)
library(ggplot2)
library(tidyr)
library(reshape2)
library(NLMR)
library(sf)
library(gdistance)
library(oSCR)

set.seed(1)

select = dplyr::select
extract = raster::extract
mutate = dplyr::mutate

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

alpha2 <- 2

cost <- exp(alpha2 * landscape0)
tr <- transition(cost,function(x) 1/mean(x),16)

# distance:
ecol_d <- costDistance(tr,
                       fromCoords = as.matrix(ss[,1:2]),
                       toCoords = as.matrix(landscape[,1:2]))

## so finally, lets use a 'sigma scaling factor

sf <- 4 #arbitrary choice
sigma_cost <- sigma*sf

#encounter prob:
ecol_p2 <- exp(-ecol_d^2 / (2*sigma_cost^2))

#is a pixel [y=1]>=0.05
ecol_p2.1 <- ecol_p2 >= 0.05

#a S x L matrix, so the number of TRUE's in a row
# is the number of landscape pixels in that ac's '95% HR'

cost2_ecol_hr2 <- apply(ecol_p2.1,2,sum)

plot(rasterFromXYZ(), 
     col = viridis(100))

r_cost2 <- cbind(landscape[,1:2], cost2_ecol_hr2)
colnames(r_cost2)[3] <- "z"

rm(cost, tr, ecol_d, ecol_p2, ecol_p2.1, cost2_ecol_hr2)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

alpha2 <- 1.4

cost <- exp(alpha2 * landscape0)
tr <- transition(cost,function(x) 1/mean(x),16)

# distance:
ecol_d <- costDistance(tr,
                       fromCoords = as.matrix(ss[,1:2]),
                       toCoords = as.matrix(landscape[,1:2]))

## so finally, lets use a 'sigma scaling factor

sf <- 4 #arbitrary choice
sigma_cost <- sigma*sf

#encounter prob:
ecol_p2 <- exp(-ecol_d^2 / (2*sigma_cost^2))

#is a pixel [y=1]>=0.05
ecol_p2.1 <- ecol_p2 >= 0.05

#a S x L matrix, so the number of TRUE's in a row
# is the number of landscape pixels in that ac's '95% HR'

cost1.4_ecol_hr2 <- apply(ecol_p2.1,2,sum)

plot(rasterFromXYZ(cbind(landscape[,1:2], cost1.4_ecol_hr2)))

r_cost1.4 <- cbind(landscape[,1:2], cost1.4_ecol_hr2)
colnames(r_cost1.4)[3] <- "z"

rm(cost, tr, ecol_d, ecol_p2, ecol_p2.1, cost1.4_ecol_hr2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


r_costDiff <- cbind(landscape[,1:2], r_cost2$z - r_cost1.4$z)
colnames(r_costDiff)[3] <- "z"
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


p1 <- ggplot(data=r_cost2, aes(x=x,y=y,fill=z)) +
  geom_tile() +
  scale_fill_viridis(option = "C", limits = c(0,2000)) +
  theme_void() +
  coord_equal() +
  ggtitle("Telem (truth)") +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust=0.5, face="bold"))

p2 <- ggplot(data=r_cost1.4, aes(x=x,y=y,fill=z)) +
  geom_tile() +
  scale_fill_viridis(option = "C", limits = c(0,2000)) +
  theme_void() +
  coord_equal() +
  ggtitle("No telem (biased)") +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust=0.5, face="bold"))


p3 <- ggplot(data=r_costDiff, aes(x=x,y=y,fill=z)) +
  geom_tile() +
  scale_fill_viridis(option = "D") +
  theme_void() +
  coord_equal() +
  ggtitle("Difference") +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust=0.5, face="bold"))


final <- p1+p2+p3

final

