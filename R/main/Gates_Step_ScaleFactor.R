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
upsilon0 <- 0.5
upsilon <- upsilon0
#upsilon <- 0.5 #bears: 0.6
sigma0 <- 2
sigma <- sigma0
#sigma <- 2  #bears: 4.3

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



#----Calculate pixels in every possibly step kernel for 10 individuals----
ids <- sample(1:nrow(ss), size = 10)

id_euc_hr_ALL <- list()
id_ecol_hr_ALL <- list()
id_ecol2_hr_ALL <- list()

for(i in 1:length(ids)){
  
  #----All steps in home range----
  
  # For a single home range / activity center
  id_px <- ids[i]
  # plot(rasterFromXYZ(cbind(landscape[,1:2], euc_p[id_px,])))
  
  sigma <- sigma0
  upsilon <- upsilon0
  
  # Select all home range pixels
  id_hr <- cbind(landscape[,1:2], hr = euc_p[id_px,]) %>%
    mutate(pixel = 1:n()) %>%
    filter(hr == TRUE) %>%
    select(x,y, pixel)
  # points(id_hr[,1:2], col = "red")
  
  # Calculate step length from first home range pixel
  id_euc_d_step <- e2dist(id_hr[,1:2],landscape[,1:2])
  # plot(rasterFromXYZ(cbind(landscape[,1:2], id_euc_d_step[1,])))
  # points(id_hr %>% filter(pixel == id_hr[1,3]) %>% select(x,y))
  
  # Calculate step length from activity center
  id_euc_d_ac <- e2dist(ss[id_px,1:2], landscape[,1:2]) %x% rep(1, nrow(id_euc_d_step))
  # plot(rasterFromXYZ(cbind(landscape[,1:2], id_euc_d_ac[1,])))
  
  # Calculate probability of stepping in each landscape pixel
  id_euc_p <- exp((-id_euc_d_step*id_euc_d_step/(2*upsilon*upsilon))  - (id_euc_d_ac*id_euc_d_ac/(2*sigma*sigma)) )
  cbind(landscape[,1:2], probs = id_euc_p[200,]) %>%
    as.data.frame() %>%
    ggplot(data=.) +
    geom_tile(aes(x=x, y=y, fill = probs)) +
    #geom_point(data = id_hr[200,1:2], aes(x=x, y=y)) +
    #scale_fill_viridis_c() +
    coord_equal() +
    theme_minimal()
  
  # Convert probaiblities to T/F
  id_euc_p <- id_euc_p >= 0.05
  # plot(rasterFromXYZ(cbind(landscape[,1:2], id_euc_p[1,])))
  id_euc_hr <- apply(id_euc_p,1,sum)
  id_euc_hr
  id_euc_hr_ALL[[i]] <- id_euc_hr
  
  # boxplot(id_euc_hr)
  
  
  #----All steps in a home range with ecol dist----
  
  sf <- 4
  sigma <- sigma0*sf
  upsilon <- upsilon0
  
  # For a single home range / activity center
  id_px
  # plot(rasterFromXYZ(cbind(landscape[,1:2], euc_p[id_px,])))
  
  # Get home range with eco dist
  alpha2 <- 2
  cost <- exp(alpha2 * landscape0)
  tr <- transition(cost,function(x) 1/mean(x),16)
  ecol_d <- costDistance(tr,
                         fromCoords = as.matrix(ss[id_px,1:2]),
                         toCoords = as.matrix(landscape[,1:2]))
  ecol_p <- exp(-ecol_d^2 / (2*sigma^2))
  ecol_p <- ecol_p >= 0.05
  ecol_hr <- apply(ecol_p,1,sum)
  
  # Select all home range pixels
  id_hr <- cbind(landscape[,1:2], hr = ecol_p[1,]) %>%
    mutate(pixel = 1:n()) %>%
    filter(hr == TRUE) %>%
    select(x,y, pixel)
  #plot(rasterFromXYZ(cbind(landscape[,1:2], ecol_p[1,])))
  #points(id_hr[,1:2], col = "red")
  
  
  # Calculate step length from first home range pixel
  id_ecol_d_step <- costDistance(tr,
                                 fromCoords = as.matrix(id_hr[,1:2]),
                                 toCoords = as.matrix(landscape[,1:2]))
  #plot(rasterFromXYZ(cbind(landscape[,1:2], id_ecol_d_step[1,])))
  # points(id_hr %>% filter(pixel == id_hr[1,3]) %>% select(x,y))
  
  # Calculate step length from activity center
  id_ecol_d_ac <- costDistance(tr,
                               fromCoords = as.matrix(ss[id_px,1:2]),
                               toCoords = as.matrix(landscape[,1:2])) %x% rep(1, nrow(id_ecol_d_step))
  #plot(rasterFromXYZ(cbind(landscape[,1:2], id_ecol_d_ac[1,])))
  
  # Calculate probability of stepping in each landscape pixel
  id_ecol_p <- exp((-id_ecol_d_step*id_ecol_d_step/(2*upsilon*upsilon))  - (id_ecol_d_ac*id_ecol_d_ac/(2*sigma*sigma)) )
  # cbind(landscape[,1:2], probs = id_ecol_p[50,]) %>%
  #   as.data.frame() %>%
  #   ggplot(data=.) +
  #   geom_tile(aes(x=x, y=y, fill = probs)) +
  #   scale_fill_viridis_c() +
  #   coord_equal() +
  #   theme_minimal()
  
  
  # Convert probaiblities to T/F
  id_ecol_p <- id_ecol_p >= 0.05
  # plot(rasterFromXYZ(cbind(landscape[,1:2], id_ecol_p[1,])))
  id_ecol_hr <- apply(id_ecol_p,1,sum)
  id_ecol_hr
  
  id_ecol_hr_ALL[[i]] <- id_ecol_hr
  
  # boxplot(id_euc_hr, id_ecol_hr)
  
  
  
  
  #----All steps in a home range with ecol dist + scale factor----
  
  sf <- 4
  sigma <- sigma0*sf
  
  
  sf2 <- 4 #arbitrary choice
  upsilon <- upsilon0*sf2
  
  
  # For a single home range / activity center
  id_px
  # plot(rasterFromXYZ(cbind(landscape[,1:2], euc_p[id_px,])))
  
  # Get home range with eco dist
  alpha2 <- 2
  cost <- exp(alpha2 * landscape0)
  tr <- transition(cost,function(x) 1/mean(x),16)
  ecol2_d <- costDistance(tr,
                          fromCoords = as.matrix(ss[id_px,1:2]),
                          toCoords = as.matrix(landscape[,1:2]))
  ecol2_p <- exp(-ecol2_d^2 / (2*sigma^2))
  ecol2_p <- ecol2_p >= 0.05
  ecol2_hr <- apply(ecol2_p,1,sum)
  
  # Select all home range pixels
  id_hr <- cbind(landscape[,1:2], hr = ecol2_p[1,]) %>%
    mutate(pixel = 1:n()) %>%
    filter(hr == TRUE) %>%
    select(x,y, pixel)
  # plot(rasterFromXYZ(cbind(landscape[,1:2], ecol2_p[1,])))
  # points(id_hr[,1:2], col = "red")
  
  
  # Calculate step length from first home range pixel
  id_ecol2_d_step <- costDistance(tr,
                                  fromCoords = as.matrix(id_hr[,1:2]),
                                  toCoords = as.matrix(landscape[,1:2]))
  # plot(rasterFromXYZ(cbind(landscape[,1:2], id_ecol2_d_step[1,])))
  # points(id_hr %>% filter(pixel == id_hr[1,3]) %>% select(x,y))
  
  # Calculate step length from activity center
  id_ecol2_d_ac <- costDistance(tr,
                                fromCoords = as.matrix(ss[id_px,1:2]),
                                toCoords = as.matrix(landscape[,1:2])) %x% rep(1, nrow(id_ecol2_d_step))
  # plot(rasterFromXYZ(cbind(landscape[,1:2], id_ecol2_d_ac[1,])))
  
  # Calculate probability of stepping in each landscape pixel
  id_ecol2_p <- exp((-id_ecol2_d_step*id_ecol2_d_step/(2*upsilon*upsilon))  - (id_ecol2_d_ac*id_ecol2_d_ac/(2*sigma*sigma)) )
  # cbind(landscape[,1:2], probs = id_ecol2_p[50,]) %>%
  #   as.data.frame() %>%
  #   ggplot(data=.) +
  #   geom_tile(aes(x=x, y=y, fill = probs)) +
  #   scale_fill_viridis_c() +
  #   coord_equal() +
  #   theme_minimal()
  
  
  # Convert probaiblities to T/F
  id_ecol2_p <- id_ecol2_p >= 0.05
  # plot(rasterFromXYZ(cbind(landscape[,1:2], id_ecol2_p[1,])))
  id_ecol2_hr <- apply(id_ecol2_p,1,sum)
  id_ecol2_hr
  
  id_ecol2_hr_ALL[[i]] <- id_ecol2_hr
  
}


#----Bring it all together----

euc_steps <- unlist(id_euc_hr_ALL)
eco_steps <- unlist(id_ecol_hr_ALL)
ecosf_steps <- unlist(id_ecol2_hr_ALL)

df0 <- data.frame(
  Scenario = c(rep("Euc", length(euc_steps)), 
               rep("Eco", length(eco_steps)),
               rep("Eco+sf", length(ecosf_steps))),
  Pixels = c(euc_steps, eco_steps, ecosf_steps))

df <- df0 %>%
  mutate(Scenario = factor(Scenario, levels = c("Euc", "Eco", "Eco+sf")))

ggplot(data=df, aes(y = Pixels, x = Scenario, fill = Scenario)) +
  geom_boxplot()+
  labs(y = "Pixels used", title = "Upsilon scenarios") +
  theme_minimal() +
  theme(aspect.ratio = 1, legend.position = "none")






