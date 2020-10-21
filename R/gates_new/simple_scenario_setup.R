library(dplyr)
library(raster)
library(ggplot2)
library(tidyr)
library(reshape2)
library(NLMR)
library(sf)
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
upsilon <- 0.5 #0.6
sigma <- 2  #4.3

# SCR
N <- 50

# Statespace
ncol <- nrow <- 249 #175
rr <- upsilon/2
autocorr <- 6

# Derived 
hr95 <- sqrt(5.99) * sigma
hr95_lim <- (3*sigma) #+ (3*upsilon)

landscape0 <- nlm_gaussianfield(
  ncol = ncol, nrow = nrow, 
  resolution = rr, autocorr_range = autocorr)
landscape_r <- landscape0^2 # I STILL DON'T LIKE THIS! RESTRICTS LOW COST
landscape <- as.data.frame(landscape_r, xy=T)

ss <- landscape %>%
  filter(x >= (min(x)+hr95_lim) & x < (max(x)-hr95_lim)) %>%
  filter(y >= (min(y)+hr95_lim) & y < (max(y)-hr95_lim))
nrow(ss)


#----Activity centers----

acs <- ss[sample(1:nrow(ss), size = N, replace = T),c("x","y")]
rownames(acs) <- NULL
acs_df <- as.data.frame(acs)


#----Traps----

# Trap array (space)
trap_array <- ss %>%
  filter(x >= (min(x)+3*sigma) & x < (max(x)-3*sigma)) %>%
  filter(y >= (min(y)+3*sigma) & y < (max(y)-3*sigma))

# Make traps
trap_array_sf <- st_as_sf(trap_array, coords = c('y', 'x'))
traps <- st_make_grid(trap_array_sf, n=c(10,10)) %>% # alternatively could use cellsize
  as_Spatial() %>%
  coordinates() %>%
  as.data.frame() %>%
  select(x = V2, y = V1)
rownames(traps) <- NULL
table(diff(traps$y))


#----Plot----
ggplot() +
  geom_tile(data = landscape, aes(x=x, y=y, fill=layer)) +
  scale_fill_viridis_c("Landscape") +
  geom_rect(aes(xmin=min(ss$x), xmax=max(ss$x), ymin=min(ss$y), ymax=max(ss$y)),
            fill=alpha("white",0.25)) +
  geom_circle(data = acs_df, aes(x0=x, y0=y, r = hr95),
              lwd = 0.2, color = alpha("white", 0.65)) +
  geom_point(data = traps, aes(x=x, y=y), fill = "black", color = "black", pch = 21) +
  coord_equal() +
  theme_minimal()
