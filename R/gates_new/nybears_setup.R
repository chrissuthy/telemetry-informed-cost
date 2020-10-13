library(NLMR)
library(dplyr)
library(plotrix)
library(raster)
library(ggplot2)
library(ggforce)
library(patchwork)

#----Simulation settings----

# Movement model
psi <- 0.3
upsilon <- 600
sigma <- 4300

# SCR
N <- 30

# Statespace
ncol <- nrow <- 175
rr <- upsilon
autocorr <- 7

# Derived 
hr95 <- sqrt(5.99) * sigma
hr95_lim <- 1.5*hr95


#----Landscape----

landscape0 <- nlm_gaussianfield(
  ncol = ncol, nrow = nrow, 
  resolution = rr, autocorr_range = autocorr)
landscape <- landscape0^2 # I STILL DON'T LIKE THIS! RESTRICTS LOW COST
landscape_df <- as.data.frame(landscape, xy=T)

ss <- coordinates(landscape0)
nrow(ss)/6400


#----Activity centers----

ss_avail <- ss %>%
  as_tibble() %>%
  filter(x > (min(x)+hr95_lim) & x < (max(x)-hr95_lim)) %>%
  filter(y > (min(y)+hr95_lim) & y < (max(y)-hr95_lim)) %>%
  as.data.frame()
  
acs <- ss_avail[sample(1:nrow(ss_avail), size = N),]
acs_df <- as.data.frame(acs)


#----Plotting----

p1 <- ggplot() +
  geom_tile(data = landscape_df, aes(x=x, y=y, fill=layer)) +
  scale_fill_viridis_c("Conductance") +
  geom_circle(data = acs_df, aes(x0=x, y0=y, r = hr95), 
              lwd = 0.2, color = alpha("white", 0.65)) +
  coord_equal() +
  theme_minimal()

text1 <- paste0(
  "N = ", N, "\n",
  "Psi = ", psi, "\n",
  "Upsilon = ", upsilon, "\n",
  "Sigma = ", sigma, "\n\n",
  "pixels = ", nrow(ss), "\n",
  "resolution =", "upsilon"
)

p2 <- ggplot() +
  annotate("text", x = 4, y = 25, size=8, label = text1) + 
  theme_void()

p1+p2
  



