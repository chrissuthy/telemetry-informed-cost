library(dplyr)
library(raster)
library(ggplot2)
library(tidyr)
library(reshape2)
library(NLMR)
library(sf)
library(gdistance)
library(oSCR)
library(viridis)
library(patchwork)
#devtools::install_github("zeehio/facetscales")
library(facetscales)

select = dplyr::select
extract = raster::extract
mutate = dplyr::mutate

set.seed(1)


#----Starting parameters----

# Get surface and ss
load("~/R/telemetry-informed-cost/output/nybears/nomove_forest.RData")
forest
ss
forest_df <- as.data.frame(forest, xy=T)


#----Landscape----

# Make the landscape
# landscape0 <- nlm_gaussianfield(
#   ncol = ncol, nrow = nrow, 
#   resolution = rr, autocorr_range = autocorr)
# landscape_r <- landscape0^2 # I STILL DON'T LIKE THIS! RESTRICTS LOW COST
# landscape <- as.data.frame(landscape_r, xy=T)


#----Unbiased cost scenario----

alpha2 <- -1.172

cost2 <- exp(alpha2 * forest)
tr <- transition(cost2,function(x) 1/mean(x),16)

# distance:
ecol_d <- costDistance(tr,
                       fromCoords = as.matrix(forest_df[,1:2]),
                       toCoords = as.matrix(forest_df[,1:2]))

## so finally, lets use a 'sigma scaling factor
sigma_cost <- exp(0.8116)

#encounter prob:
ecol_p2 <- exp(-ecol_d^2 / (2*sigma_cost^2))

#is a pixel [y=1]>=0.05
ecol_p2.1 <- ecol_p2 >= 0.05

#a S x L matrix, so the number of TRUE's in a row
# is the number of landscape pixels in that ac's '95% HR'

cost2_ecol_hr2 <- apply(ecol_p2.1,2,sum)

plot(rasterFromXYZ(cbind(forest_df[,1:2], cost2_ecol_hr2)))

r_cost2 <- cbind(forest_df[,1:2], cost2_ecol_hr2)
colnames(r_cost2)[3] <- "z"

rm(tr, ecol_d, ecol_p2, ecol_p2.1, cost2_ecol_hr2)


#----Biased cost scenario----

alpha2 <- 0.2261

cost1.4 <- exp(alpha2 * forest)
tr <- transition(cost1.4,function(x) 1/mean(x),16)

# distance:
ecol_d <- costDistance(tr,
                       fromCoords = as.matrix(forest_df[,1:2]),
                       toCoords = as.matrix(forest_df[,1:2]))

## so finally, lets use a 'sigma scaling factor

sigma_cost <- exp(1.2845)

#encounter prob:
ecol_p2 <- exp(-ecol_d^2 / (2*sigma_cost^2))

#is a pixel [y=1]>=0.05
ecol_p2.1 <- ecol_p2 >= 0.05

#a S x L matrix, so the number of TRUE's in a row
# is the number of landscape pixels in that ac's '95% HR'

cost1.4_ecol_hr2 <- apply(ecol_p2.1,2,sum)

plot(rasterFromXYZ(cbind(forest_df[,1:2], cost1.4_ecol_hr2)))

r_cost1.4 <- cbind(forest_df[,1:2], cost1.4_ecol_hr2)
colnames(r_cost1.4)[3] <- "z"

rm(tr, ecol_d, ecol_p2, ecol_p2.1, cost1.4_ecol_hr2)


#----Compile the data----

# Cost surfaces
cost1.4_df <- as.data.frame(cost1.4, xy=T)
colnames(cost1.4_df)[3] <- "z"

cost2_df <- as.data.frame(cost2, xy=T)
colnames(cost2_df)[3] <- "z"

# Full data frame
df <- rbind(
  r_cost2, r_cost1.4,
  cost2_df, cost1.4_df) %>%
  mutate(surface = rep(c("Potential\nconnectivity", "Cost"), 
                       each = nrow(cost2_df)*2)) %>%
  mutate(result = rep(c("SCR[eco]+MM[eco]", "SCR[eco]",
                        "SCR[eco]+MM[eco]", "SCR[eco]"),
                      each = nrow(cost2_df))) %>%
  mutate(result = factor(result, levels=c("SCR[eco]+MM[eco]","SCR[eco]")))


#----Final plot----

p1 <- df %>%
  filter(surface == "Cost") %>%
  ggplot(data=., aes(x=x,y=y,fill=z)) +
  geom_tile() +
  facet_grid(result~surface) +
  scale_fill_viridis(option = "D") + #, limits = c(0, max(cost2_df$z))) +
  theme_void() + coord_equal() +
  labs(title = "Cost") +
  theme(legend.position = "none",
        #plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        #plot.background = element_rect(fill = "yellow"),
        strip.text.y = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(hjust=0.5, face="bold", size=15))
        #strip.text.x = element_text(face="bold", size=15))

p2 <- df %>%
  filter(surface == "Potential\nconnectivity") %>%
  ggplot(data=., aes(x=x,y=y,fill=z)) +
  geom_tile() +
  facet_grid(result~surface, labeller = label_parsed) +
  scale_fill_viridis(option = "C") + #, limits = c(0, max(r_cost1.4$z))) +
  theme_void() + coord_equal() +
  labs(title = "Potential\nconnectivity") +
  theme(legend.position = "none",
        #plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        #plot.background = element_rect(fill = "yellow"),
        strip.text.y = element_text(face="bold", size=15),
        strip.text.x = element_blank(),
        plot.title = element_text(hjust=0.5, face="bold", size=15),
        # strip.text.x = element_text(hjust=0.5, face="bold", size=15),
        text = element_text(face="bold"))



final <- p1+p2
final




