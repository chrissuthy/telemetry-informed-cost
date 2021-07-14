library(FedData)
library(magrittr)
library(raster)
library(sf)
library(ggplot2)
library(oSCR)
library(dplyr)
library(gdistance)
library(cowplot)
library(RColorBrewer)

source("R/likelihoods/scr_move_cost_predict.R")

load("output/nybears/ntel=3_share=F_forest.RData")

for_df <- data.frame(X = coordinates(forest)[,1],
                     Y = coordinates(forest)[,2],
                     forest = values(forest))
traps <- as.data.frame(nybears$traplocs/1000)

tel <- nybears$teldata %>%
  select(id=animalid, x1 = X_UTM, y1 = Y_UTM, fix = fixnum) %>%
  mutate(x1 = x1/1000, y1=y1/1000) %>%
  arrange(id, fix)

pbears <- ggplot()+
  geom_raster(data=for_df, aes(x=X,y=Y,fill=forest)) +
  scale_fill_gradientn(colors = alpha(c(brewer.pal(5, "Greens")), alpha = 0.75),
                       "% Forest ", limits=c(0,1)) +
  geom_point(data=traps, aes(x=X_UTM, y=Y_UTM), 
             pch=3, color="gray30", stroke=0.7) +
  labs(y = NULL, x = NULL, title = expression("New York black bear data")) +
  geom_path(data=tel, aes(x=x1, y=y1, color=id)) +
  scale_color_manual(values=c("purple2", "blue2", "red2")) +
  theme_bw() + guides(color = F) + 
  theme(legend.position="bottom", 
        plot.title = element_text(hjust=0.5),
        aspect.ratio = 1)
  