library(FedData)
library(magrittr)
library(raster)
library(sf)
library(ggplot2)
library(oSCR)
library(dplyr)
library(gdistance)
library(cowplot)
library(patchwork)

source("R/likelihoods/scr_move_cost_predict.R")

load("output/nybears/ntel=3_share=F_forest.RData")
load("output/nybears/nomove_forest.RData")

#--------SCR + Move--------

p <-  mm_forest$estimate

dwc_bears <- scr_move_cost_predict(param = p, mod = "gauss", share_sigma = FALSE, 
                                   teldata = bears_teldata, spatdata  = spatdata, 
                                   landscape = forest, scr_ss = ss, K = nybears$K, 
                                   scr_y = nybears$y2d, trap_locs = traps, dist = "lcp", 
                                   popcost=TRUE, popmove=TRUE, fixcost=FALSE, use.sbar=TRUE, 
                                   prj=NULL, predict=TRUE)

r_dens <- data.frame(X=dwc_bears[,1],
                     Y=dwc_bears[,2],
                     rd = apply(dwc_bears[,-c(1,2)],1,sum))
for_df <- data.frame(X = coordinates(forest)[,1],
                     Y = coordinates(forest)[,2],
                     forest = values(forest))
ggplot(r_dens,aes(x=X,y=Y,fill=rd))+
  geom_raster() +
  theme_bw() +
  scale_fill_distiller(direction = 1)

#DWC

cost <- exp(p[1]*forest)
tr1 <- transition(cost,transitionFunction=function(x) 1/mean(x),directions=16)
tr1CorrC <- geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
D <- costDistance(tr1CorrC,coordinates(ss),coordinates(ss)) # Cost distance
P <- plogis(p[5]) * exp(- D^2/(2*exp(p[4])^2))
r_dens$dwc <- apply(P*r_dens$rd,2,sum)


p0 <- ggplot(for_df,aes(x=X,y=Y,fill=forest))+
  geom_raster() +
  labs(title = "a. Forest", fill=NULL) +
  theme_bw() +
  scale_fill_distiller(palette = "Greens",direction = 1) +
  theme(legend.position="bottom", aspect.ratio = 1,
        axis.text = element_text(size=12),
        plot.title = element_text(hjust=0, face="bold"))

p1 <- ggplot(r_dens,aes(x=X,y=Y,fill=rd))+
  geom_raster() +
  labs(title = "b. Realized density", fill=NULL) +
  theme_bw() +
  scale_fill_distiller(palette = "Blues",direction = 1, limits = c(0,0.14)) +
  theme(legend.position="bottom", aspect.ratio = 1,
        axis.text = element_text(size=12),
        plot.title = element_text(hjust=0, face="bold"))

p2 <- ggplot(r_dens,aes(x=X,y=Y,fill=dwc))+
  geom_raster() +
  labs(title = "c. Density-weighted Connectivity", fill=NULL) +
  theme_bw() +
  scale_fill_distiller(palette = "Reds",direction = 1, limits = c(0,0.19)) +
  theme(legend.position="bottom", aspect.ratio = 1,
        axis.text = element_text(size=12),
        plot.title = element_text(hjust=0, face="bold"))

#plot_grid(p0,p1,p2,nrow = 1)



#--------SCR--------

source("R/likelihoods/scr_cost_predict.R")

p_nomove <- m_forest$estimate

dwc_bears <- scr_cost_predict(
  param = p_nomove, dist=c("euc","circ","lcp")[3], 
  scr_ss = ss, scr_y = nybears$y2d, K = nybears$K, trap_locs = traps, landscape = forest,
  mod=c("exp","gauss")[2], prj = NULL, predict=T)

r_dens <- data.frame(X=dwc_bears[,1],
                     Y=dwc_bears[,2],
                     rd = apply(dwc_bears[,-c(1,2)],1,sum))


#DWC

cost <- exp(p_nomove[1]*forest)
tr1 <- transition(cost,transitionFunction=function(x) 1/mean(x),directions=16)
tr1CorrC <- geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
D <- costDistance(tr1CorrC,coordinates(ss),coordinates(ss)) # Cost distance
P <- plogis(p_nomove[3]) * exp(- D^2/(2*exp(p_nomove[2])^2))
r_dens$dwc <- apply(P*r_dens$rd,2,sum)

p0_bottom <- ggplot(for_df,aes(x=X,y=Y,fill=forest))+
  geom_raster() +
  labs(title = "a. Forest", fill=NULL) +
  theme_bw() +
  scale_fill_distiller(palette = "Greens",direction = 1) +
  theme(legend.position="bottom", aspect.ratio = 1,
        axis.text = element_text(size=12),
        plot.title = element_text(hjust=0, face="bold"))

p3 <- ggplot(r_dens,aes(x=X,y=Y,fill=rd))+
  geom_raster() +
  labs(title = "b. Realized density", fill=NULL) +
  theme_bw() +
  scale_fill_distiller(palette = "Blues",direction = 1, limits = c(0,0.14)) +
  theme(legend.position="bottom", aspect.ratio = 1,
        axis.text = element_text(size=12),
        plot.title = element_text(hjust=0, face="bold"))

p4 <- ggplot(r_dens,aes(x=X,y=Y,fill=dwc))+
  geom_raster() +
  labs(title = "c. Density-weighted Connectivity", fill=NULL) +
  theme_bw() +
  scale_fill_distiller(palette = "Reds",direction = 1, limits = c(0,0.19)) +
  theme(legend.position="bottom", aspect.ratio = 1,
        axis.text = element_text(size=12),
        plot.title = element_text(hjust=0, face="bold"))

#(p0/p0_bottom)
#(p1/p3) + (p2/p4)

(p0+p3+p4)/(p0+p1+p2)
