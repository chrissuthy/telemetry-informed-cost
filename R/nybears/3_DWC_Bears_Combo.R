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
library(viridis)
library(RColorBrewer)

source("R/likelihoods/scr_move_cost_predict.R")

load("output/nybears/ntel=3_share=F_forest.RData")
load("output/nybears/nomove_forest.RData")


for_df <- data.frame(X = coordinates(forest)[,1],
                     Y = coordinates(forest)[,2],
                     forest = values(forest))

#----SCR+move----

p <-  mm_forest$estimate

dwc_bears0 <- scr_move_cost_predict(param = p, mod = "gauss", share_sigma = FALSE, 
                                    teldata = bears_teldata, spatdata  = spatdata, 
                                    landscape = forest, scr_ss = ss, K = nybears$K, 
                                    scr_y = nybears$y2d, trap_locs = traps, dist = "lcp", 
                                    popcost=TRUE, popmove=TRUE, fixcost=FALSE, use.sbar=TRUE, 
                                    prj=NULL, predict=TRUE)

r_dens0 <- data.frame(X=dwc_bears0[,1],
                      Y=dwc_bears0[,2],
                      rd = apply(dwc_bears0[,-c(1,2)],1,sum))

cost <- exp(p[1]*forest)
tr1 <- transition(cost,transitionFunction=function(x) 1/mean(x),directions=16)
tr1CorrC <- geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
D <- costDistance(tr1CorrC,coordinates(ss),coordinates(ss)) # Cost distance
P <- plogis(p[5]) * exp(- D^2/(2*exp(p[4])^2))
r_dens0$dwc <- apply(P*r_dens0$rd,2,sum)




#----SCR----

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


dwc_diff <- data.frame(
  X = dwc_bears[,1],
  Y = dwc_bears[,2],
  diff = r_dens$dwc - r_dens0$dwc,
  percent_forest = extract(x = raster::aggregate(forest, 1), y = dwc_bears[,1:2])
)

par(pty="s")
plot(diff~percent_forest, dwc_diff, 
     ylab = "(DWC_scr) – (DWC_scr+move)", xlab = "% Forest",
     pch=20, cex = 0.5, col = alpha(viridis(5)[2], 0.1))


boxplot(r_dens$dwc - r_dens0$dwc, ylab = "(DWC_scr) – (DWC_scr+move)")
abline(h=0)


div.pal <- colorRampPalette(brewer.pal(n = 8, name = "RdBu"))

par(mfrow=c(1,2))

plot(rasterFromXYZ(dwc_diff[,c(1,2,3)]), col = div.pal(1000),
     main = "(DWC_scr) – (DWC_scr+move)",
     zlim=c(-max(abs(range(dwc_diff$diff))), max(abs(range(dwc_diff$diff)))))


plot(forest, "% Forest")


