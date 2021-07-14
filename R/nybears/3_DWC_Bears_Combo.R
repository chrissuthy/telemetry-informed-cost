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
library(colorspace)

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
  diff = 100 * ((r_dens$dwc - r_dens0$dwc)/r_dens0$dwc),
  percent_forest = extract(x = raster::aggregate(forest, 2), y = dwc_bears[,1:2])
)

par(pty="s")
plot(diff ~ cut(percent_forest, breaks=seq(0,1,0.05)), dwc_diff, 
     xlab = "% Forest", ylab = "% Difference (M - iM)", 
     col = alpha(viridis(5)[2], 0.75))


#----Plots----

p2 <- ggplot(r_dens,aes(x=X,y=Y,fill=dwc))+
  geom_raster() +
  labs(title = expression("(b) Density-weighted Connectivity" [M[lcp]]), 
       x = NULL, y = NULL, fill=NULL) +
  theme_bw() +
  scale_fill_distiller(palette = "Purples",direction = 1, limits = c(0,0.19)) +
  theme(legend.position="none", aspect.ratio = 1,
        axis.text = element_text(size=12),
        plot.title = element_text(hjust=0.5, face="bold"))

p3 <- ggplot(r_dens0,aes(x=X,y=Y,fill=dwc))+
  geom_raster() +
  labs(title = expression("(c) Density-weighted Connectivity" [iM[lcp]]), 
       x = NULL, y = NULL, fill=NULL) +
  theme_bw() +
  scale_fill_distiller(palette = "Purples",direction = 1, limits = c(0,0.19)) +
  theme(legend.position="bottom", aspect.ratio = 1,
        axis.text = element_text(size=12),
        plot.title = element_text(hjust=0.5, face="bold"))

p4 <- ggplot(dwc_diff,aes(x=X,y=Y,fill=diff))+
  geom_raster() +
  labs(title = expression("(e) % Difference" ~ (DWC[M[lcp]] - DWC[iM[lcp]])), 
       x = NULL, y = NULL, fill=NULL) +
  theme_bw() +
  # scale_fill_gradient2(high =  "#0080FF", mid = alpha("white",1), 
  #                      low = "#FF0000", midpoint = 0, breaks = c(-0.05, 0, 0.05)) +
  scale_fill_gradient2(high =  "#0080FF", mid = alpha("white",1), 
                       low = "#FF0000", midpoint = 0) +
  #scale_fill_distiller(palette = "RdBu", limits = c(-0.065,0.065)) +
  theme(legend.position="bottom", aspect.ratio = 1,
        axis.text = element_text(size=12),
        plot.title = element_text(hjust=0.5, face="bold"))

p5 <- dwc_diff %>%
  mutate(bin=cut_width(percent_forest, width=0.1, boundary=0)) %>%
  ggplot(., aes(x=percent_forest, y=diff, group = bin)) +
  geom_hline(yintercept = 0, lwd=1) +
  geom_boxplot(color = "black", fill = lighten("black", 0.85), pch = 16,
               outlier.alpha = 0.3) +
  labs(y = expression("% Difference" ~ (DWC[M[lcp]] - DWC[iM[lcp]])),
       x = "% Forest", title = expression("(d) % Difference" ~ (DWC[M[lcp]] - DWC[iM[lcp]]))) +
  theme_bw() +
  theme(legend.position="bottom", aspect.ratio = 1,
        axis.text = element_text(size=12),
        plot.title = element_text(hjust=0.5, face="bold"))

pdwc <- (p2/p3)

pdiff <- p5 / p4


#----Bears----

load("output/nybears/ntel=3_share=F_forest.RData")

for_df <- data.frame(X = coordinates(forest)[,1],
                     Y = coordinates(forest)[,2],
                     forest = values(forest))
traps <- as.data.frame(nybears$traplocs/1000)

tel <- nybears$teldata %>%
  dplyr::select(id=animalid, x1 = X_UTM, y1 = Y_UTM, fix = fixnum) %>%
  mutate(x1 = x1/1000, y1=y1/1000) %>%
  arrange(id, fix)

pbears0 <- ggplot()+
  geom_raster(data=for_df, aes(x=X,y=Y,fill=forest)) +
  scale_fill_gradientn(colors = alpha(c(brewer.pal(5, "Greens")), alpha = 0.75),
                       "% Forest ", limits=c(0,1)) +
  geom_point(data=traps, aes(x=X_UTM, y=Y_UTM), 
             pch=3, color="gray30", stroke=0.7) +
  labs(y = NULL, x = NULL, title = expression("(a) New York black bear data")) +
  geom_path(data=tel, aes(x=x1, y=y1, color=id)) +
  scale_color_manual(values=c("purple2", "blue2", "red2")) +
  theme_bw() + guides(color = F) + 
  theme(legend.position="bottom", aspect.ratio = 1,
        axis.text = element_text(size=12),
        plot.title = element_text(hjust=0.5, face="bold"))

pbears <- (pbears0+theme(legend.position = "none"))/pbears0


#----FINAL PLOT----
    "FINAL PLOT"

#pbears | pdwc | pdiff


#----Final plot 2----

library(gridExtra)

# Convert to grobs
#plots <- list(pbears0, pdwc, pdiff)

# plots <- list(pbears0, p2, p3, p5, p4)
# lst_p <- lapply(plots, ggplotGrob)
# 
# # Plot using gridExtra and grid
# layout_mat <-  matrix(c(1,2,2,3,4,4,6,6,5,5,7,7), byrow = F, ncol = 3, nrow = 4)
# 
# test <- gridExtra::grid.arrange(
#   grobs = list(
#     grid::nullGrob(),
#     lst_p[[1]],  
#     grid::nullGrob(),
#     lst_p[[2]], lst_p[[4]],
#     lst_p[[3]], lst_p[[5]]),
#   layout_matrix = layout_mat)


# devtools::install_github('baptiste/egg')
library(egg)

# test_p1 <- ggarrange(ggplot(), pbears0, ggplot(), ncol=1)
# test_p2 <- ggarrange(p2, p5, p3, p4, nrow = 2, byrow = T)
# 
# grid.arrange(grid::nullGrob(), pbears0, grid::nullGrob(), 
#              test_p2, ncol=2, #widths = c(1,2),
#              layout_matrix = matrix(c(1,2,2,3,4,4,4,4), ncol = 2) )


# grid.arrange(
#   grid::nullGrob(), 
#   set_panel_size(width = unit(1, "in"), height = unit(1, "in"), pbears0), 
#   grid::nullGrob(), 
#   test_p2, 
#   ncol=2, layout_matrix = matrix(c(1,2,2,3,4,4,4,4), ncol = 2) )

ap1 <- set_panel_size(pbears0, width = unit(3, "in"), height = unit(3, "in"))
ap2  <- set_panel_size(p2, width = unit(3, "in"), height = unit(3, "in"))
ap3  <- set_panel_size(p5, width = unit(3, "in"), height = unit(3, "in"))
ap4  <- set_panel_size(p3, width = unit(3, "in"), height = unit(3, "in"))
ap5  <- set_panel_size(p4, width = unit(3, "in"), height = unit(3, "in"))

out <- grid.arrange(
  grid::nullGrob(), 
  ap1,
  grid::nullGrob(), 
  ap2, ap4, ap3, ap5,
  layout_matrix = matrix(c(1,2,2,3,4,4,5,5,6,6,7,7), nrow = 4))

# ggplot2::ggsave(plot = out, filename = "test.pdf", device = "pdf",
#                 path = "/Users/gatesdupont/Desktop", scale = 0.75,
#                 width = 15, height = 11, units = "in")


# ggplot2::ggsave(plot = out, filename = "test.pdf", device = "pdf",
#                 path = "/Users/gatesdupont/Desktop", scale = 1.6,
#                 width = 7.2, height = 5.25, units = "in")

ggplot2::ggsave(plot = out, filename = "Figure3.pdf", device = "pdf", dpi = 600,
                path = "/Users/gatesdupont/Desktop", scale = 1.94,
                width = 6, height = 4.375, units = "in")


