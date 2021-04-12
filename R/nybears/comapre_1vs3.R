library(RColorBrewer)
library(ggplot2) 

tel1a <- as.matrix(read.table("output/nybears/ntel=1w1_share=F_MSE_SE_forest.txt")[,-1])
tel1b <- as.matrix(read.table("output/nybears/ntel=1w2_share=F_MSE_SE_forest.txt")[,-1])
tel1c <- as.matrix(read.table("output/nybears/ntel=1w3_share=F_MSE_SE_forest.txt")[,-1])
tel3 <- as.matrix(read.table("output/nybears/ntel=3_share=F_MSE_SE_forest.txt")[,-1])


params <- c("alpha[1]", "sigma[step]", "psi", "sigma[det]", "p[0]", "log(lambda)", "sigma[home]")

out.df <- data.frame(
  Collared = rep(c("A","B","C","All"),each=7),
  #Parameter = colnames(tel1a),
  Parameter = params,
  Estimate = c(tel1a[1,],tel1b[1,],tel1c[1,],tel3[1,]),
  SE = c(tel1a[2,],tel1b[2,],tel1c[2,],tel3[2,]))
out.df$Collared <- factor(out.df$Collared, levels = c("A","B","C","All"))
out.df$Upper <- out.df$Estimate + 1.96*out.df$SE
out.df$Lower <- out.df$Estimate - 1.96*out.df$SE

p <- ggplot(subset(out.df,!(Parameter %in% c("p[0]","psi"))),
       aes(x=Collared,y=Estimate, color=Collared, shape = Collared)) +
  geom_errorbar(aes(ymin=Lower,ymax=Upper,width=0)) +
  geom_point(size=3) +
  scale_color_manual(values = brewer.pal(4, "Set1"), guide=F) +
  scale_shape_manual(values = c(16,16,16,1), guide = FALSE) +
  labs(x = NULL) +
  facet_wrap(.~Parameter, scales = "free",nrow=2, labeller = label_parsed) + 
  theme_minimal() +
  theme(aspect.ratio = 1,
        text = element_text(size = 14),
        strip.text.y = element_text(angle = 0, hjust = 0),
        #axis.text.x = element_blank(),
        panel.border = element_rect(fill = NA, color = "gray70", size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank());p

ggsave(filename = "AppendixS2FigS2.pdf", plot = p, device="pdf",
       dpi = 600, scale = 0.8, height = 5.25, width = 8.75,
       path = "/Users/gatesdupont/Desktop")
  
