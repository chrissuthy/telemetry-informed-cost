tel1a <- as.matrix(read.table("output/nybears/ntel=1w1_share=F_MSE_SE_forest.txt")[,-1])
tel1b <- as.matrix(read.table("output/nybears/ntel=1w2_share=F_MSE_SE_forest.txt")[,-1])
tel1c <- as.matrix(read.table("output/nybears/ntel=1w3_share=F_MSE_SE_forest.txt")[,-1])
tel3 <- as.matrix(read.table("output/nybears/ntel=3_share=F_MSE_SE_forest.txt")[,-1])


out.df <- data.frame(
  Collared = rep(c("1","2","3","All"),each=7),
  Parameter = colnames(tel1a),
  Estimate = c(tel1a[1,],tel1b[1,],tel1c[1,],tel3[1,]),
  SE = c(tel1a[2,],tel1b[2,],tel1c[2,],tel3[2,]))
out.df$Upper <- out.df$Estimate + 1.96*out.df$SE
out.df$Lower <- out.df$Estimate - 1.96*out.df$SE

library(ggplot2) 

ggplot(subset(out.df,!(Parameter %in% c("p0","psi"))),
       aes(x=Collared,y=Estimate, color=Collared)) +
  geom_errorbar(aes(ymin=Lower,ymax=Upper,width=0)) +
  geom_point(size=3) +
  facet_wrap(.~Parameter, scales = "free",nrow=2) + theme_bw()
  
