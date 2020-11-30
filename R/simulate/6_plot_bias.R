library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(patchwork)

#----True Values----

N        <- 100
sspixels <- 625
t_alpha2 <- 1
#t_ups    <- 0.25*2.5
t_psi    <- 0.9
t_sig    <- 2.5
t_d0     <- N/sspixels
t_sig_mm <- 2.5


#----Custom functions----

# 50% quantile
ggquant <- function(y, int){
  Lint <- 0.5-(int/2)
  Uint <- 0.5+(int/2)
  return(
    data.frame(
      ymin = quantile(y, Lint), 
      ymax = quantile(y, Uint)))
}

meantrim <- function(y){
  result <- mean(y, trim = 0.1)
  return(result)
}

# Analyze results files
analyze <- function(df, sigma, scenario, scenario_ups, t_ups){
  if(sigma == "unshared"){
    result <- df %>% 
      as.data.frame() %>%
      select(-p0) %>%
      gather() %>%
      mutate(true = rep(c(t_alpha2, t_ups, t_psi, t_sig, t_d0, t_sig_mm), each = nrow(.)/6)) %>%
      mutate(prbias = 100*((value-true)/true)) %>%
      mutate(key = recode(key, alpha2 = "cost", d0 = "density", sig = "sigma", 
                          upsilon = "upsilon", psi = "pr(moved)", sig_mm = "sigma_move")) %>%
      mutate(key = factor(key, levels = c("cost", "density", "sigma", 
                                          "sigma_move", "upsilon", "pr(moved)"))) %>%
      mutate(Scenario = scenario) %>%
      mutate(Upsilon = scenario_ups)
  } else if(sigma == "shared") {
    result <- df %>%
      as.data.frame() %>%
      select(-p0) %>%
      gather() %>%
      mutate(true = rep(c(t_alpha2, t_ups, t_psi, t_sig, t_d0), each = nrow(.)/5)) %>%
      mutate(prbias = 100*((value-true)/true)) %>%
      mutate(key = recode(key, alpha2 = "cost", d0 = "density", sig = "sigma", 
                          upsilon = "upsilon", psi = "pr(moved)")) %>%
      mutate(key = factor(key, levels = c("cost", "density", "sigma", 
                                          "upsilon", "pr(moved)"))) %>%
      mutate(Scenario = scenario) %>%
      mutate(Upsilon = scenario_ups)
  }
  return(result)
}


#----Analyze results----

# No movement
noMove <- read.table("./output/small ups/est_noMove/screco_results.txt") %>%
  as.data.frame() %>%
  select(-p0) %>%
  gather() %>%
  mutate(true = rep(c(t_alpha2, t_sig, t_d0), each = nrow(.)/3)) %>%
  mutate(prbias = 100*((value-true)/true)) %>%
  mutate(key = recode(key, alpha2 = "cost", d0 = "density", sig = "sigma")) %>%
  mutate(key = factor(key, levels = c("cost", "density", "sigma"))) %>%
  mutate(Scenario = "No movement") %>%
  mutate(Upsilon = "small")
  
# NTEL = 1
ntel1_unshared <- read.table("./output/small ups/est_ntel=1_share=F/results.txt") %>% 
  analyze(., sigma = "unshared", scenario = "ntel=1, unshared", scenario_ups = "small", t_ups = 0.25*2.5)
ntel1_shared <- read.table("./output/small ups/est_ntel=1_share=T/results.txt") %>% 
  analyze(., sigma = "shared", scenario = "ntel=1, shared", scenario_ups = "small", t_ups = 0.25*2.5)

# NTEL = 3
ntel3_unshared <- read.table("./output/small ups/est_ntel=3_share=F/results.txt") %>%
  analyze(., sigma = "unshared", scenario = "ntel=3, unshared", scenario_ups = "small", t_ups = 0.25*2.5)
ntel3_shared <- read.table("./output/small ups/est_ntel=3_share=T/results.txt") %>% 
  analyze(., sigma = "shared", scenario = "ntel=3, shared", scenario_ups = "small", t_ups = 0.25*2.5)

# NTEL = 5
ntel5_unshared <- read.table("./output/small ups/est_ntel=5_share=F/results.txt") %>% 
  analyze(., sigma = "unshared", scenario = "ntel=5, unshared", scenario_ups = "small", t_ups = 0.25*2.5)
ntel5_shared <- read.table("./output/small ups/est_ntel=5_share=T/results.txt") %>% 
  analyze(., sigma = "shared", scenario = "ntel=5, shared", scenario_ups = "small", t_ups = 0.25*2.5)


#----Combine data----

# Combine
df <- rbind(
  noMove,
  ntel1_unshared, ntel1_shared,
  ntel3_unshared, ntel3_shared,
  ntel5_unshared, ntel5_shared) %>%
  mutate(Scenario = factor(Scenario, 
                           levels = c("No movement", 
                           "ntel=1, shared", "ntel=1, unshared", 
                           "ntel=3, shared", "ntel=3, unshared",  
                           "ntel=5, shared", "ntel=5, unshared")))


#----Plot every scenario----

# Color palettes
orange  <- brewer.pal(n = 9, name = "Oranges")
green   <- brewer.pal(n = 9, name = "Greens")
blue    <- brewer.pal(n = 9, name = "Blues")
purple <- brewer.pal(n = 9, name = "Purples")

# Plots
ggplot(data=df,  aes(x = key, y = prbias, color = Scenario)) +
  geom_rect(ymin = -5, ymax = 5, xmin = 0, xmax = 10, 
            fill = "gray91", color = "gray91") +
  geom_hline(yintercept = 0, color = "gray40", size = 0.7) +
  stat_summary(fun.data=ggquant, fun.args = list(int=0.5), lwd=1, 
               geom="errorbar", position = position_dodge(0.85), width = 0) +
  stat_summary(fun = meantrim, geom="point", position = position_dodge(0.85), size=2.0, aes(shape=Scenario)) +
  scale_color_manual(values = c(orange[5], green[6], green[8], blue[6], blue[8], purple[6], purple[8])) +
  scale_shape_manual(values = c(15, 16, 17, 16, 17, 16, 17)) +
  facet_wrap(~key, scales = "free_x", nrow=1) +
  theme_minimal() + labs(y = "% Relaive bias", x = NULL) +
  coord_cartesian(ylim = c(-60, 60)) +
  theme(#aspect.ratio = 1,
    strip.text.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    plot.title = element_blank(),
    text = element_text(size = 14))


#----Plot cost and density comparisons----

# COST
p1 <- df %>%
  filter(key == "cost") %>%
  ggplot(data=.,  aes(x = key, y = prbias, color = Scenario)) +
  geom_rect(ymin = -5, ymax = 5, xmin = 0, xmax = 10, 
            fill = "gray91", color = "gray91") +
  geom_hline(yintercept = 0, color = "gray40", size = 0.7) +
  stat_summary(fun.data=ggquant, fun.args = list(int=0.5), lwd=1,
               geom="errorbar", position = position_dodge(0.85), width = 0) +
  stat_summary(fun = meantrim, geom="point", position = position_dodge(0.85), size=2.0, aes(shape=Scenario)) +
  scale_color_manual(values = c(orange[5], green[6], green[8], blue[6], blue[8], purple[6], purple[8])) +
  scale_shape_manual(values = c(15, 16, 17, 16, 17, 16, 17)) +
  theme_minimal() + labs(y = "% Relaive bias", x = NULL) +
  coord_cartesian(ylim=c(-150, 150)) +
  theme(aspect.ratio = 1, legend.position = "none",
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(color = "gray80"),
        strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        plot.title = element_blank(),
        text = element_text(size = 14))

# DENSITY
p2 <- df %>%
  filter(key == "density") %>%
  ggplot(data=.,  aes(x = key, y = prbias, color = Scenario)) +
  geom_rect(ymin = -5, ymax = 5, xmin = 0, xmax = 10, 
            fill = "gray91", color = "gray91") +
  geom_hline(yintercept = 0, color = "gray40", size = 0.7) +
  stat_summary(fun.data=ggquant, fun.args = list(int=0.5), lwd=1,
               geom="errorbar", position = position_dodge(0.85), width = 0) +
  stat_summary(fun = meantrim, geom="point", position = position_dodge(0.85), size=2.0, aes(shape=Scenario)) +
  scale_color_manual(values = c(orange[5], green[6], green[8], blue[6], blue[8], purple[6], purple[8])) +
  scale_shape_manual(values = c(15, 16, 17, 16, 17, 16, 17)) +
  theme_minimal() + labs(y = "% Relaive bias", x = NULL) +
  coord_cartesian(ylim=c(-25, 25)) +
  theme(aspect.ratio = 1,
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(color = "gray80"),
        strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        plot.title = element_blank(),
        text = element_text(size = 14))

# Plot both together
p1+p2

