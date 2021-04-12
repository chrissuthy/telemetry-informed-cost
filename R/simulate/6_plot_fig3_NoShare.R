library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(patchwork)
library(facetscales)


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
  result <- mean(y, trim = 0.0)
  return(result)
}

# Analyze results files
analyze <- function(df, sigma, scenario, scenario_ups, t_ups, ntel){
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
      mutate(Upsilon = scenario_ups) %>%
      mutate(ntel = ntel)
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
      mutate(Upsilon = scenario_ups) %>%
      mutate(ntel = ntel)
  }
  return(result)
}


#----Analyze results -- Small ups----

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
  mutate(Upsilon = "high-res") %>%
  mutate(ntel = 0)

# NTEL = 1
ntel1_unshared <- read.table("./output/small ups/est_ntel=1_share=FALSE/results.txt") %>% 
  analyze(., sigma = "unshared", scenario = "ntel=1, unshared", scenario_ups = "high-res", 
          t_ups = 0.25*2.5, ntel = 1)
ntel1_shared <- read.table("./output/small ups/est_ntel=1_share=TRUE/results.txt") %>% 
  analyze(., sigma = "shared", scenario = "ntel=1, shared", scenario_ups = "high-res", 
          t_ups = 0.25*2.5, ntel = 1)

# NTEL = 3
ntel3_unshared <- read.table("./output/small ups/est_ntel=3_share=FALSE/results.txt") %>%
  analyze(., sigma = "unshared", scenario = "ntel=3, unshared", scenario_ups = "high-res", 
          t_ups = 0.25*2.5, ntel = 3)
ntel3_shared <- read.table("./output/small ups/est_ntel=3_share=TRUE/results.txt") %>% 
  analyze(., sigma = "shared", scenario = "ntel=3, shared", scenario_ups = "high-res", 
          t_ups = 0.25*2.5, ntel = 3)

# NTEL = 5
ntel5_unshared <- read.table("./output/small ups/est_ntel=5_share=FALSE/results.txt") %>% 
  analyze(., sigma = "unshared", scenario = "ntel=5, unshared", scenario_ups = "high-res", 
          t_ups = 0.25*2.5, ntel = 5)
ntel5_shared <- read.table("./output/small ups/est_ntel=5_share=TRUE/results.txt") %>% 
  analyze(., sigma = "shared", scenario = "ntel=5, shared", scenario_ups = "high-res", 
          t_ups = 0.25*2.5, ntel = 5)


#----Analyze results -- Big ups----

# No movement
bigups_noMove <- read.table("./output/big ups/est_noMove/results.txt") %>%
  as.data.frame() %>%
  filter(sig < 3) %>%
  select(-p0) %>%
  gather() %>%
  mutate(true = rep(c(t_alpha2, t_sig, t_d0), each = nrow(.)/3)) %>%
  mutate(prbias = 100*((value-true)/true)) %>%
  mutate(key = recode(key, alpha2 = "cost", d0 = "density", sig = "sigma")) %>%
  mutate(key = factor(key, levels = c("cost", "density", "sigma"))) %>%
  mutate(Scenario = "No movement") %>%
  mutate(Upsilon = "low-res") %>%
  mutate(ntel = 0)

# NTEL = 1
bigups_ntel1_unshared <- read.table("./output/big ups/est_ntel=1_share=FALSE/results.txt") %>% 
  analyze(., sigma = "unshared", scenario = "ntel=1, unshared", scenario_ups = "low-res", 
          t_ups = 1*2.5, ntel = 1)
bigups_ntel1_shared <- read.table("./output/big ups/est_ntel=1_share=TRUE/results.txt") %>% 
  analyze(., sigma = "shared", scenario = "ntel=1, shared", scenario_ups = "low-res", 
          t_ups = 1*2.5, ntel = 1)

# NTEL = 3
bigups_ntel3_unshared <- read.table("./output/big ups/est_ntel=3_share=FALSE/results_PRELIM.txt") %>%
  analyze(., sigma = "unshared", scenario = "ntel=3, unshared", scenario_ups = "low-res", 
          t_ups = 1*2.5, ntel = 3)
bigups_ntel3_shared <- read.table("./output/big ups/est_ntel=3_share=TRUE/results.txt") %>%
  analyze(., sigma = "shared", scenario = "ntel=3, shared", scenario_ups = "low-res",
          t_ups = 1*2.5, ntel = 3)

# NTEL = 5
bigups_ntel5_unshared <- read.table("./output/big ups/est_ntel=5_share=FALSE/results.txt") %>%
  analyze(., sigma = "unshared", scenario = "ntel=5, unshared", scenario_ups = "low-res",
          t_ups = 1*2.5, ntel = 5)
bigups_ntel5_shared <- read.table("./output/big ups/est_ntel=5_share=TRUE/results.txt") %>% 
  analyze(., sigma = "shared", scenario = "ntel=5, shared", scenario_ups = "low-res", 
          t_ups = 1*2.5, ntel = 5)


#----Combine data----

# missing <- expand.grid(
#   key = c("cost", "density", "sigma", "sigma_move"),
#   value = c(NA),
#   true = c(NA),
#   prbias = c(NA),
#   Scenario = c("ntel=3, shared", "ntel=5, unshared"),
#   Upsilon = "low-res",
#   ntel = c(3, 5)
# )

# Combine
df <- rbind(
  noMove,
  ntel1_unshared, ntel1_shared,
  ntel3_unshared, ntel3_shared,
  ntel5_unshared, ntel5_shared,
  bigups_noMove,
  bigups_ntel1_unshared, 
  bigups_ntel1_shared,
  bigups_ntel3_unshared, 
  bigups_ntel3_shared,
  bigups_ntel5_unshared,
  bigups_ntel5_shared) %>%
  filter(key != "upsilon") %>%
  filter(key != "pr(moved)") %>%
  #rbind(., missing) %>%
  mutate(Scenario = factor(Scenario, 
                           levels = c("No movement", 
                                      "ntel=1, shared", "ntel=1, unshared", 
                                      "ntel=3, shared", "ntel=3, unshared",  
                                      "ntel=5, shared", "ntel=5, unshared"))) %>%
  mutate(Upsilon = factor(Upsilon, levels = c("high-res", "low-res")))


#----Summarize results----

results0 <- df %>%
  group_by(key, Scenario, Upsilon, ntel) %>%
  #na.omit() %>%
  summarise(p.mean = mean(prbias), # we no longer use trim mean
            bias_upper = quantile(prbias, 0.75),
            bias_lower = quantile(prbias, 0.25))


results <- results0 %>%
  #rbind(., missing_results) %>%
  mutate(key = recode(key,
                      cost = "Cost ~ (alpha[1])",
                      density = "Density ~ (lambda)",
                      sigma = "sigma[SCR]",
                      sigma_move = "sigma[MM]")) %>%
  mutate(Upsilon = recode(Upsilon,
                          `high-res` = "sigma < sigma[det]",
                          `low-res` = "sigma == sigma[det]"))


#----NEW plot----

# Color palettes
ibm <- c("#648fff", "#785ef0", "#dc267f")
pal <- c("black", ibm[1], ibm[1], ibm[2], ibm[2], ibm[3], ibm[3])
pal <- c("black", ibm[3], ibm[3], ibm[3])

# Facet custom scales
scales_y <- list(
  "Cost ~ (alpha[1])" = scale_y_continuous(
    limits = c(-60,60), breaks = c(-50, -25, 0, 25, 50)),
  "Density ~ (lambda)" = scale_y_continuous(
    limits = c(-15,15), breaks = c(-10, 0, 10)))

# Color palettes
ibm <- c("#648fff", "#785ef0", "#dc267f")
pal <- ibm[c(1,3)]

fig_dat <- results %>%
  filter(Scenario %in% as.character(unique(results$Scenario)[c(1,3,5,7)])) %>%
  filter(key %in% c(
    expression(Cost ~ (alpha[1])),
    expression(Density ~ (lambda)))) %>%
  mutate(Upsilon = plyr::revalue(
    Upsilon, 
    c("sigma < sigma[det]" = "sigma[step] < sigma[home]", 
      "sigma == sigma[det]" = "sigma[step] == sigma[home]")))

p1 <- ggplot(data = fig_dat %>% filter(key %in% expression(Cost ~ (alpha[1]))), 
             aes(x = Scenario, group = Upsilon, shape = Scenario, color = Upsilon)) +
  # Bias bar
  geom_rect(ymin = -5, ymax = 5, xmin = -2, xmax = 10, 
            fill = "gray91", color = "gray91") +
  geom_hline(yintercept = 0, color = "gray40", size = 0.7) +
  
  # Main results
  geom_pointrange(size = 0.7, position = position_dodge(0.6),
                  aes(y = p.mean, ymin = bias_lower, ymax = bias_upper)) +
  # Color
  scale_color_manual(values = pal) +
  scale_shape_manual(values = c(17, 15, 15, 15)) +
  # Ntel text
  geom_text(aes(y = p.mean, label = ntel, x = Scenario, group = Upsilon), 
            color = "white", size = 2.5, fontface = "bold",
            position = position_dodge(0.6)) +
  # Facet
  ylim(-60,60) +
  # Theme
  labs(y = "% Relative bias", x = NULL, title = Cost ~ (alpha[1])) +
  theme_minimal() +
  theme(aspect.ratio = 1,
        legend.position = "none",
        text = element_text(size = 13),
        plot.title = element_text(hjust=0.5, size=13),
        strip.text.y = element_text(angle = 0, hjust = 0),
        panel.border = element_rect(fill = NA, color = "gray70", size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_blank())

p2 <- ggplot(data = fig_dat %>% filter(key %in% expression(Density ~ (lambda))), 
             aes(x = Scenario, group = Upsilon, shape = Scenario, color = Upsilon)) +
  # Bias bar
  geom_rect(ymin = -5, ymax = 5, xmin = -2, xmax = 10, 
            fill = "gray91", color = "gray91") +
  geom_hline(yintercept = 0, color = "gray40", size = 0.7) +
  
  # Main results
  geom_pointrange(size = 0.7, position = position_dodge(0.6),
                  aes(y = p.mean, ymin = bias_lower, ymax = bias_upper)) +
  # Color
  scale_color_manual(values = pal) +
  scale_shape_manual(values = c(17, 15, 15, 15)) +
  # Ntel text
  geom_text(aes(y = p.mean, label = ntel, x = Scenario, group = Upsilon), 
            color = "white", size = 2.5, fontface = "bold",
            position = position_dodge(0.6)) +
  # Facet
  ylim(-15,15) +
  # Theme
  labs(y = NULL, x = NULL, title = Density ~ (lambda)) +
  theme_minimal() +
  theme(aspect.ratio = 1,
        legend.position = "none",
        text = element_text(size = 13),
        plot.title = element_text(hjust=0.5, size=13),
        strip.text.y = element_text(angle = 0, hjust = 0),
        panel.border = element_rect(fill = NA, color = "gray70", size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_blank())

p <- p1+p2;p

ggsave(filename = "Figure2.pdf", plot = p, device="pdf",
       dpi = 600, scale = 0.8, height = 5.25, width = 8.75,
       path = "/Users/gatesdupont/Desktop")








parse.labels <- function(x) parse(text = x)

p1 <- fig_dat %>% 
  filter(key %in% expression(Cost ~ (alpha[1]))) %>%
  mutate(ntel = as.character(ntel)) %>%
  ggplot(data = ., aes(x = ntel, group = Upsilon, shape = ntel, color = Upsilon)) +
  # Bias bar
  geom_rect(ymin = -5, ymax = 5, xmin = -2, xmax = 10, 
            fill = "gray91", color = "gray91") +
  geom_hline(yintercept = 0, color = "gray40", size = 0.7) +
  
  # Main results
  geom_pointrange(size = 0.7, position = position_dodge(0.6),
                  aes(y = p.mean, ymin = bias_lower, ymax = bias_upper)) +
  # Color
  scale_color_manual(values = pal) +
  scale_shape_manual(values = c(17, 15, 15, 15)) +
  # Ntel text
  geom_text(aes(y = p.mean, label = ntel, x = ntel, group = Upsilon), 
            color = "white", size = 2.5, fontface = "bold",
            position = position_dodge(0.6)) +
  # Facet
  ylim(-60,60) +
  # Theme
  labs(y = "% Relative bias", x = NULL, title = Cost ~ (alpha[1])) +
  theme_minimal() +
  theme(aspect.ratio = 1,
        legend.position = "none",
        text = element_text(size = 13),
        plot.title = element_text(hjust=0.5, size=13),
        strip.text.y = element_text(angle = 0, hjust = 0),
        panel.border = element_rect(fill = NA, color = "gray70", size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_blank());p1

p2 <- fig_dat %>% 
  filter(key %in% expression(Density ~ (lambda))) %>%
  mutate(ntel = as.character(ntel)) %>%
  ggplot(data = ., aes(x = ntel, group = Upsilon, shape = ntel, color = Upsilon)) +
  # Bias bar
  geom_rect(ymin = -5, ymax = 5, xmin = -2, xmax = 10, 
            fill = "gray91", color = "gray91") +
  geom_hline(yintercept = 0, color = "gray40", size = 0.7) +
  
  # Main results
  geom_pointrange(size = 0.7, position = position_dodge(0.6),
                  aes(y = p.mean, ymin = bias_lower, ymax = bias_upper)) +
  # Color
  scale_color_manual(values = pal, NULL, labels = parse.labels) +
  scale_shape_manual(values = c(17, 15, 15, 15), guide=F) +
  # Ntel text
  geom_text(aes(y = p.mean, label = ntel, x = ntel, group = Upsilon), 
            color = "white", size = 2.5, fontface = "bold",
            position = position_dodge(0.6)) +
  guides(color = guide_legend(override.aes = list(shape = 15, linetype = 0, size = 1.1))) +
  # Facet
  ylim(-15,15) +
  # Theme
  labs(y = NULL, x = NULL, title = Density ~ (lambda)) +
  theme_minimal() +
  theme(aspect.ratio = 1,
        #legend.position = "none",
        text = element_text(size = 13),
        plot.title = element_text(hjust=0.5, size=13),
        strip.text.y = element_text(angle = 0, hjust = 0),
        panel.border = element_rect(fill = NA, color = "gray70", size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_blank());p2

p <- p1+p2;p

ggsave(filename = "Figure2.pdf", plot = p, device="pdf",
       dpi = 600, scale = 0.8, height = 3.75, width = 8.75,
       path = "/Users/gatesdupont/Desktop")
