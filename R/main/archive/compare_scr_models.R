# Gates Dupont
# October 2020

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

#----Setup----

# Density
density <- 300/6400

# Tests
true_alpha2 <- rep(c(0.1, 1:5), each = 50)


#----Load sim data----

# Load SCR(-T)
attach('output/Valpha2_scrNOtelem_10032020.RData')
out_scr <- out %>%
  as.data.frame() %>%
  mutate(true = true_alpha2) %>%
  mutate(model = "scr-telem")
detach('file:output/Valpha2_scrNOtelem_10032020.RData')

# Load SCR(+T)
attach('output/Valpha2_test1_10012020.RData')
out_scrT <- out %>%
  as.data.frame() %>%
  mutate(true = true_alpha2) %>%
  mutate(model = "scr+telem")
detach('file:output/Valpha2_test1_10012020.RData')


#----Combine----

# Select common and combine
out <- rbind(
  out_scr %>%
    select(alpha2, d0, true, model),
  out_scrT %>%
    select(alpha2, d0, true, model)
)


#----Relative bias----

bias_df <- out %>% 
  group_by(true, model) %>%
  mutate(
    b_cost = 100*(alpha2 - true),
    b_density = 100*(d0 - density),
    rb_cost = 100*((alpha2 - true)/true),
    rb_density = 100*((d0 - density)/density))


"SCR-TELEM"

#----Plot: Bias in COST----

# Color palettes
blue = brewer.pal(n = 9, name = "Blues")
red = brewer.pal(n = 9, name = "Reds")
pal1 = c(red[5:9])
pal2 = c(red[4:8])

# % Bias cost
bias_df %>%
  filter(true != 0.1) %>%
  filter(model == "scr-telem") %>%
  ggplot(data = ., aes(x = as.factor(true), 
                       y = b_cost,
                       color = as.factor(true),
                       fill = as.factor(true),
                       group = as.factor(true))) +
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  geom_hline(yintercept = c(-5,5), lty = 2, size = 0.2) +
  scale_fill_manual(expression(paste("\u03b1"[2])), values = pal2) +
  scale_color_manual(expression(paste("\u03b1"[2])), values = pal1) +
  theme_minimal() + ylim(-320, 320) +
  labs(x = expression(paste("\u03b1"[2])), 
       y = expression(paste("% Bias in ", hat("\u03b1"[2]))),
       title = "Cost estimates | 50 sims") +
  theme(aspect.ratio = 1)


#----Plot: Bias in DENSITY----

# Color palettes
blue = brewer.pal(n = 9, name = "Blues")
red = brewer.pal(n = 9, name = "Reds")
pal1 = c(red[5:9])
pal2 = c(red[4:8])

# % Bias cost
bias_df %>%
  filter(true != 0.1) %>%
  filter(model == "scr-telem") %>%
  ggplot(data = ., aes(x = as.factor(true), 
                       y = b_density,
                       color = as.factor(true),
                       fill = as.factor(true),
                       group = as.factor(true))) +
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  geom_hline(yintercept = c(-5,5), lty = 2, size = 0.2) +
  scale_fill_manual(expression(paste("\u03b1"[2])), values = pal2) +
  scale_color_manual(expression(paste("\u03b1"[2])), values = pal1) +
  theme_minimal() + ylim(-10, 10) +
  labs(x = expression(paste("\u03b1"[2])), 
       y = expression(paste("% Bias in ", hat("d"[0]))),
       title = "Density estimates | 50 sims") +
  theme(aspect.ratio = 1)




"SCR+TELEM"

#----Plot: Bias in COST----

# Color palettes
blue = brewer.pal(n = 9, name = "Blues")
red = brewer.pal(n = 9, name = "Reds")
pal1 = c(blue[5:9])
pal2 = c(blue[4:8])

# % Bias cost
bias_df %>%
  filter(true != 0.1) %>%
  filter(model == "scr+telem") %>%
  ggplot(data = ., aes(x = as.factor(true), 
                       y = b_cost,
                       color = as.factor(true),
                       fill = as.factor(true),
                       group = as.factor(true))) +
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  geom_hline(yintercept = c(-5,5), lty = 2, size = 0.2) +
  scale_fill_manual(expression(paste("\u03b1"[2])), values = pal2) +
  scale_color_manual(expression(paste("\u03b1"[2])), values = pal1) +
  theme_minimal() + ylim(-25, 25) +
  labs(x = expression(paste("\u03b1"[2])), 
       y = expression(paste("% Bias in ", hat("\u03b1"[2]))),
       title = "Cost estimates | 50 sims") +
  theme(aspect.ratio = 1)


#----Plot: Bias in DENSITY----

# Color palettes
blue = brewer.pal(n = 9, name = "Blues")
red = brewer.pal(n = 9, name = "Reds")
pal1 = c(blue[5:9])
pal2 = c(blue[4:8])

# % Bias cost
bias_df %>%
  filter(true != 0.1) %>%
  filter(model == "scr+telem") %>%
  ggplot(data = ., aes(x = as.factor(true), 
                       y = b_density,
                       color = as.factor(true),
                       fill = as.factor(true),
                       group = as.factor(true))) +
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  geom_hline(yintercept = c(-5,5), lty = 2, size = 0.2) +
  scale_fill_manual(expression(paste("\u03b1"[2])), values = pal2) +
  scale_color_manual(expression(paste("\u03b1"[2])), values = pal1) +
  theme_minimal() + ylim(-10, 10) +
  labs(x = expression(paste("\u03b1"[2])), 
       y = expression(paste("% Bias in ", hat("d"[0]))),
       title = "Density estimates | 50 sims") +
  theme(aspect.ratio = 1)


"COMPARE: COST"

#----Plot: COMPARE, Bias in COST----

# Color palettes
blue = brewer.pal(n = 9, name = "Blues")
red = brewer.pal(n = 9, name = "Reds")
pal1 = c(red[5:9], blue[5:9])
pal2 = c(red[4:8], blue[4:8])

# % Bias cost
bias_df %>%
  filter(true != 0.1) %>%
  ggplot(data = ., aes(x = as.factor(true), y = b_cost,
                       color = interaction(true, model, sep = ", "),
                       fill = interaction(true, model, sep = ", "),
                       group = interaction(true, model, sep = ", "))) +
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  geom_hline(yintercept = c(-5,5), lty = 2, size = 0.2) +
  scale_fill_manual(expression(paste("\u03b1"[2], ", model")), values = pal2) +
  scale_color_manual(expression(paste("\u03b1"[2], ", model")), values = pal1) +
  theme_minimal() +
  labs(x = expression(paste("\u03b1"[2])), 
       y = expression(paste("% Bias in ", hat("\u03b1"[2]))),
       title = "Cost estimates | 50 sims") +
  theme(aspect.ratio = 1)



"COMPARE: DENSITY"

#----Plot: COMPARE, Bias in DENSITY----

# Color palettes
blue = brewer.pal(n = 9, name = "Blues")
red = brewer.pal(n = 9, name = "Reds")
pal1 = c(red[5:9], blue[5:9])
pal2 = c(red[4:8], blue[4:8])

# % Bias cost
bias_df %>%
  filter(true != 0.1) %>%
  ggplot(data = ., aes(x = as.factor(true), y = b_density,
                       color = interaction(true, model, sep = ", "),
                       fill = interaction(true, model, sep = ", "),
                       group = interaction(true, model, sep = ", "))) +
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  geom_hline(yintercept = c(-5,5), lty = 2, size = 0.2) +
  scale_fill_manual(expression(paste("\u03b1"[2], ", model")), values = pal2) +
  scale_color_manual(expression(paste("\u03b1"[2], ", model")), values = pal1) +
  theme_minimal() + ylim(-8, 8) +
  labs(x = expression(paste("\u03b1"[2])), 
       y = expression(paste("% Bias in ", hat("d"[0]))),
       title = "Density estimates | 50 sims") +
  theme(aspect.ratio = 1)