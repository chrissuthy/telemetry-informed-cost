library(dplyr)
library(ggplot2)

summary.SCR <- function(y_in){
  
  n_captured <- c()
  n_captured_2plus <- c()
  mean_captures <- c()
  
  for(i in 1:length(y_in)){
    y <- y_in[[i]]
    
    # Number of individuals captures
    n_captured[i] <- nrow(y)
    
    # Number of captures on more than 1 trap
    y_tmp <- apply(y, 1:2, sum)
    y_tmp[y_tmp > 0] <- 1
    n_captured_2plus[i] <- sum(rowSums(y_tmp)>1)
    
    # Mean number of captures per individual
    y_tmp <- apply(y, 1:2, sum)
    mean_captures[i] <- mean(rowSums(y_tmp))
  }
  
  return(list(n_captured, n_captured_2plus, mean_captures))
}

# Analyze small-upsilon scenario
file <- paste0("./output/", "small ups", "/model_data/y.RData")
y1 <- readRDS(file)
y1_summary <- summary.SCR(y1)

# Analyze small-upsilon scenario
file <- paste0("./output/", "big ups", "/model_data/y.RData")
y2 <- readRDS(file)
y2_summary <- summary.SCR(y1)


# Compile data
metrics <- c("n[captured]", "n[captured~`>1`~trap]", "mean[captures]")
y_summary <- c(
  unlist(y1_summary),
  unlist(y2_summary)) %>%
  data.frame(value = .) %>%
  mutate(key = rep(c(metrics, metrics), each = 100)) %>%
  mutate(scenario = rep(
    c("sigma < sigma[det]", "sigma == sigma[det]"), 
    each = 300)) %>%
  mutate(key = factor(key, levels = metrics))


# Plot
ggplot(data = y_summary, aes(x = value, fill = key)) +
  geom_histogram() +
  facet_grid(scenario~key, labeller = label_parsed, scales = "free") +
  labs(y="Frequency", x=NULL) +
  scale_fill_manual(values=c("#003f5c","#bc5090","#ffa600")) +
  theme_minimal() +
  theme(aspect.ratio = 1, 
        text = element_text(size=16),
        axis.text.x = element_text(size=10),
        legend.position = "none",
        panel.border = element_rect(fill = NA),
        panel.grid = element_blank())



