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
  result <- mean(y, trim = 0.1)
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
  #filter(key != "upsilon") %>%
  #filter(key != "pr(moved)") %>%
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
  summarise(trim.mean = mean(prbias, trim = 0.1),
            bias_upper = quantile(prbias, 0.975),
            bias_lower = quantile(prbias, 0.025))

# missing_results <- expand.grid(
#   key = c("cost", "density", "sigma", "sigma_move"),
#   Scenario = c("ntel=3, shared", "ntel=5, unshared"),
#   Upsilon = "low-res",
#   ntel = c(3, 5),
#   trim.mean = NA,
#   bias_upper = NA,
#   bias_lower = NA
# )


results <- results0 %>%
  #rbind(., missing_results) %>%
  mutate(key = recode(key,
                      cost = "Cost ~ (alpha)",
                      density = "Density ~ (D)",
                      sigma = "sigma[SCR]",
                      sigma_move = "sigma[MM]")) %>%
  mutate(Upsilon = recode(Upsilon,
                          `high-res` = "upsilon < sigma",
                          `low-res` = "upsilon == sigma"))



#-----Table----

AppendixS1 <- results %>%
  select(key, Scenario, Upsilon, ntel, trim.mean, bias_lower, bias_upper) %>%
  mutate_if(is.numeric, round, 2) %>%
  mutate(LU = paste0("(",
    sprintf("%.2f", bias_lower), ", ",
    sprintf("%.2f", bias_upper), ")")) %>%
  select(-bias_lower, -bias_upper) %>%
  arrange(Upsilon) %>%
  mutate(key = recode(key,
                      `Cost ~ (alpha)` = "alpha_1",
                      `Density ~ (D)` = "d_0",
                      `sigma[SCR]` = "sigma_{det}'",
                      `sigma[MM]` = "sigma_{det}''",
                      `upsilon` = "sigma",
                      `pr(moved)` = "psi")) %>%
  mutate(Scenario = recode(Scenario,
                           `No movement` = "standard",
                           `ntel=1, shared` = "shared",
                           `ntel=1, unshared` = "independent",
                           `ntel=3, shared` = "shared",
                           `ntel=3, unshared` = "independent",
                           `ntel=5, shared` = "shared",
                           `ntel=5, unshared` = "independent")) %>%
  mutate(formulation = Scenario) %>%
  mutate(scenario = Upsilon) %>%
  mutate(scenario = recode(scenario,
                           `upsilon < sigma` = "small-sigma_{det}",
                           `upsilon == sigma` = "large-sigma_{det}")) %>%
  ungroup() %>%
  select(-Scenario, -Upsilon) %>%
  select(scenario, 
         parameter = key, 
         `sigma_{det}` = formulation, 
         `n_{tel}` = ntel, 
         mean = trim.mean, 
         bounds = LU)

write.csv(AppendixS1, file = "/Users/gatesdupont/Desktop/AppendixS1.csv")




