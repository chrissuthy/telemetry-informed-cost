library(raster)
library(ggplot2)
library(viridis)
library(patchwork)
library(oSCR)
library(dplyr)
library(tidyr)
library(lubridate)
library(amt)

# Get data
data("nybears")

# Select telem
df0 <- nybears$teldata

# Plot the data
ggplot(data = df0, aes(x = Lat, y = Long)) +
  geom_path() +
  facet_wrap(~animalid, scales = "free") +
  theme_minimal() +
  theme(aspect.ratio = 1)

"BASIC EXPLORATORY"

# Clean up df
df <- df0 %>%
  mutate(fix_time = dmy_hms(paste(Date, Time))) %>%
  select(X_UTM, Y_UTM, id = animalid, fixnum, lat = Lat, long = Long, time = fix_time)

# Calculate the time difference between fixes per individual
diffs <- df %>%
  group_by(id) %>%
  summarise(diffs = as.numeric(diff(time), units = "mins")) %>%
  mutate(Index_ind = 1:n()) %>%
  ungroup() %>%
  mutate(Index_total = 1:nrow(.))

# Plot the diffs per individual
ggplot(data = diffs, aes(y = diffs, x = Index_ind)) +
  geom_point(aes(color = id)) +
  ylim(0, NA) +
  facet_wrap(~id) +
  theme_minimal() +
  theme(legend.position = "none")

# Record diffs into groups of sampling freqs
diffs.summary <- diffs %>%
  mutate(diffs.summary0 = if_else(diffs < 70 & diffs > 50, 0, NULL)) %>%
  mutate(diffs.summary1 = if_else(diffs < 130 & diffs > 110, 1, NULL)) %>%
  mutate(diffs.summary2 = if_else(diffs < 190 & diffs > 170, 2, NULL)) %>%
  mutate(diffs.summary3 = if_else(diffs > 190, 3, NULL)) %>%
  rowwise() %>%
  mutate(diffs.summary = sum(diffs.summary0, diffs.summary1, diffs.summary2, diffs.summary3, na.rm = T)) %>%
  select(id, diffs.summary)

# Create table of grouped sampling freqs
diffs.tab <- table(diffs.summary$id, diffs.summary$diffs.summary)

# Calculate the percent of fixes with no misses
apply(diffs.tab, 1, function(x) x[1]/(sum(x[1:4])))

# Do this row-wise for each cell
prop.table(table(diffs.summary$id, diffs.summary$diffs.summary), 1)


"SPLIT TRACKS"

# Loop through each individual, split tracks.
tracks_list <- list()
sl_list <- list()
ny_crs <- CRS(sf::st_crs("EPSG:32618")$proj4string)

ids <- df %>% pull(id) %>% unique %>% as.character
for(i in 1:length(ids)){
  
  # Split tracks
  tmp_df <- df %>% 
    filter(id == ids[i]) %>%
    arrange(time) %>%
    mk_track(., .x = X_UTM, .y = Y_UTM, .t = time, crs = ny_crs) %>%
    track_resample(., rate = minutes(60), tolerance = minutes(7)) %>%
    filter_min_n_burst(min_n = 2) %>%
    mutate(id = ids[i])
  
  # (Also) record step length
  tmp_sl <- tmp_df %>%
    select(-id) %>%
    steps_by_burst() %>%
    pull(sl_)

  tracks_list[[i]] <- tmp_df
  sl_list[[i]] <- tmp_sl
  
}

# Compile data from loops
tracks_df <- do.call(rbind, tracks_list)
sl_df <- data.frame(
  id = c(
    rep(ids[1], length(sl_list[[1]])),
    rep(ids[2], length(sl_list[[2]])),
    rep(ids[3], length(sl_list[[3]]))),
  sl = unlist(sl_list)
)


# Plot the tracks now that they're split
p1 <- ggplot(data = tracks_df, aes(x = x_, y = y_, group = burst_, color = burst_)) +
  geom_path() +
  scale_color_viridis("Track ID") +
  facet_wrap(~id, scales = "free") +
  labs(x = "X", y = "Y") +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(aspect.ratio = 1)
p1

# Histogram of step length
p2 <- ggplot(data = sl_df, aes(sl)) +
  geom_histogram() +
  facet_wrap(~id, scales = "free_y") +
  xlab("Step length (m)") +
  theme_minimal() +
  theme(aspect.ratio = 1)
p2

# Histogram of step length -- Zoomed
p3 <- ggplot(data = sl_df, aes(sl)) +
  geom_histogram() +
  facet_wrap(~id, scales = "free_y") +
  xlab("Step length (m)") +
  theme_minimal() +
  xlim(0,50) +
  theme(aspect.ratio = 1)
p3


"COMPARE STEP LENGTH TO RASTER"

# Create data frame of ss with elevation
r_df <- cbind(nybears$ss, nybears$elevation)
colnames(r_df) <- c("x", "y", "z")
r <- rasterFromXYZ(r_df)


# Make a separate plot for each individual because of geom_tile scale issue
plot_r_tracks <- function(df, id_x){
  
  tracks_df <- df
  
  a1 <- tracks_df %>%
    filter(id == ids[id_x])
  
  tmp_b <- ggplot(data = a1) +
    geom_tile(data=r_df, aes(x=x, y=y, fill=z)) +
    geom_path(aes(x=x_, y=y_, group = burst_), color = "black") +
    scale_fill_gradientn(colours=c("red", "yellow")) +
    xlim(min(a1$x_), max(a1$x_)) + 
    ylim(min(a1$y_), max(a1$y_)) +
    theme_minimal() + ggtitle(paste0("N fixes = ", nrow(a1))) +
    theme(aspect.ratio = 1, legend.position = "none")
  
  return(tmp_b)
}

b1 <- plot_r_tracks(tracks_df, 1)
b2 <- plot_r_tracks(tracks_df, 2)
b3 <- plot_r_tracks(tracks_df, 3)
  
# Final plot
(b1+b2+b3)/p2 / p3


"CHECK FOR CIRCADIAN STEP LENGTH"

# Function to split tracks and get hourly step length
df_sl_hr_all <- function(df, id_x){
  
  tmp_df <- df %>% 
    filter(id == ids[id_x]) %>%
    arrange(time) %>%
    mk_track(., .x = X_UTM, .y = Y_UTM, .t = time, crs = ny_crs) %>%
    track_resample(., rate = minutes(60), tolerance = minutes(7)) %>%
    filter_min_n_burst(min_n = 2) %>%
    steps_by_burst() %>%
    mutate(hr = hour(t1_), 
           id = ids[id_x]) %>%
    select(sl_, hr, id)
  
  return(tmp_df)
  
}

# Bring all the data together
sl_hr_all <- rbind(
  df_sl_hr_all(df, 1),
  df_sl_hr_all(df, 2),
  df_sl_hr_all(df, 3)
)

# Final plot
ggplot(data = sl_hr_all, aes(x = hr, y = sl_, color = id)) +
  geom_point(pch = 1, alpha = 0.5) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 4)) +
  scale_y_log10() +
  facet_grid(~id) +
  xlab("Hour") +
  ylab("Step length (m)") +
  theme_minimal() +
  theme(aspect.ratio = 1, legend.position = "none")

# Check step length for more active hours
sl_hr_all %>%
  filter( (hr<3 | hr>22) | (hr>11 | hr<14) ) %>%
  ggplot(data = ., aes(x = sl_, fill = id)) +
  geom_histogram() +
  facet_wrap(~id, scales = "free") +
  theme_minimal() +
  xlim(0,50) +
  theme(aspect.ratio = 1, legend.position = "none")


"AVG HOME RANGE RADIUS"

# Calculate sbar
s_bars <- df %>%
  group_by(id) %>%
  summarise(
    mean_x = mean(X_UTM),
    mean_y = mean(Y_UTM))

# Calculate distance from sbar for each fix
c_list1 <- list()
c_list2 <- list()
c_list3 <- list()
c_list <- list(c_list1, c_list2, c_list3)
for(i in 1:length(ids)){
  
  s_bar <- s_bars %>%
    filter(id == ids[i])
  
  tmp_df <- df %>%
    filter(id == ids[i]) %>%
    select(x = X_UTM, y = Y_UTM)
  
  s1 <- tmp_df
  s2 <- as.numeric(s_bar[,c(2,3)])
  
  
  for(j in 1:nrow(s1)){
    c_list[[i]][[j]] <- sqrt(sum(abs(s2-s1[j,])^2))
  }

}

# Plot the output
c <- unlist(c_list)
c_df <- data.frame(c)
ggplot(data = c_df, aes(x = c)) +
  geom_histogram(fill = "purple1", color = "white") +
  theme_minimal() + 
  xlab("Distance from sbar (m)") +
  theme(aspect.ratio = 1)
  
# Home range radius
ggplot(data = sl_df, aes(sl)) +
  geom_histogram() +
  #facet_wrap(~id, scales = "free_y") +
  xlab("Step length (m)") +
  theme_minimal() +
  xlim(101,NA) +
  theme(aspect.ratio = 1)


#----Approximating parameters----

cutoff <- 1000

"CALCULATE UPSILON"

ups_test <- sl_df %>% filter(sl > cutoff) %>% pull(sl)
ups_quant <- as.numeric(quantile(ups_test , prob = 0.95))
upsilon <- (ups_quant)/sqrt(5.99)
upsilon


"CALCULATE SIGMA"

sig_test <- c[c > cutoff]
sig_quant <- as.numeric(quantile(sig_test , prob = 0.95))
sigma <- (sig_quant)/sqrt(5.99)
sigma


"CALCULATE PSI"

psi_num <- sl_df %>% filter(sl > cutoff) %>% nrow()
psi_denom <- nrow(sl_df)
psi <- psi_num/psi_denom
psi



# Relationship btwn step length cutoff and psi
out_res <- seq(10, 1000, by = 10)
out_psi <- c()

for(i in 1:length(out_res)){
  
  psi_num <- sl_df %>% filter(sl > out_res[i]) %>% nrow()
  psi_denom <- nrow(sl_df)
  psi <- psi_num/psi_denom
  
  out_psi[i] <- psi
}


# Relationship btwn step length cutoff and upsilon
out_res <- seq(10, 1000, by = 10)
out_ups <- c()

for(i in 1:length(out_res)){
  
  ups_test <- sl_df %>% filter(sl > out_res[i]) %>% pull(sl)
  ups_quant <- as.numeric(quantile(ups_test , prob = 0.95))
  upsilon <- (ups_quant)/sqrt(5.99)
  
  out_ups[i] <- upsilon
}

# Relationship btwn step length cutoff and sigma
out_res <- seq(10, 1000, by = 10)
out_sig <- c()

for(i in 1:length(out_res)){
  
  sig_test <- c[c > out_res[i]]
  sig_quant <- as.numeric(quantile(sig_test , prob = 0.95))
  sigma <- (sig_quant)/sqrt(5.99)
  
  out_sig[i] <- sigma
}





par(mfrow=c(1,3))

par(pty="s")
plot(out_psi~out_res,
     type ="l", col = "red",
     xlab = "Step length cutoff (m)", 
     ylab = "Approximated psi")

par(pty="s")
plot(out_ups~out_res,
     type ="l", col = "red",
     xlab = "Step length cutoff (m)", 
     ylab = "Approximated upsilon", asp = 1)
abline(a = 0, b=1, lty = 2)


par(pty="s")
plot(out_sig~out_res,
     type ="l", col = "red",
     xlab = "Step length cutoff (m)", 
     ylab = "Approximated sigma")


par(mfrow=c(1,1))
plot(out_sig~out_ups,
     type ="l", col = "red",
     xlab = "Approximated upsilon", 
     ylab = "Approximated sigma")

