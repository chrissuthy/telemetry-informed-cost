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
tracks_list <- list()
sl_list <- list()

ids <- df %>% pull(id) %>% unique %>% as.character
for(i in 1:length(ids)){
  
  tmp_df <- df %>% 
    filter(id == ids[i]) %>%
    arrange(time) %>%
    mk_track(., .x = X_UTM, .y = Y_UTM, .t = time) %>%
    track_resample(., rate = minutes(60), tolerance = minutes(7)) %>%
    filter_min_n_burst(min_n = 2) %>%
    mutate(id = ids[i])
  
  tmp_sl <- tmp_df %>%
    select(-id) %>%
    steps_by_burst() %>%
    pull(sl_)
  
  tracks_list[[i]] <- tmp_df
  sl_list[[i]] <- tmp_sl
  
}

tracks_df <- do.call(rbind, tracks_list)
sl_df <- data.frame(
  id = c(
    rep(ids[1], length(sl_list[[1]])),
    rep(ids[2], length(sl_list[[2]])),
    rep(ids[3], length(sl_list[[3]]))),
  sl = unlist(sl_list)
)

p1 <- ggplot(data = tracks_df, aes(x = x_, y = y_, group = burst_, color = burst_)) +
  geom_path() +
  scale_color_viridis("Track ID") +
  facet_wrap(~id, scales = "free") +
  labs(x = "X", y = "Y") +
  theme_minimal() +
  theme(aspect.ratio = 1, legend.position = "none", axis.text = element_blank())

p2 <- ggplot(data = sl_df, aes(sl)) +
  geom_histogram() +
  facet_wrap(~id, scales = "free_y") +
  xlab("Step length") +
  theme_minimal() +
  theme(aspect.ratio = 1)

p3 <- ggplot(data = sl_df, aes(sl)) +
  geom_histogram() +
  facet_wrap(~id, scales = "free_y") +
  xlab("Step length") +
  theme_minimal() +
  xlim(0,50) +
  theme(aspect.ratio = 1)

p1 / p2 / p3












