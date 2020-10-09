library(oSCR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(raster)
library(viridis)
library(moveHMM)
select = dplyr::select

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

# Add elevation to ss
elev <- data.frame(
  x = nybears$ss$X_UTM, 
  y = nybears$ss$Y_UTM, 
  elev = nybears$elevation)

# Plot tracks on raster just to get a sense for resolution
par(mfrow=c(1,3), oma=c(0,0,0,0),mar=c(1,1,1,1))
ids <- unique(df0$animalid)
for(i in 1:3){
  tmp_df <- df0 %>%
    filter(animalid == ids[i]) %>%
    select(x = X_UTM, y = Y_UTM)
  
  plot(y~x, data = tmp_df, asp = 1)
  plot(rasterFromXYZ(elev), viridis(1000), add = T, legend = F)
  lines(y~x, data = tmp_df)
}

dev.off()


"STEP LENGTH"

# Calculate step ratings saving all data
sldf <- df %>%
  group_by(id) %>%
  mutate(diffs = c(as.numeric(diff(time), units = "mins"), NA)) %>%
  mutate(Index_ind = 1:n()) %>%
  ungroup() %>%
  mutate(Index_total = 1:nrow(.)) %>%
  mutate(diffs.summary0 = if_else(diffs < 70 & diffs > 50, 0, NULL)) %>%
  mutate(diffs.summary1 = if_else(diffs < 130 & diffs > 110, 1, NULL)) %>%
  mutate(diffs.summary2 = if_else(diffs < 190 & diffs > 170, 2, NULL)) %>%
  mutate(diffs.summary3 = if_else(diffs > 190, 3, NULL)) %>%
  rowwise() %>%
  mutate(diffs.summary = sum(
    diffs.summary0, diffs.summary1, diffs.summary2, diffs.summary3, na.rm = T)) %>%
  select(X_UTM, Y_UTM, id, fixnum, lat, long, time, diffs, diffs.summary) %>%
  as.data.frame()

# Plot tracks with step ratings
sldf %>%
  ggplot(data = ., aes(x = long, y = lat, color = as.factor(diffs.summary))) +
  geom_path() +
  facet_wrap(~id, scales = "free") +
  labs(color = "Missed fixes") +
  theme_minimal() +
  theme(aspect.ratio = 1)

# Parse contiguous tracks
sldf_trackid <- sldf %>%
  arrange(id, time) %>%
  group_by(id) %>%
  mutate(track_id = if_else(diffs.summary == 0, 0, 1))

id <- 1
j = 0
for(i in 1:nrow(sldf_trackid)){
  
  j = j + 1
  
  if(as.data.frame(sldf_trackid)[i,"track_id"] == 0){
    
    sldf_trackid[i,"track_id"] <- id
    
  }else{
    
    # Break the track, start new
    sldf_trackid[i,"track_id"] <- -1
    id <- id + 1
    
  }
  
} # NOTE: THIS LOSES A STEP ON THE BACK END OF EACH TRACK ID

bears <- as.character(unique(df$id))
split_tracks <- list()

# -
for(b in 1:length(bears)){
  
  tmp_df_b <- sldf_trackid %>%
    as.data.frame() %>%
    filter(id == bears[b]) %>%
    arrange(time)
  
  # - -
  for(i in 2:(nrow(tmp_df_b))){
    
    # - - -
    tmp_df <- as.data.frame(tmp_df_b)
    
    # - - - If previous fix is part of a track and next fix is is -1, then...
    if(tmp_df[i-1,"track_id"] != -1 & tmp_df[i,"track_id"] == -1){
      
      # - - - - Add that final piece
      tmp_df_b[i,"track_id"] <- tmp_df[i-1,"track_id"]
    }
    
  }
  
  # - -
  split_tracks[[b]] <- tmp_df_b
  
}

split_tracks_df <- do.call(rbind, split_tracks)

sldf_trackid <- split_tracks_df %>%
  filter(track_id != -1) %>%
  as.data.frame() %>%
  arrange(track_id, time)

# Plot contiguous tracks
sldf_trackid %>%
  ggplot(data = ., aes(x = long, y = lat, group = track_id, color = track_id)) +
  scale_color_viridis() +
  geom_path(lwd=0.5) +
  theme_minimal() +
  facet_wrap(~id, scales = "free") +
  theme(aspect.ratio = 1)

# Calc step length function
dist <- function(df){
  result <- prepData(
    df %>% select(ID = id, long,lat),
    type="LL",coordNames=c("long","lat"))
  
  return(result$step)
}

# Calculate step length across tracks
trx <- unique(sldf_trackid$track_id)
trx_upsilon <- list()
for(i in 1:length(trx)){
  
  tmp_df <- sldf_trackid %>%
    filter(track_id == trx[i]) %>%
    select(ID = id, long, lat) %>%
    as.data.frame()
  
  result <- prepData(tmp_df, type="LL", coordNames = c("long","lat"))
  
  trx_upsilon[[i]] <- result$step
}

upsilon <- unlist(trx_upsilon)



data.frame()

prepData(data.frame(x = 1:3, y = c(1,2,2)), coordNames = c("x","y"))




dist <- function(df){
  result <- prepData(
    df %>% select(ID = id, X_UTM,Y_UTM),
    type="UTM",coordNames=c("X_UTM","Y_UTM"))
  
  return(result$step)
}

x <- dist(sldf_trackid)

upsilon <- round(mean(x[x<60], na.rm=T))

r_df <- cbind(nybears$ss, z= nybears$elevation) %>% as.data.frame()
r <- rasterFromXYZ(r_df)

ggplot() +
  geom_tile(data = r_df, aes(x=X_UTM, y=Y_UTM, fill = z)) +
  geom_path(data = sldf_trackid, 
            aes(
              x = X_UTM, y = Y_UTM, 
              group = track_id, color = as.factor(track_id))) +
  coord_equal() +
  theme(legend.position = "none")


#----Lat/Long----

dist <- function(df){
  result <- prepData(
    df %>% select(ID = id, long,lat),
    type="LL",coordNames=c("long","lat"))
  
  return(result$step)
}

x <- dist(sldf_trackid)

hist(x[x<0.05], breaks = 40)

mean(na.rm = T, x[x<0.05])

