library(oSCR)
library(dplyr)
library(ggplot2)
library(lubridate)
library(raster)
library(viridis)
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

# Clean up df
df <- df0 %>%
  mutate(fix_time = dmy_hms(paste(Date, Time))) %>%
  select(id = animalid, fixnum, lat = Lat, long = Long, time = fix_time)

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










