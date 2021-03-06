# # # # # # # # # # # # # # #
#                           #
#  ____HEY!!_DO THIS!!____  #
#     Parameterize sims     #
#                           #
# # # # # # # # # # # # # # #

# Right here:
select_ups  <- c("small ups", "big ups")[2]

# # # # # # # # # # # # # # #
#                           #
#       GREAT, THANKS!      #
#                           #
# # # # # # # # # # # # # # #

library(gdistance)
library(dplyr)
select = dplyr::select

file <- paste0("output/", select_ups, "/model_data/cost_data.RData")
spatdata_old  <- readRDS(file) # colnames? # fine


# Re-construct spatdata
spatdata <- list()
for(sim in 1:length(spatdata_old)){
  
  spatdata[[sim]] <- list()
  for(ind in 1:length(spatdata_old[[sim]])){
    
    tmp_df <- spatdata_old[[sim]][[ind]] %>%
      as.data.frame()
    
    sbar <- tmp_df %>%
      select(x, y) %>%
      colMeans() %>%
      as.numeric() %>%
      matrix(ncol = 2)
    
    tmp_r <- raster::rasterFromXYZ(tmp_df)
    
    sbar_indx <- raster::extract(x = tmp_r, y = sbar, cellnumbers=T)[,1]
    sbar_on_r <- tmp_df[sbar_indx,c("x", "y")]
    
    tmp_result <- tmp_df %>%
      select(x,y) %>%
      mutate(sbar = ifelse(
        (x == sbar_on_r[,1]) & (y == sbar_on_r[,2]),
        1,0)) %>%
      as.matrix()
    
    spatdata[[sim]][[ind]] <- tmp_result
    
  }
  
  file <- paste0("output/", select_ups, "/model_data/cost_data_light/", "cost_data_", sim, ".RData")
  saveRDS(spatdata[[sim]], file)
  
}
