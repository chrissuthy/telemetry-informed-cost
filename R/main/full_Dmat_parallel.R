# This only takes 40 seconds! Worth adding to speed up the code!

df <- data.frame(
  slice1 = floor(seq(1, nrow(landscape), length.out = 11))[1:10],
  slice2 = floor(seq(1, nrow(landscape), length.out = 11))[2:11]
)

df[2:10, 1] <- df[2:10, 1]+1

dmat <- list()

# Parallel setup
ncores = detectCores()-1 # Number of available cores -1 to leave for computer
cl = makeCluster(ncores) # Make the cluster with that many cores
registerDoParallel(cl)  
clusterExport(cl, varlist = c("e2dist"), envir = environment()) # Export required function to the cores

t0 <- Sys.time()
dmat <- foreach(j=1:10, .packages = c(.packages())) %dopar% {

  m <- costDistance(
    tr1Corr,  
    from = landscape %>% select(x,y) %>% slice(df[j,1]:df[j,2]) %>% as.matrix(),
    to = landscape %>% select(x,y) %>% as.matrix())
  
}
stopCluster(cl)
tf <- Sys.time()

dmat <- do.call(rbind, dmat)
