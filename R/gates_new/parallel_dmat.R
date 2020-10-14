# Parallel setup
ncores = detectCores()-1 # Number of available cores -1 to leave for computer
cl = makeCluster(ncores) # Make the cluster with that many cores
registerDoParallel(cl)  
clusterExport(cl, varlist = c("e2dist"), envir = environment()) # Export required function to the cores

# Data-collection list
XLdmat_pieces = list()


# Cluster!
t1 <- Sys.time()
XLdmat_pieces <- foreach(i=1:3, .packages = c(.packages())) %dopar% {
 
  chunk_size = c(1, 5720, 11440, 17161)
  
  tmp_m <- landscape %>% 
    select(x,y) %>%
    slice(chunk_size[i]:chunk_size[i+1]) %>%
    as.matrix
  
  XLdmat_i <- costDistance(
    tr1Corr,  
    fromCoords = tmp_m,
    toCoords = landscape[,c("x","y")] %>% as.matrix())
  
  XLdmat_pieces[[i]] <- XLdmat_i
   
}
t2 <- Sys.time()

t2-t1



