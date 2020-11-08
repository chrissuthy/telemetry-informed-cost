

sims = 500

c1 <- sample(1:nrow(landscape), size = sims)
c2 <- sample(1:nrow(landscape), size = sims)

out <- c()
for(i in 1:sims){
  cat(i)
  
  a = costDistance(
    tr1Corr, 
    as.matrix(landscape[c1[i],1:2]), 
    as.matrix(landscape[c2[i],1:2])) %>%
    as.numeric()
  
  b = ref.Dmat_i(
    from = as.matrix(landscape[c1[i],1:2]),
    to = as.matrix(landscape[c2[i],1:2]),
    local_ss_r = landscape_r,
    Dmat_i = dmat) %>%
    as.numeric()
  
  out[i] = (a==b)
}

sum(out)
