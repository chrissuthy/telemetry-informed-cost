scr_move_cost_like <- function(
  param, share_sigma = TRUE, 
  dist=c("circ","lcp")[3], mod=c("exp","gauss")[2],
  teldata, spatdata = NULL, landscape = NULL, use.sbar=FALSE, 
  scr_y = NULL, K = NULL, trap_locs = NULL, scr_ss = NULL, 
  prj = NULL){

  #alpha2: cost parameter
  #upsilon: spatial scale (steps)
  #psi: pr(move cell)
  #sigma: spatial scale (range)
  #p0: baseline encounter probability
  #d0: population density
  
  # Basics
  nguys <- length(teldata)
  ll <- rep(NULL,nguys)
  np <- ncol(spatdata[[1]])-2
  
  # Some starting parameters
  alpha2 <- matrix(param[1], nrow = nguys, ncol = np, byrow = TRUE)
  upsilon <- rep(exp(param[2]), nguys)
  psi <- rep(plogis(param[3]), nguys)
  sigma <- exp(param[4])
  
  #----S C R----
  
  # SCR starting parameters
  a2_scr <- as.numeric(unique(alpha2))
  p0 <- plogis(param[5]) # Probability, on link scale
  d0 <- exp(param[6]) # Shouldn't be negative, so on link scale
  
  # Expanding K for later use
  if(length(K)==1) K<- rep(K,nrow(trap_locs)) # Generalized for irregular sampling periods
  
  G <- coordinates(scr_ss) # Pixels
  nG <- nrow(G) # Number of pixels
  
  # Cost distance pieces
  cost <- exp(a2_scr*landscape) # Cost surface w/ proposed parameter
  tr1 <- transition(cost,transitionFunction=function(x) 1/mean(x),directions=16)
  tr1CorrC <- geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
  D <- costDistance(tr1CorrC,trap_locs,G) # Cost distance
  
  # Half-normal encounter model
  probcap <- p0 * exp(- D^2/(2*sigma^2))  # rows = traps, col = pixels
  
  # Detection matrix, rows = traps, cols = pixels
  Pm <- matrix(NA,nrow=nrow(probcap),ncol=ncol(probcap))
  
  # Encounter historiess, augemented with 0 row
  # to estimate probability of an uncaptured individual
  # being in one of those pixels
  if(!is.na(dim(scr_y)[3])){
    scr_y <- apply(scr_y, 1:2, sum)
  }
  
  ymat <- rbind(scr_y,rep(0,ncol(scr_y)))
  
  # Loop through encounter histories (inidividuals + extra)
  # Calculate likelihood of that individuals activity center being at each pixel
  lik.marg <- rep(NA,nrow(ymat))
  for(i in 1:nrow(ymat)){
    Pm[1:length(Pm)] <-  (dbinom(rep(ymat[i,],nG),rep(K,nG),probcap[1:length(Pm)],log=TRUE))
    lik.cond <-  exp(colSums(Pm)) # Likelihood of number of captures conditional on s_i at pixel G_i
    lik.marg[i] <-  sum( lik.cond*(1/nG) ) # Likelihood averaged across all pixels
  }
  
  # Remaining likelihood pieces
  nv <- c(rep(1, length(lik.marg) - 1), 1)
  atheta <- 1 - lik.marg[nrow(ymat)]
  nind <- nrow(ymat) - 1
  part1 <- nind * log(sum(d0 * nG)) - sum(d0 * nG) * atheta
  part2 <- sum(nv[1:nind] * log(lik.marg[1:nind]))
  scr_out <- part1 + part2
  
  #---- Movement model ----
  
  sigma_mm <- ifelse(share_sigma == TRUE, sigma, exp(param[7]))
  
  #loop over individuals
  ll <- numeric(nguys)
  for(ind in 1:nguys){
    
    ss <- spatdata[[ind]][,1:2] # SPATDATA AS COORDINATES OF RASTER (still the bigger box)
    
    tel.locs <- as.matrix(teldata[[ind]])
    pixels <- raster::extract(rasterFromXYZ(cbind(ss, z=1)), tel.locs, cellnumbers=T)[,1]
    moved <- 1 - as.numeric(diff(pixels) == 0)
    
    if(dist=="circ"){
      D <- as.matrix(gdistance::commuteDistance(tr1CorrC, ss[1:(nrow(ss)-1),]))
      D <- D/ncell(cost)
      D <- D[pixels,]
    }
    if(dist=="lcp"){
      D_ss <- gdistance::costDistance(tr1CorrC,  ss, ss)
      D <- D_ss[pixels,]
    }
    
    if(use.sbar){
      
      dsbar <- D_ss[which(spatdata[[ind]][,3] == 1),]
      
      dsbar <- matrix(dsbar, nrow=nrow(D), ncol=length(dsbar), byrow=TRUE)
      
      if(mod=="gauss"){
        kern<- exp(-D*D/(2*upsilon[ind]*upsilon[ind]) - dsbar*dsbar/(2*sigma_mm*sigma_mm))
      }
      if(mod=="exp"){ 
        kern<- exp(-D/(2*upsilon[ind]*upsilon[ind]) - dsbar*dsbar/(2*sigma_mm*sigma_mm))
      }
    }else{
      if(mod=="gauss"){
        kern <- exp(-D*D/(2*upsilon[ind]*upsilon[ind]))
      }
      if(mod=="exp"){   
        kern <- exp(-D/(2*upsilon[ind]*upsilon[ind]))
      }
    }
    kern[cbind(1:nrow(kern), pixels )] <- 0 # *cannot move to same pixel bc CONDITIONAL on moved
    kern <- kern/rowSums(kern)
    probs <- kern[cbind(1:(nrow(kern)-1), pixels[-1] )] # remove terminus because it has no movement outcome 
    part1 <- moved*log(psi[ind]) + (1-moved)*log(1-psi[ind])
    part2 <- log(probs[moved==1])
    ll[ind]<- sum(part1) + sum(part2) 
  }
  nll <- -1*(sum(ll)+scr_out)
  nll
}