scr_cost_like <- function(
  param, dist=c("euc","circ","lcp")[3], 
  scr_y = NULL, K = NULL, trap_locs = NULL, landscape = NULL,
  mod=c("exp","gauss")[2], prj = NULL){
  
  # Debugging
  # browser()
  
  #alpha2: cost parameter
  #sigma: spatial scale (range)
  #p0: baseline encounter probability
  #d0: population density
  
  #if(is.null(prj)) stop("Must provide a projection")

  np <- 1
  
  # Old param tracking stuff
  {
    
    # This should really be limited to...
    # - fixcost == F
    # - popcost == T
    # - popmove == T
    # - use.sbar == T
    
    # Just commenting it out, shifted to below section
    # Kept in case we want to apply it elsewhere in the future
    # if(!fixcost){
    #   if(popcost){
    #     alpha2 <- matrix(param[1:np], nrow = nguys, ncol = np, byrow = TRUE)
    #     upsilon <- rep(exp(param[np+1]), nguys)
    #     if(popmove){
    #       psi <- rep(plogis(param[np+2]), nguys)
    #     }else{
    #       #psi <- plogis(param[(np+2):(np+2+nguys)])
    #     }
    #     if(use.sbar){
    #       sigma <- exp(param[length(param)])
    #     }
    #   }else{
    #     #alpha2 <- matrix(param[1:(nguys*np)], nrow = nguys, ncol = np)
    #     #upsilon <-  exp(param[(nguys*np+1):(nguys*np+nguys)])
    #     if(popmove){
    #       #psi <- rep(plogis(param[(nguys + nguys*np)+1]), nguys)
    #     }else{
    #       #psi <- plogis(param[(nguys+nguys*np+1):(nguys+nguys*np+nguys)] )
    #     }
    #     if(use.sbar){
    #       #sigma <- exp(param[length(param)])
    #     }
    #   }
    # }else{
    #   if(use.sbar){
    #     #sigma <- exp(param[length(param)]) 
    #   }
    #   #alpha2 <- matrix(rep(0, np))
    #   #upsilon <- exp(param[1])
    #   #psi <- plogis(param[2])
    # }
    
  }
  
  alpha2 <- matrix(param[1], nrow = 1, ncol = np, byrow = TRUE) # need to fix this for SCR
  sigma <- exp(param[2])
  
  ################################################## S C R ##################################################
  
  # SCR starting parameters
  a2_scr <- as.numeric(unique(alpha2))
  p0 <- plogis(param[3]) # Probability, on link scale
  d0 <- exp(param[4]) # Shouldn't be negative, so on link scale
  
  # Expanding K for later use
  if(length(K)==1) K<- rep(K,nrow(trap_locs)) # Generalized for irregular sampling periods
  
  G <- coordinates(landscape) # Pixels
  nG <- nrow(G) # Number of pixels
  
  # Cost distance pieces
  cost <- exp(a2_scr*landscape) # Cost surface w/ proposed parameter
  tr1 <- transition(cost,transitionFunction=function(x) 1/mean(x),directions=8)
  tr1CorrC <- geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
  D <- costDistance(tr1CorrC,trap_locs,G) # Cost distance
  
  # Half-normal encounter model
  probcap <- p0 * exp(- D^2/(2*sigma^2))  # rows = traps, col = pixels
  
  # Detection matrix, rows = traps, cols = pixels
  Pm <- matrix(NA,nrow=nrow(probcap),ncol=ncol(probcap))
  
  # Encounter historiess, augemented with 0 row
  # to estimate probability of an uncaptured individual
  # being in one of those pixels
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
  
  nll <- -1*scr_out
  
}