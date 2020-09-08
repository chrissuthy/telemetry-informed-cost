scr_tic <- function(start, 
                    add_telem = FALSE,
                    scr_y, 
                    rsf_y,
                    K, 
                    trap_locs,
                    landscape, #is scr_statespace and rsf_statespace
                    direction=16, 
                    dist="ecol"){

    
    ## parameter set up
    p0 <- plogis(start[1])
    sigma <- exp(start[2])
    d0 <- exp(start[3])
    alpha_cost <- exp(start[4])
    ss_locs <- as.matrix(coordinates(landscape))
    n_pixels <- nrow(ss_locs)
    if(!is.matrix(trap_locs)) trap_locs <- as.matrix(as.data.frame(trap_locs))
  
    ## general cost surface items 
    cost_surface <- exp(alpha_cost * landscape)
    trLayer <- geoCorrection(
                             transition(cost_surface,
                                        transitionFunction = function(x) (1/(mean(x))),
                                        direction = direction),
                                        scl = F)


    ## SCR likelihood
    scr_dmat <- costDistance(trLayer, trap_locs, ss_locs)
    scr_pmat <- p0 * exp(-scr_dmat^2/(2*sigma^2))
    
    Pm <- matrix(NA, nrow = nrow(scr_pmat), ncol = ncol(scr_pmat))
    scr_y <- rbind(scr_y, rep(0, ncol(scr_y)))
    lik.marg <- rep(NA, nrow(ymat))
    
    for (i in 1:nrow(scr_y)){
        Pm[1:length(Pm)] <- (dbinom(x = rep(ymat[i, ], nG), 
                                    size = rep(K,nG),
                                    prob = probcap[1:length(Pm)], 
                                    log = TRUE))
        lik.cond <- exp(colSums(Pm))
        lik.marg[i] <- sum(lik.cond * (1/n_pixels))
    }
    

    nv <- c(rep(1, length(lik.marg) - 1), 1)
    atheta <- 1 - lik.marg[nrow(scr_y)]
    nind <- nrow(scr_y) - 1
    part1 <- nind * log(sum(d0 * n_pixels)) - sum(d0 * n_pixels) * atheta
    part2 <- sum(nv[1:nind] * log(lik.marg[1:nind]))
    scr_out <- part1 + part2
    
    ## RSF likelihood (DAN, you probably want to step in here!)
    if(add_telem){
      rsf_dmat <- costDistance(trLayer, ss_locs, ss_locs)
      rsf_pmat <- exp(1-rsf_dmat^2/(2*sigma^2))
      rsf_pmat <- rsf_pmat/apply(rsf_pmat,1,sum)
  
     for (i in 1:nrow(rsf_y)){
      
       log.probs <- log((1-2*1e-25)*rsf_pmat + 1e-25)
       lik.marg.tel[i] <- sum( exp(rsf_y[i,,drop=F] %*% log.probs) * as.vector(pi.s) )
       #if (telemetry.type == "dep"){
       #if (i <= length(cap.tel)){
       # combine conditional likelihoods if some collared ind were captured
       lik.cond.tot <- (rsf_y[i,,drop=F] %*% log.probs) + lik.cond.tel[i,]
       #lik.cond.tot[trimC[[s]][[cap.tel[i]]]] <- exp(lik.cond.tot[trimC[[s]][[cap.tel[i]]]])
       lik.cond.tot[is.na(lik.cond.tot)] <- 0
       lik.cond.tot[lik.cond.tot != 0] <- exp(lik.cond.tot[lik.cond.tot != 0])
                     
       lik.marg[cap.tel[i]] <- sum(lik.cond.tot * as.vector(pi.s)) 
       lik.marg.tel[i] <- 1
     }
    
    rsf_out <- sum(log(lik.marg.tel))
    
    }else{

      rsf_out <- 0
      
    }

    ## combine the ll components    
    tic_out <- -1*(scr_out + rsf_out)
}

