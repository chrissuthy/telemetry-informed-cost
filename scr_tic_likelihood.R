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
    if(!is.matrix(trap_locs)) trap_locs <- as.matrix(as.data.frame(trap_locs))
  
    ## general cost surface items 
    cost_surface <- exp(alpha_cost * landscape)
    transistion_surface <- geoCorrection(
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
        lik.marg[i] <- sum(lik.cond * (1/nG))
    }
    nv <- c(rep(1, length(lik.marg) - 1), n0)
    part1 <- lgamma(nrow(y) + n0 + 1) - lgamma(n0 + 1)
    part2 <- sum(nv * log(lik.marg))
    out <- -1 * (part1 + part2)
    out

    ## RSF likelihood
    if(add_telem){
      rsf_dmat <- costDistance(trLayer, ss_locs, ss_locs)
      rsf_pmat <- exp(1-telem_rsf_dmat^2/(2*sigma^2))
      rsf_pmat <- rsf_pmat/apply(rsf_pmat,1,sum)
    }
    
}

