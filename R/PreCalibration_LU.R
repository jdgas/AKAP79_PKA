# Uncertainty Quantification: Precalibration for ABC-MCMC
# Copyright (C) 2018 Alexandra Jauhiainen (alexandra.jauhiainen@gmail.com)
# Modified slightly 2021 by Joao Antunes (joaodgantunes@gmail.com) and Olivia Eriksson (olivia@kth.se)


# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

preCalibration <- function(input, parIdx, parDef, xtarget, ytarget, ytarget_min, ytarget_max, yy_min, yy_max, npc, U, Z, copula, rInd, amounts, datatype, nChains){
 
  np <- length(parIdx) #nr of parameters
  R <- RVineSim(npc, copula) #gets different params values from hyper structure 'copula'
  prePar <- matrix(0, npc, np)
  for(i in 1:np){
    prePar[,i] = spline(Z[,i],U[,i],xout=R[,i])$y #'fusion' of different params, i.e., join information from normal distribution & CDF
  }
  
  fun <- function(i) {
    tpar <- prePar[i,]
    invisible(capture.output(out <- withTimeout(runModel(tpar, parIdx, parDef, input, rInd, amounts, datatype), timeout=25, onTimeout="silent"))) #simulates runModel for npc samples
    tmp <- ifelse(is.null(out), NA, getScore(ytarget, out$yy, ytarget_min, ytarget_max, yy_min, yy_max)) #Evaluates how "far" the simulated data is from the exp. dataset
    return(tmp)
  }
  
  #preDelta <- lapply(1:npc, fun) #single chain
  preDelta <- mclapply(1:npc, fun, mc.preschedule = FALSE, mc.cores = nChains) 
  preDelta <- unlist(preDelta)
  
  return(list(preDelta=preDelta, prePar=prePar))
}

getMCMCPar <- function(prePar, preDelta, p, sfactor, delta, nChains){
  prePar <- prePar[!is.na(preDelta),]
  preDelta <- preDelta[!is.na(preDelta)]
  nk <- nrow(prePar)*p
  pick1  <- grep(TRUE, (preDelta <= delta)) # pick set of params that are under threshold 'delta'
  pick2 <- order(preDelta, decreasing = F)[1:nk] # if cannot find above, pick top p percent 
  
  if(length(pick1)>length(pick2)){
    pick <- pick1
  }else{
    pick <- pick2
  }
  
  Scorr <- cor(prePar[pick,])*0.8 #tone down corrs
  diag(Scorr) <- 1
  sdv <- apply(prePar[pick,], 2, sd)
  Sigma <- sfactor * Scorr * tcrossprod(sdv)
  #browser()
  startPar <- prePar[sample(pick, nChains, replace = FALSE),] #pics startting point
  list(Sigma=Sigma, startPar=startPar)
}
