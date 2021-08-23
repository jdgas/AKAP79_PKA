# Uncertainty Quantification: ABC-MCMC with copulas
# Copyright (C) 2018 Alexandra Jauhiainen (alexandra.jauhiainen@gmail.com)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

ABCMCMC <- function(input, parIdx, parDef, nSims, xtarget,ytarget, startPar, Sigma0, delta, U, Z, Y, copula, rInd, ll, ul, amounts, datatype){
  
  cat("Started chain.\n")
  Sigma1 <- 0.25*diag(diag(Sigma0))
  curDelta <- Inf
  np <- length(startPar)
  scount <- 1  

  curPar  <- startPar
  #browser()
  invisible(capture.output(out <- runModel(curPar, parIdx, parDef, input, rInd, amounts, datatype))) #generate simulated data with model params
  #browser()
  curDelta <- getMaxScore(xtarget, ytarget, out$xx, out$yy) #get current delta distance
  curPrior <- dprior(curPar, U, Z, Y, copula, ll, ul)
  
  draws <- matrix(0, nSims,np)
  
  n <- 0
  
  while (n < nSims){
    #browser()
    if(runif(1)<=0.95){ #generates random deviates - paper supp P2.3
      canPar <- mvrnorm(n=1, curPar, Sigma0)
    }else{
      canPar <- mvrnorm(n=1, curPar, Sigma1)
    }
    out <- parUpdate(input, parDef, curPar, canPar, curDelta, curPrior, xtarget, ytarget, delta, U, Z, Y, copula, rInd, ll, ul, amounts, datatype) #propose delta transition
    curPar <- out$curPar
    curDelta <- out$curDelta
    curPrior <- out$curPrior
    
    scount <- ifelse(!all(curPar == canPar), scount+1, 1)
    #browser()
    if (curDelta <= delta & all(curPar == canPar)){ #p2.5 & P2.8
      n=n+1
      draws[n,] = curPar
    }
    
    if(scount>500){ #terminate chain if stuck
      cat('Aborted chain.\n')
      return(draws)
    }		
  }
  cat("Finished chain.\n")
  return(draws)
}



dprior <- function(inx, U, Z, Y, copula, ll, ul){
  np <- length(inx)
  
  ed <- sapply(1:np, function(i) approx(U[,i], Z[,i], xout=inx[i])$y)
  mpdf <- sapply(1:np, function(i) approx(U[,i], Y[,i], xout=inx[i])$y)
  
  if(any(is.na(ed))|any(is.na(mpdf))){ # outside of copula defined limits
    jpdf <- 0
  }else if(!(all(inx >=ll) & all(inx<=ul))){ # outside of prior
    jpdf <- 0
  }else{
    jpdf <- RVinePDF(ed, copula, verbose = TRUE)*prod(mpdf)
  }
  
  return(jpdf)
}

parUpdate <- function(input, parDef, curPar, canPar, curDelta, curPrior, xtarg, ytarg, delta, U, Z, Y, copula, rInd, ll, ul, amounts, datatype){
  #browser()
  invisible(capture.output(out <- withTimeout(runModel(canPar, parIdx, parDef, input, rInd, amounts, datatype), timeout=25, onTimeout="silent"))) #simulate model with candidate params/delta
  if(is.null(out)){
    canDelta <- Inf
    canPrior <- 0
  }else{  	
    canDelta <- getMaxScore(xtarg, ytarg, out$xx, out$yy)
    canPrior <- dprior(canPar, U, Z, Y, copula, ll, ul)
  }
  
  if (canDelta <= max(delta,curDelta)){
    if (canPrior==0){
      h <- 0
    }else{
      h <- min(1,canPrior/curPrior) #P2.6
    }
    
    if (runif(1) <= h){ #P2.7
      curDelta <- canDelta
      curPrior <- canPrior
      curPar <- canPar
    }
  }
  list(curPar=curPar, curDelta=curDelta, curPrior=curPrior)
}

