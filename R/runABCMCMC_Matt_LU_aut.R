# Uncertainty Quantification: Running ABC-MCMC with copulas
# Copyright (C) 2018 Alexandra Jauhiainen (alexandra.jauhiainen@gmail.com)
# Modified 2021 by Joao Antunes (joaodgantunes@gmail.com) and Olivia Eriksson (olivia@kth.se)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.


## Delete/Comment Working Directory specification
#setwd("/Users/olivia/GitHub/AKAP79_PKA")
#setwd("C:/Users/Joao/Desktop/Paper SWE/Only Joao & Olivia/ABC - All cAMP levels")


## Load R Scripts and Libraries
source('ScoringFunction_LU.R')
source('LoadTargets_AKAP79.R')
source('runModel_AKAP79.R')
source('copulaFunctions_LU.R')
source('ABCMCMCFunctions_LU.R')
source('PreCalibration_LU.R')
source('PKAModel.R')
library(deSolve)
library(ks)
library(VineCopula)
library(MASS)
library(R.utils)
library(R.matlab) 
library(stringr)


## Use if Running on a Cluster. Changes to be made in the following
## scripts: copulaFunctions.R, preCalibration.R and in this script (runABCMCMC_Matt_LU_aut.R)
library(doMC)
registerDoMC(19) # example with 19 cores throughout
getDoParWorkers()


# Load Parameters and Rii Amounts
matdata_Parms <<- read.table("pkaParms_Restri.txt", header = TRUE);
parNames = matdata_Parms$Name.;
parVal = matdata_Parms$Value.;


# Indexes of Experiments to Run
exp_idxs<-c(7, 13, 19, 22, 8, 14, 20, 23, 9, 15, 21, 24)
# 7: 0 No_CaN
# 8: 0 With_CaN
# 9: 0 CaN_AKAP
# 10: 0.1 No_CaN
# 11: 0.1 With_CaN
# 12: 0.1 CaN_AKAP
# 13: 0.2 No_CaN *
# 14: 0.2 With_CaN *
# 15: 0.2 CaN_AKAP *
# 16: 0.5  No_CaN
# 17: 0.5 With_CaN
# 18: 0.5 CaN_AKAP
# 19:1  No_CaN *
# 20:1 With_CaN *
# 21:1 CaN_AKAP *
# 22:2  No_CaN *
# 23:2 With_CaN *
# 24:2 CaN_AKAP *


# Input (cAMP) Values for Each Experimental Setting 
xAll <- c(-1,-1, -1,-1,-1,-1, 0, 0, 0, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.5, 0.5, 0.5, 1, 1, 1, 2, 2, 2)


# Exp type
exp_types<-c("","","","","","","SS_No_CaN", "SS_with_CaN", "SS_CaN_AKAP","SS_No_CaN", "SS_with_CaN", "SS_CaN_AKAP","SS_No_CaN", "SS_with_CaN", "SS_CaN_AKAP",
                               "SS_No_CaN", "SS_with_CaN", "SS_CaN_AKAP","SS_No_CaN", "SS_with_CaN", "SS_CaN_AKAP","SS_No_CaN", "SS_with_CaN", "SS_CaN_AKAP")




# Define Initial Concentration/Amount of Each Specie
R_Amount = 6.93 # Use to calculate the species' amount 
species_iconc_TS <- c(Rii = 6.3, cAMP = 0, RiiP = 0, Rii_C = 0.63,RiiP_cAMP = 0,RiiP_C = 0,RiiP_C_cAMP = 0,C = 0,
                      Rii_cAMP = 0,Rii_C_cAMP = 0, Total_RII = R_Amount, CaN = 1.5, RiiP_CaN = 0, RiiP_cAMP_CaN  = 0,Total_C = 0.63, cycle_Rii = 0,
                      Cycle_RiiP = 0,Thr_unphos_Rii  = 1,Thr_phos_Rii = 1, b_AKAP = 0, AKAR4 = 0.2, AKAR4_C = 0, AKAR4p = 0);



# Set ll and up for standard parameters
scale <- 1000


#Define Lower and Upper Limits for which the Parameters will be Allowed to Vary 
ll = c( parVal[1:6]/scale, parVal[7]/scale, parVal[8:9]/scale, parVal[10]/scale, parVal[11:20]/scale, parVal[21]/1.9, parVal[22:24]/scale,  parVal[25:26]/1.25, parVal[27]/1.25,  parVal[28:29]/1.5,  parVal[30]/2);
ul = c( parVal[1:6]*scale, parVal[7]*scale, parVal[8:9]*scale, parVal[10]*scale, parVal[11:20]*scale, parVal[21]*1.9, parVal[22:24]*scale,  parVal[25:26]*1.25, parVal[27]*1.25,  parVal[28:29]*1.5,  parVal[30]*2);


#Transform ll ul to a Logarithmic Scale
ll = log10(ll);
ul = log10(ul);


#Indexes of Parameters Used
parIdx <- c(1:30); 


# Load Experimental Values
out <- loadTargets_TSCC()
xtarget <- out$xtarget
ytarget <- out$ytarget
rm(out)


ytarget_min<-0 #Zero phosphorylation
ytarget_max<-0.2 #Full phosphorylation
yy_min<-0 #Zero phosphorylation
yy_max<-0.2  #Full phosphorylation


# Define Number of Samples for the Precalibration (npc) and each ABC-MCMC chain (ns)
ns <- 1000 # no of samples required from each ABC-MCMC chain #WAS 1000
npc <- 50000 # pre-calibration  WAS 50.000


# Define ABC-MCMC Settings
p <- 0.01 # For the Pre-Calibration: Choose Top 1% Samples with Shortest Distance to the Experimental Values
nChains <- 19 # For the ABC-MCMC Process: Nr of Parallel Chains; 
delta <- c(-1,-1, -1,-1,-1,-1, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)


# For Reproducibility of the ABC-MCMC Findings
set.seed(7619201)


# Loop through the Different Experimental Settings
for (i in 1:length(exp_idxs)) {
  
  exp_idx=exp_idxs[i]
  cat(sprintf("#####Starting run for Exp. Dataset %i #####\n", exp_idx))
  
  
  ## If First Experimental Setting, Create an Independente Colupla
  if(i==1){
    cat(sprintf("-Fitting independent Copula \n"))
    out <- makeIndepCopula_JA(ll, ul)
  ## Otherwise, Take Copula from the Previous Exp Setting and Use as a Prior
  } else {
    cat(sprintf("-Fitting Copula based on previous MCMC runs\n"))
    out <- fitCopula(draws, ll, ul, nChains)
  }
  copula <- out$copula
  U <- out$U
  Z <- out$Z
  Y <- out$Y
  
  ## Run Pre-Calibration Sampling
  cat(sprintf("-Precalibration \n"))
  if (exp_idx >=7 & exp_idx <= 24) {
    out1 <- preCalibration(xAll[[exp_idx]], parIdx, parVal, xtarget[[exp_idx]], ytarget[[exp_idx]],ytarget_min, ytarget_max, yy_min, yy_max, npc, U, Z, copula, exp_types[exp_idx], species_iconc_TS, "T_S", nChains);
  }
  else {
    cat(sprintf("OBS:Index %i is not an experiment (yet)\n", exp_idx))
  }
  sfactor <- 0.1 # scaling factor 

  ## Get Starting Parameters from Pre-Calibration
  out2 <- getMCMCPar(out1$prePar, out1$preDelta, p, sfactor, delta[exp_idx], nChains)
  Sigma <- out2$Sigma
  startPar <- out2$startPar
  
  ## Run ABC-MCMC Sampling
  cat(sprintf("-Running MCMC chains \n"))
  # If Running on a Cluster 
  draws <- mclapply(1:nChains, function(k) ABCMCMC(xAll[exp_idx], parIdx, parVal, ns, xtarget[[exp_idx]], ytarget[[exp_idx]], ytarget_min, ytarget_max, yy_min, yy_max, startPar[k,], Sigma, delta[exp_idx], U, Z, Y, copula, exp_types[exp_idx], ll, ul, species_iconc_TS, "T_S"), mc.preschedule = FALSE, mc.cores = nChains );
  # Otherwise 
  #draws <- lapply(1:nChains, function(k) ABCMCMC(xAll[exp_idx], parIdx, parVal, ns, xtarget[[exp_idx]], ytarget[[exp_idx]], ytarget_min, ytarget_max, yy_min, yy_max, startPar[k,], Sigma, delta[exp_idx], U, Z, Y, copula, exp_types[exp_idx], ll, ul, species_iconc_TS, "T_S"));
  
  # put together
  draws <- do.call("rbind", draws)
  pick <- !apply(draws, 1, function(rw) all(rw==0))
  draws <- draws[pick,]
 
  for(j in 1:i){
    filt_idx=exp_idxs[j]
    cat(sprintf("-Checking fit with dataset %i \n", exp_idxs[j]))
    pick <- apply(draws,1, function(u){out <- runModel(u, parIdx = parIdx, parDef = parVal, input = xAll[filt_idx], 
                                                     rInd = exp_types[filt_idx], species_iconc_TS, "T_S");
                 getScore(ytarget[[filt_idx]], out$yy, ytarget_min, ytarget_max, yy_min, yy_max)}) <= delta[filt_idx]
    n1=nrow(draws)
    draws <- draws[pick,];
    n2=nrow(draws)
    nonFits=n1-n2;
    cat(sprintf("-- %i samples of posterior after dataset %i did notfit dataset %i \n",nonFits,exp_idx, filt_idx))
  }
                 
  # Save Resulting Samples to MATLAB and R files.
  cat(sprintf("-Saving sample \n"))
  outFile=paste(exp_idxs[1:i], collapse="-")
  outFileR=paste0("DrawsNoThermoScale1000_",outFile,".RData",collapse="_")
  outFileM=paste0("DrawsNoThermoScale1000_",outFile,".mat",collapse="_")
  save(draws, parNames, file=outFileR)
  writeMat(outFileM, samples=10^draws)
  
}



