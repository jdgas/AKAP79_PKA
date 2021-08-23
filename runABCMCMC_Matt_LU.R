# Uncertainty Quantification: Running ABC-MCMC with copulas
# Copyright (C) 2018 Alexandra Jauhiainen (alexandra.jauhiainen@gmail.com)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.


setwd("C:/Users/Joao/Desktop/Paper SWE/Only Joao & Olivia/ABC - 1 and 2 cAMP/delta 0.01")

source('ScoringFunction_LU.R')
source('LoadTargets_Matt.R')
source('runModel_LU.R')
source('copulaFunctions_LU.R')
source('ABCMCMCFunctions_LU.R')
source('PreCalibration_LU.R')
source("PKA_AKAP_Restri.R")
library(deSolve)
library(ks)
library(VineCopula)
library(MASS)
library(R.utils)
library(R.matlab) #closeAllConnections()
library(stringr)

# # use if on a cluster, changes to be made in
# # copulaFunctions.R, preCalibration.R and in this file
# library(doMC)
# registerDoMC(20) # example with 20 cores throughout
# getDoParWorkers()

# Load Parameters and Rii Amounts
matdata_Parms <<- read.table("pkaParms_Restri.txt", header = TRUE);
parNames = matdata_Parms$Name.;
parVal = matdata_Parms$Value.;
R_Amount = 6.93
species_iconc_TS <- c(Rii = 6.3, cAMP = 0, RiiP = 0, Rii_C = 0.63,RiiP_cAMP = 0,RiiP_C = 0,RiiP_C_cAMP = 0,C = 0,
                      Rii_cAMP = 0,Rii_C_cAMP = 0, Total_RII = R_Amount, CaN = 1.5, RiiP_CaN = 0, RiiP_cAMP_CaN  = 0, Total_C = 0.63, cycle_Rii = 0,
                      Cycle_RiiP = 0,Thr_unphos_Rii  = 1,Thr_phos_Rii = 1, b_AKAP = 0, AKAR4 = 0.2, AKAR4_C = 0, AKAR4p = 0, p_AKAP = 0);

species_iconc_CC <- c(Rii =0, cAMP = 0, RiiP = 0, Rii_C = 0, RiiP_cAMP = 0,RiiP_C = 0,RiiP_C_cAMP = 0,C = 0,
                      Rii_cAMP = 0,Rii_C_cAMP = 0, Total_RII = R_Amount, CaN = 0, RiiP_CaN = 0, RiiP_cAMP_CaN  = 0, Total_C = 0.63, cycle_Rii = 0,
                      Cycle_RiiP = 0,Thr_unphos_Rii  = 0,Thr_phos_Rii = 0, b_AKAP = 0, AKAR4 = 0.2, AKAR4_C = 0, AKAR4p = 0, p_AKAP = 0);

species_iconc_MutA <- c(Rii = 6.3, cAMP = 1, RiiP = 0, Rii_C = 0.63,RiiP_cAMP = 0,RiiP_C = 0,RiiP_C_cAMP = 0,C = 0,
                        Rii_cAMP = 0,Rii_C_cAMP = 0, Total_RII = R_Amount, CaN = 0, RiiP_CaN = 0, RiiP_cAMP_CaN  = 0, Total_C = 0.63, cycle_Rii = 0,
                        Cycle_RiiP = 0,Thr_unphos_Rii  = 1,Thr_phos_Rii = 1, b_AKAP = 0, AKAR4 = 0.2, AKAR4_C = 0, AKAR4p = 0, p_AKAP = 0);

species_iconc_MutE <- c(Rii = 0, cAMP = 1, RiiP = 6.3, Rii_C = 0,RiiP_cAMP = 0,RiiP_C = 0.63,RiiP_C_cAMP = 0,C = 0,
                        Rii_cAMP = 0,Rii_C_cAMP = 0, Total_RII = R_Amount, CaN = 0, RiiP_CaN = 0, RiiP_cAMP_CaN  = 0, Total_C = 0.63, cycle_Rii = 0,
                        Cycle_RiiP = 0,Thr_unphos_Rii  = 1,Thr_phos_Rii = 1, b_AKAP = 0, AKAR4 = 0.2, AKAR4_C = 0, AKAR4p = 0, p_AKAP = 0);


# Set ll and up for standard parameters
scale <- 100

ll = c( parVal[1:6]/scale, parVal[7]/2.5, parVal[8:9]/scale, parVal[10]/2.5, parVal[11:13]/scale, parVal[14]/1000 , parVal[15:20]/scale, parVal[21]/1.9, parVal[22:24]/scale,  parVal[25:27]/1.5, parVal[28:29]/10,  parVal[30]/2);
ul = c( parVal[1:6]*scale, parVal[7]*2.5, parVal[8:9]*scale, parVal[10]*2.5, parVal[11:20]*scale, parVal[21]*1.9, parVal[22:23]*scale, parVal[24]*1000,   parVal[25:27]*1.5, parVal[28:29]*10,  parVal[30]*2);
ll = log10(ll);
ul = log10(ul);
parIdx <- c(1:30); 


## input values - as up now, only simulate specific concentrations/timings
xAll <- c(0, 0.1, 0.2, 0.5, 1, 2);

# load targets
out <- loadTargets_TSCC()
xtarget <- out$xtarget
ytarget <- out$ytarget
rm(out)

# no of samples
ns <- 1000 # no of samples required from each ABC-MCMC chain #WAS 1000
npc <- 40000 # pre-calibration  WAS 50.000

# settings
p <- 0.1 # was 0.01
nChains <- 1 #was 20
delta <- 0.01 #was 0.1



# 19 (cAMP = 1) ####
##############################
## Run MCMC for Dataset   19##
##############################
out <- makeIndepCopula_JA(ll, ul)
copula <- out$copula
U <- out$U
Z <- out$Z
Y <- out$Y


# precalibration all at once
out1 <- preCalibration(1, parIdx, parVal, xtarget[[19]], ytarget[[19]], npc, U, Z, copula, 'SS_CaN_AKAP', species_iconc_TS, 'T_S', nChains);
sfactor <- 0.05 #scaling factor


# Get Starting Parameters
out2 <- getMCMCPar(out1$prePar, out1$preDelta, p, sfactor, delta, nChains)
Sigma <- out2$Sigma
startPar <- out2$startPar

# Run Chains
draws19 <- mclapply(1:nChains, function(i) ABCMCMC(1, parIdx, parVal, ns, xtarget[[19]], ytarget[[19]], startPar[i,], Sigma, delta, U, Z, Y, copula, 'SS_CaN_AKAP', ll, ul, species_iconc_TS,'T_S'), mc.preschedule = FALSE, mc.cores = nChains);

# put together
draws19 <- do.call("rbind", draws19)
pick <- !apply(draws19, 1, function(rw) all(rw==0))
draws19 <- draws19[pick,]

draws <- draws19
outFile <- sprintf('Draws-Dataset19-ABCMCMC-Scale%d.RData', scale)
save(draws, parNames, file=outFile)



# 20 (cAMP = 1) ####
##############################
## Run MCMC for Dataset   20##
##############################

out <- fitCopula(draws19, ll, ul, nChains)
copula <- out$copula
U <- out$U
Z <- out$Z
Y <- out$Y

# precalibration all at once
out1 <- preCalibration(1, parIdx, parVal, xtarget[[20]], ytarget[[20]], npc, U, Z, copula, 'SS_CaN_AKAP', species_iconc_TS, 'T_S', nChains);
sfactor <- 0.05 #scaling factor


# Get Starting Parameters
out2 <- getMCMCPar(out1$prePar, out1$preDelta, p, sfactor, delta, nChains)
Sigma <- out2$Sigma
startPar <- out2$startPar

# Run Chains
draws20 <- mclapply(1:nChains, function(i) ABCMCMC(1, parIdx, parVal, ns, xtarget[[20]], ytarget[[20]], startPar[i,], Sigma, delta, U, Z, Y, copula, 'SS_CaN_AKAP', ll, ul, species_iconc_TS,'T_S'), mc.preschedule = FALSE, mc.cores = nChains);

# put together
draws20 <- do.call("rbind", draws20)
pick <- !apply(draws20, 1, function(rw) all(rw==0))
draws20 <- draws20[pick,]

# filter
pick19 <- apply(draws20,1, function(u){out <- runModel(u, parIdx = parIdx, parDef = parVal, input = 1, rInd = 'SS_No_CaN', species_iconc_TS, 'T_S');getMaxScore(xtarget[[19]], ytarget[[19]], out$xx, out$yy)}) <= delta
draws20 <- draws20[pick19,];


draws <- draws20
outFile <- sprintf('Draws-Dataset20-ABCMCMC-Scale%d.RData', scale)
save(draws, parNames, file=outFile)



# 21 (cAMP = 1) ####
##############################
## Run MCMC for Dataset   21##
##############################

out <- fitCopula(draws20, ll, ul, nChains)
copula <- out$copula
U <- out$U
Z <- out$Z
Y <- out$Y

# precalibration all at once
out1 <- preCalibration(1, parIdx, parVal, xtarget[[21]], ytarget[[21]], npc, U, Z, copula, 'SS_CaN_AKAP', species_iconc_TS, 'T_S', nChains);
sfactor <- 0.05 #scaling factor


# Get Starting Parameters
out2 <- getMCMCPar(out1$prePar, out1$preDelta, p, sfactor, delta, nChains)
Sigma <- out2$Sigma
startPar <- out2$startPar

# Run Chains
draws21 <- mclapply(1:nChains, function(i) ABCMCMC(1, parIdx, parVal, ns, xtarget[[21]], ytarget[[21]], startPar[i,], Sigma, delta, U, Z, Y, copula, 'SS_CaN_AKAP', ll, ul, species_iconc_TS,'T_S'), mc.preschedule = FALSE, mc.cores = nChains);

# put together
draws21 <- do.call("rbind", draws21)
pick <- !apply(draws21, 1, function(rw) all(rw==0))
draws21 <- draws21[pick,]

# filter
pick19 <- apply(draws21,1, function(u){out <- runModel(u, parIdx = parIdx, parDef = parVal, input = 1, rInd = 'SS_No_CaN', species_iconc_TS, 'T_S');getMaxScore(xtarget[[19]], ytarget[[19]], out$xx, out$yy)}) <= delta
pick20 <- apply(draws21,1, function(u){out <- runModel(u, parIdx = parIdx, parDef = parVal, input = 1, rInd = 'SS_with_CaN', species_iconc_TS, 'T_S');getMaxScore(xtarget[[20]], ytarget[[20]], out$xx, out$yy)}) <= delta
draws21 <- draws21[pick19 & pick20,];


draws <- draws21
outFile <- sprintf('Draws-Dataset21-ABCMCMC-Scale%d.RData', scale)
save(draws, parNames, file=outFile)




# 22 (cAMP = 2) ####
##############################
## Run MCMC for Dataset   22##
##############################

out <- fitCopula(draws21, ll, ul, nChains)
copula <- out$copula
U <- out$U
Z <- out$Z
Y <- out$Y

# precalibration all at once
out1 <- preCalibration(2, parIdx, parVal, xtarget[[22]], ytarget[[22]], npc, U, Z, copula, 'SS_CaN_AKAP', species_iconc_TS, 'T_S', nChains);
sfactor <- 0.05 #scaling factor


# Get Starting Parameters
out2 <- getMCMCPar(out1$prePar, out1$preDelta, p, sfactor, delta, nChains)
Sigma <- out2$Sigma
startPar <- out2$startPar

# Run Chains
draws22 <- mclapply(1:nChains, function(i) ABCMCMC(2, parIdx, parVal, ns, xtarget[[22]], ytarget[[22]], startPar[i,], Sigma, delta, U, Z, Y, copula, 'SS_CaN_AKAP', ll, ul, species_iconc_TS,'T_S'), mc.preschedule = FALSE, mc.cores = nChains);

# put together
draws22 <- do.call("rbind", draws22)
pick <- !apply(draws22, 1, function(rw) all(rw==0))
draws22 <- draws22[pick,]

# filter
pick19 <- apply(draws22,1, function(u){out <- runModel(u, parIdx = parIdx, parDef = parVal, input = 1, rInd = 'SS_No_CaN', species_iconc_TS, 'T_S');getMaxScore(xtarget[[19]], ytarget[[19]], out$xx, out$yy)}) <= delta
pick20 <- apply(draws22,1, function(u){out <- runModel(u, parIdx = parIdx, parDef = parVal, input = 1, rInd = 'SS_with_CaN', species_iconc_TS, 'T_S');getMaxScore(xtarget[[20]], ytarget[[20]], out$xx, out$yy)}) <= delta
pick21 <- apply(draws22,1, function(u){out <- runModel(u, parIdx = parIdx, parDef = parVal, input = 1, rInd = 'SS_CaN_AKAP', species_iconc_TS, 'T_S');getMaxScore(xtarget[[21]], ytarget[[21]], out$xx, out$yy)}) <= delta
draws22 <- draws22[pick19 & pick20 & pick21,];


draws <- draws22
outFile <- sprintf('Draws-Dataset22-ABCMCMC-Scale%d.RData', scale)
save(draws, parNames, file=outFile)



# 23 (cAMP = 2) ####
##############################
## Run MCMC for Dataset   23##
##############################

out <- fitCopula(draws22, ll, ul, nChains)
copula <- out$copula
U <- out$U
Z <- out$Z
Y <- out$Y

# precalibration all at once
out1 <- preCalibration(2, parIdx, parVal, xtarget[[23]], ytarget[[23]], npc, U, Z, copula, 'SS_CaN_AKAP', species_iconc_TS, 'T_S', nChains);
sfactor <- 0.05 #scaling factor


# Get Starting Parameters
out2 <- getMCMCPar(out1$prePar, out1$preDelta, p, sfactor, delta, nChains)
Sigma <- out2$Sigma
startPar <- out2$startPar

# Run Chains
draws23 <- mclapply(1:nChains, function(i) ABCMCMC(2, parIdx, parVal, ns, xtarget[[23]], ytarget[[23]], startPar[i,], Sigma, delta, U, Z, Y, copula, 'SS_CaN_AKAP', ll, ul, species_iconc_TS,'T_S'), mc.preschedule = FALSE, mc.cores = nChains);

# put together
draws23 <- do.call("rbind", draws23)
pick <- !apply(draws23, 1, function(rw) all(rw==0))
draws23 <- draws23[pick,]

# filter
pick19 <- apply(draws23,1, function(u){out <- runModel(u, parIdx = parIdx, parDef = parVal, input = 1, rInd = 'SS_No_CaN', species_iconc_TS, 'T_S');getMaxScore(xtarget[[19]], ytarget[[19]], out$xx, out$yy)}) <= delta
pick20 <- apply(draws23,1, function(u){out <- runModel(u, parIdx = parIdx, parDef = parVal, input = 1, rInd = 'SS_with_CaN', species_iconc_TS, 'T_S');getMaxScore(xtarget[[20]], ytarget[[20]], out$xx, out$yy)}) <= delta
pick21 <- apply(draws23,1, function(u){out <- runModel(u, parIdx = parIdx, parDef = parVal, input = 1, rInd = 'SS_CaN_AKAP', species_iconc_TS, 'T_S');getMaxScore(xtarget[[21]], ytarget[[21]], out$xx, out$yy)}) <= delta
pick22 <- apply(draws23,1, function(u){out <- runModel(u, parIdx = parIdx, parDef = parVal, input = 2, rInd = 'SS_No_CaN', species_iconc_TS, 'T_S');getMaxScore(xtarget[[22]], ytarget[[22]], out$xx, out$yy)}) <= delta
draws23 <- draws23[pick19 & pick20 & pick21 & pick22,];


draws <- draws23
outFile <- sprintf('Draws-Dataset23-ABCMCMC-Scale%d.RData', scale)
save(draws, parNames, file=outFile)



# 24 (cAMP = 2) ####
##############################
## Run MCMC for Dataset   24##
##############################

out <- fitCopula(draws23, ll, ul, nChains)
copula <- out$copula
U <- out$U
Z <- out$Z
Y <- out$Y

# precalibration all at once
out1 <- preCalibration(2, parIdx, parVal, xtarget[[24]], ytarget[[24]], npc, U, Z, copula, 'SS_CaN_AKAP', species_iconc_TS, 'T_S', nChains);
sfactor <- 0.05 #scaling factor


# Get Starting Parameters
out2 <- getMCMCPar(out1$prePar, out1$preDelta, p, sfactor, delta, nChains)
Sigma <- out2$Sigma
startPar <- out2$startPar

# Run Chains
draws24 <- mclapply(1:nChains, function(i) ABCMCMC(2, parIdx, parVal, ns, xtarget[[24]], ytarget[[24]], startPar[i,], Sigma, delta, U, Z, Y, copula, 'SS_CaN_AKAP', ll, ul, species_iconc_TS,'T_S'), mc.preschedule = FALSE, mc.cores = nChains);

# put together
draws24 <- do.call("rbind", draws24)
pick <- !apply(draws24, 1, function(rw) all(rw==0))
draws24 <- draws24[pick,]

# filter
pick19 <- apply(draws24,1, function(u){out <- runModel(u, parIdx = parIdx, parDef = parVal, input = 1, rInd = 'SS_No_CaN', species_iconc_TS, 'T_S');getMaxScore(xtarget[[19]], ytarget[[19]], out$xx, out$yy)}) <= delta
pick20 <- apply(draws24,1, function(u){out <- runModel(u, parIdx = parIdx, parDef = parVal, input = 1, rInd = 'SS_with_CaN', species_iconc_TS, 'T_S');getMaxScore(xtarget[[20]], ytarget[[20]], out$xx, out$yy)}) <= delta
pick21 <- apply(draws24,1, function(u){out <- runModel(u, parIdx = parIdx, parDef = parVal, input = 1, rInd = 'SS_CaN_AKAP', species_iconc_TS, 'T_S');getMaxScore(xtarget[[21]], ytarget[[21]], out$xx, out$yy)}) <= delta
pick22 <- apply(draws24,1, function(u){out <- runModel(u, parIdx = parIdx, parDef = parVal, input = 2, rInd = 'SS_No_CaN', species_iconc_TS, 'T_S');getMaxScore(xtarget[[22]], ytarget[[22]], out$xx, out$yy)}) <= delta
pick23 <- apply(draws24,1, function(u){out <- runModel(u, parIdx = parIdx, parDef = parVal, input = 2, rInd = 'SS_with_CaN', species_iconc_TS, 'T_S');getMaxScore(xtarget[[23]], ytarget[[23]], out$xx, out$yy)}) <= delta
draws24 <- draws24[pick19 & pick20 & pick21 & pick22 & pick23,];


draws <- draws24
outFile <- sprintf('Draws-Dataset24-ABCMCMC-Scale%d.RData', scale)
save(draws, parNames, file=outFile)

