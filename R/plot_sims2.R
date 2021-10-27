# plot_sims2
# Copyright (C) 2021 by Joao Antunes (joaodgantunes@gmail.com) 

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

#plot_draws=draws
data=readMat('DrawsNoThermoScale1000_7-13-19-22-8-14-20-23-9-15-21-24_large_joao_subm_prior.mat');
plot_draws=log10(data$samples)
#nsamples=dim(plot_draws)[1]
nsamples=1000

species_iconc_TS <- c(Rii = 6.3, cAMP = 0, RiiP = 0, Rii_C = 0.63,RiiP_cAMP = 0,RiiP_C = 0,
                      RiiP_C_cAMP = 0,C = 0,Rii_cAMP = 0,Rii_C_cAMP = 0, Total_RII = R_Amount, CaN = 0,
                      RiiP_CaN = 0, RiiP_cAMP_CaN  = 0, Total_C = 0.63, cycle_Rii = 0,
                      Cycle_RiiP = 0,Thr_unphos_Rii  = 1,Thr_phos_Rii = 1, b_AKAP = 0,
                      AKAR4 = 0.2, AKAR4_C = 0, AKAR4p = 0);#, p_AKAP = 0)


for(i in 1:12)
{ 
  p_idx=exp_idxs[i]
  ns=nrow(plot_draws) 
  sims <- apply(plot_draws[1:nsamples,],1, 
                function(u){out <- runModel(u, parIdx = parIdx, parDef = parVal, 
                                            input = xAll[p_idx], 
                                            rInd = exp_types[p_idx], 
                                            species_iconc_TS, 'T_S')})
  title=paste(exp_types[p_idx],"cAMP=",xAll[p_idx], collapse=" ")
  matplot(sims[[1]]$xx, sims[[1]]$yy,type="l",pch=19, col="red", xlim=c(0,605), 
          ylim=c(0,0.25), main=title, xlab="t (s)", ylab="AKAR4p")
  for (j in 1:nsamples){
    matplot(sims[[j]]$xx, sims[[j]]$yy, type="l",pch=19, col="red", add=TRUE)
  }
  matplot(xtarget[[p_idx]], ytarget[[p_idx]],pch=19,add=TRUE)
  #boxplot(plot_draws, main=title)
}



#