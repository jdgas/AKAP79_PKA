# Uncertainty Quantification: Model simulation

runModel <- function(tpar, parIdx, parDef, input, rInd, amounts, datatype){

  # Set Parameters
  parDefVal = parDef;
  
  #Assign set of parameters to simulate
  pars <- parDefVal
  pars[parIdx] <- 10^tpar 
  
  # If ABC-MCMC is fitting Time-Series Data
  if (datatype == "T_S"){

    AKAR4p = c();
    
    #funtion called to simulate 1 of 3 conditons
    ode <- function(param, rInd, amounts){
      
      if (rInd=="SS_No_CaN") {
        #Set CaN cocnentration to zero
        amounts[12] = 0;
      }else if(rInd=="SS_with_CaN") {
        #Set CaN
        amounts[12] = 1.5;
      }else if (rInd=="SS_CaN_AKAP") {
        amounts[12] = 1.5;
        # Set AKAP on
        amounts[20] = 1;
      }
      
      # Set cAMP
      amounts[2] = param[31];

      # Simulate Model
      sol_BegEnd = lsode(y = amounts, times = c(0:605), func = PKAModel, parms = param) 
      df_sol_BegEnd = as.data.frame(sol_BegEnd);
      
      # Add colum to indicate ratio of C per second
      df_sol_BegEnd = cbind(df_sol_BegEnd, df_sol_BegEnd[,9]/df_sol_BegEnd[,16]);
      
      # Add colum to indicate ratio of AKAR4p per second
      df_sol_BegEnd = cbind(df_sol_BegEnd, df_sol_BegEnd[,24]/( df_sol_BegEnd[,22] + df_sol_BegEnd[,23] + df_sol_BegEnd[,24]) );
      
      # Add normalized ratio of AKAR4p
      df_sol_BegEnd = cbind(df_sol_BegEnd, (df_sol_BegEnd[,26] - min(df_sol_BegEnd[,26])) / (max(df_sol_BegEnd[,26])- min(df_sol_BegEnd[,26])) ); 
      colnames(df_sol_BegEnd) = c("time", "Rii", "cAMP", "RiiP", "Rii_C","RiiP_cAMP","RiiP_C","RiiP_C_cAMP","C",
                                  "Rii_cAMP","Rii_C_cAMP","Total_RII","CaN","RiiP_CaN","RiiP_cAMP_CaN","Total_C","cycle_Rii",
                                  "Cycle_RiiP","Thr_unphos_Rii","Thr_phos_Rii", "b_AKAP","AKAR4", "AKAR4_C", "AKAR4p", "ratioC",  "ratioAKAR4p", "Normal_AKAR4p");

      # Select time points in 5 seconds interval
      AKAR4p = cbind(AKAR4p, as.numeric( df_sol_BegEnd$AKAR4p[c(seq(1, 606, by = 5))]) );
      
    }
    params <- t(sapply(input, function(ip) c(pars,ip))); 
    vf <- apply(params, 1, function(params) c(ode(params, rInd, amounts)));
    xx <- seq(0, 605, by = 5); 
    yy <- as.numeric( vf );  
  }
  return(list(xx=xx, yy=yy))
}
