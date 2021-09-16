# Uncertainty Quantification: Model simulation

runModel <- function(tpar, parIdx, parDef, input, rInd, amounts, datatype){

  # Set Parameters
  parDefVal = parDef;
  
  #Assign set of parameters to simulate
  pars <- parDefVal
  pars[parIdx] <- 10^tpar 
  
  if (datatype == "D_R"){
    
      ratioC = c();
      # funtion to be called to simulate 1 of 3 conditons
      ode <- function(param, rInd, amounts){
    
        if (rInd=="SS_No_CaN") {
          #Set CaN cocnentration to zero
          amounts[12] = 0;
        }else if(rInd=="SS_with_CaN") {
          #Set CaN
          amounts[12] = 1.5;
        }else if (rInd=="SS_CaN_AKAP") {
          #Set CaN
          amounts[12] = 1.5;
          # Set AKAP on
          amounts[20] = 1;
        }
        
        # Set cAMP
        amounts[2] = param[31];

        # Simulate Relaxation
        sol_BegEnd = lsode(y = amounts, times = c(0,2000), func = PKAModel, parms = param )[, -c(1)];
        df_sol_BegEnd = as.data.frame(sol_BegEnd);
        
        # Add colum to indicate ratio of C per second
        df_sol_BegEnd = cbind(df_sol_BegEnd, df_sol_BegEnd[,8]/df_sol_BegEnd[,15]); # (i.e., C divided by Total_C)
        
        colnames(df_sol_BegEnd) = c("Rii", "cAMP", "RiiP", "Rii_C","RiiP_cAMP","RiiP_C","RiiP_C_cAMP","C",
                                  "Rii_cAMP","Rii_C_cAMP","Total_RII","CaN","RiiP_CaN","RiiP_cAMP_CaN","Total_C","cycle_Rii",
                                  "Cycle_RiiP","Thr_unphos_Rii","Thr_phos_Rii", "b_AKAP","AKAR4", "AKAR4_C", "AKAR4p", "ratioC");

        # Ratio of C 
        ratioC = cbind(ratioC, df_sol_BegEnd$ratioC[2]);
      }
    
      params <- t(sapply(input, function(ip) c(pars,ip))); 
      vf <- apply(params, 1, function(params) c(ode(params, rInd, amounts)));
    
      xx <- input; 
      yy <- vf; 
    
  }else if (datatype == "T_S"){

    AKAR4p = c();
    #funtion to be called to simulate 1 of 3 conditons
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
      sol_BegEnd = lsode(y = amounts, times = c(0:605), func = PKAModel, parms = param) #, rtol = 1e300, atol = 1e00)
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

      AKAR4p = cbind(AKAR4p, as.numeric( df_sol_BegEnd$AKAR4p[c(seq(1, 606, by = 5))]) );
      
    }
    params <- t(sapply(input, function(ip) c(pars,ip))); 
    vf <- apply(params, 1, function(params) c(ode(params, rInd, amounts)));
    xx <- seq(0, 605, by = 5); 
    yy <- as.numeric( vf ); 
    
  }else if (datatype == "C_C"){
    
    AKAR4p = c();
    #funtion to be called to simulate 1 of 3 conditons
    ode <- function(param, rInd, amounts){
      
      # Set C
      amounts[8] = param[31];
      amounts[15] = param[31];
      
      # Simulate Model
      sol_BegEnd = lsode(y = amounts, times = c(0:1105), func = PKAModel, parms = param) #, rtol = 1e300, atol = 1e00)
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
      
      AKAR4p = cbind(AKAR4p, as.numeric( df_sol_BegEnd$AKAR4p[c(seq(1, 1106, by = 5))]) );
      
    }
    params <- t(sapply(input, function(ip) c(pars,ip))); 
    vf <- apply(params, 1, function(params) c(ode(params, rInd, amounts)));
    xx <- seq(0, 1105, by = 5); 
    yy <- as.numeric( vf ); 
    
  } 
  return(list(xx=xx, yy=yy))
}
