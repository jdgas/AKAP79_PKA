# Load targets
# Copyright (C) 2021 by Joao Antunes (joaodgantunes@gmail.com) 

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.


loadTargets_TSCC <- function(){
  #browser()
  
  xtarget <- list()
  ytarget <- list()
  
  
  #### Calibration Curves Unnormalized Data ####
  # Read txt file for CCurves
  A <- read.csv('../Targets/Calib_Curves.txt', sep='\t')

  
  # Normalize data and re-assign it to 'C'
  First = ( (A$C0.4 - 100 ) / ( 171.67 - 100 ) ) * 0.2  ;
  Second = ( (A$C0.2 - 100 ) / ( 171.67 - 100 ) ) * 0.2  ;
  Third = ( (A$C0.1 - 100 ) / ( 171.67 - 100 ) ) * 0.2  ;
  Fourth = ( (A$C0.05 - 100 ) / ( 171.67 - 100 ) ) * 0.2  ;
  Fith = ( (A$C0.025 - 100 ) / ( 171.67 - 100 ) ) * 0.2  ;
  Sixth = ( (A$C0 - 100 ) / ( 171.67 - 100 ) ) * 0.2  ;
  normAKAR = cbind(A$Time, First, Second, Third, Fourth, Fith, Sixth);
  colnames(normAKAR) = c("Time", "0.4C", "0.2C", "0.1C", "0.05C", "0.025C", "0C");
  A = as.data.frame(normAKAR);
  
  # 0.4 uM
  xtarget[[1]] <- A[,1]
  ytarget[[1]] <- A[,2]; 
  # 0.2 uM
  xtarget[[2]] <- A[,1]
  ytarget[[2]] <- A[,3]; 
  # 0.1 uM
  xtarget[[3]] <- A[,1]
  ytarget[[3]] <- A[,4];
  # 0.05 uM
  xtarget[[4]] <- A[,1]
  ytarget[[4]] <- A[,5]; 
  # 0.025 uM
  xtarget[[5]] <- A[,1]
  ytarget[[5]] <- A[,6]; 
  # 0 uM
  xtarget[[6]] <- A[,1]
  ytarget[[6]] <- A[,7];
  
  
  
  #### Time-Series Normalized Data ####
  # Read txt file for 0uM cAMP
  B <- read.csv('../Targets/0_cAMP.txt', sep='\t');

  
  # Normalize data between 0 and 1
  First = ((B$NoCaNnoAKAP - 100) / ( 171.67 - 100)) * 0.2 ;
  Second = ((B$CaNonly - 100) / ( 171.67 - 100)) * 0.2 ;
  Third = ((B$CaNAKAP - 100) / ( 171.67 - 100)) * 0.2  ;
  normAKAR = cbind(B$Time, First, Second, Third);
  colnames(normAKAR) = c("Time", "AKAPoffNoCaN", "AKAPoffCaN", "AKAPon");
  B = as.data.frame(normAKAR);
  
  # No CaN & No AKAP
  xtarget[[7]] <- B[,1]
  ytarget[[7]] <- B[,2]; 
  # CaN Only
  xtarget[[8]] <- B[,1]
  ytarget[[8]] <- B[,3]; 
  # With Can & AKAP
  xtarget[[9]] <- B[,1]
  ytarget[[9]] <- B[,4]; 
  
  
  
  # Read txt file for 0.1uM cAMP
  C <- read.csv('../Targets/0.1_cAMP.txt', sep='\t');

 
  First = ((C$NoCaNnoAKAP - 100) / ( 171.67 - 97)) * 0.2 ;
  Second = ((C$CaNonly - 100) / ( 171.67 - 97)) * 0.2 ;
  Third = ((C$CaNAKAP - 100) / ( 171.67 - 97)) * 0.2  ;
  normAKAR = cbind(C$Time, First, Second, Third);
  colnames(normAKAR) = c("Time", "AKAPoffNoCaN", "AKAPoffCaN", "AKAPon");
  C = as.data.frame(normAKAR);
  
  # No CaN & No AKAP
  xtarget[[10]] <- C[,1]
  ytarget[[10]] <- C[,2]; 
  # CaN Only
  xtarget[[11]] <- C[,1]
  ytarget[[11]] <- C[,3]; 
  # With Can & AKAP
  xtarget[[12]] <- C[,1]
  ytarget[[12]] <- C[,4]; 
  
  
 
  # Read txt file for 0.2uM cAMP
  D <- read.csv('../Targets/0.2_cAMP.txt', sep='\t');

 
  First = ((D$NoCaNnoAKAP - 100) / ( 171.67 - 100)) * 0.2 ;
  Second = ((D$CaNonly - 100) / ( 171.67 - 100)) * 0.2 ;
  Third = ((D$CaNAKAP - 100) / ( 171.67 - 100)) * 0.2  ;
  normAKAR = cbind(D$Time, First, Second, Third);
  colnames(normAKAR) = c("Time", "AKAPoffNoCaN", "AKAPoffCaN", "AKAPon");
  D = as.data.frame(normAKAR);
  
  # No CaN & No AKAP
  xtarget[[13]] <- D[,1]
  ytarget[[13]] <- D[,2]; 
  # CaN Only
  xtarget[[14]] <- D[,1]
  ytarget[[14]] <- D[,3]; 
  # With Can & AKAP
  xtarget[[15]] <- D[,1]
  ytarget[[15]] <- D[,4]; 
  
  
  
  # Read txt file for 0.5uM cAMP
  E <- read.csv('../Targets/0.5_cAMP.txt', sep='\t');

  
  First = ((E$NoCaNnoAKAP - 100) / ( 171.67 - 100)) * 0.2 ;
  Second = ((E$CaNonly - 100) / ( 171.67 - 100)) * 0.2 ;
  Third = ((E$CaNAKAP - 100) / ( 171.67 - 100)) * 0.2  ;
  normAKAR = cbind(E$Time, First, Second, Third);
  colnames(normAKAR) = c("Time", "AKAPoffNoCaN", "AKAPoffCaN", "AKAPon");
  E = as.data.frame(normAKAR);
  
  # No CaN & No AKAP
  xtarget[[16]] <- E[,1]
  ytarget[[16]] <- E[,2]; 
  # CaN Only
  xtarget[[17]] <- E[,1]
  ytarget[[17]] <- E[,3]; 
  # With Can & AKAP
  xtarget[[18]] <- E[,1]
  ytarget[[18]] <- E[,4]; 
  
  
  
  #### Time-Series Unnormalized Data ####
  # Read txt file for 1uM cAMP
  G <- read.csv('../Targets/1_cAMP.txt', sep='\t');

  
  First = ((G$NoCaNnoAKAP - 100) / ( 171.67 - 100)) * 0.2 ;
  Second = ((G$CaNonly - 100) / ( 171.67 - 100)) * 0.2 ;
  Third = ((G$CaNAKAP - 100) / ( 171.67 - 100)) * 0.2  ;
  normAKAR = cbind(G$Time, First, Second, Third);
  colnames(normAKAR) = c("Time", "AKAPoffNoCaN", "AKAPoffCaN", "AKAPon");
  G = as.data.frame(normAKAR);
  
  # No CaN & No AKAP
  xtarget[[19]] <- G[,1]
  ytarget[[19]] <- G[,2]; 
  # CaN Only
  xtarget[[20]] <- G[,1]
  ytarget[[20]] <- G[,3]; 
  # With Can & AKAP
  xtarget[[21]] <- G[,1]
  ytarget[[21]] <- G[,4]; 
  
  
  
  # Read txt file for 2uM cAMP
  H <- read.csv('../Targets/2_cAMP.txt', sep='\t')

  
  First = ((H$NoCaNnoAKAP - 100) / ( 171.67 - 100)) * 0.2  ;
  Second = ((H$CaNonly - 100 ) / ( 171.67 - 100 )) * 0.2  ;
  Third = ((H$CaNAKAP - 100 ) / ( 171.67 - 100 )) * 0.2  ;
  normAKAR = cbind(H$Time, First, Second, Third);
  colnames(normAKAR) = c("Time", "AKAPoffNoCaN", "AKAPoffCaN", "AKAPon");
  H = as.data.frame(normAKAR);
  
  # No CaN & No AKAP
  xtarget[[22]] <- H[,1]
  ytarget[[22]] <- H[,2]; 
  # CaN Only
  xtarget[[23]] <- H[,1]
  ytarget[[23]] <- H[,3]; 
  # With Can & AKAP
  xtarget[[24]] <- H[,1]
  ytarget[[24]] <- H[,4]; 
  
  
  return(list(xtarget=xtarget, ytarget=ytarget))
  
}
