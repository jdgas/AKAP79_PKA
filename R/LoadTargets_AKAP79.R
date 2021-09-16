loadTargets_TSCC <- function(){
  
  
  #setwd("C:/Users/Joao/Desktop/Paper SWE/Only Joao & Olivia/ABC - All cAMP levels/Targets")
  #browser()
  
  xtarget <- list()
  ytarget <- list()
  
  
  #### Calibration Curves Unnormalized Data ####
  # Read txt file for CCurves
  A <- read.csv('Targets/Calib_Curves.txt', sep='\t')
  # Normalize data with constant 1.08
  #A[,c(2:7)] = A[,c(2:7)] / 1.08
  
  # Normalize data and re-assign it to 'C'
  First = ( (A$X0.4 - 100 ) / ( 171.67 - 100 ) ) * 0.2  ;
  Second = ( (A$X0.2 - 100 ) / ( 171.67 - 100 ) ) * 0.2  ;
  Third = ( (A$X0.1 - 100 ) / ( 171.67 - 100 ) ) * 0.2  ;
  Fourth = ( (A$X0.05 - 100 ) / ( 171.67 - 100 ) ) * 0.2  ;
  Fith = ( (A$X0.025 - 100 ) / ( 171.67 - 100 ) ) * 0.2  ;
  Sixth = ( (A$X0 - 100 ) / ( 171.67 - 100 ) ) * 0.2  ;
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
  B <- read.csv('Targets/TS_0cAMP.txt', sep='\t');
  # Normalize data with constant 1.08
  #B[,c(2:4)] = B[,c(2:4)] / 1.08
  
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
  C <- read.csv('Targets/TS_01cAMP.txt', sep='\t');
  # Normalize data with constant 1.08
  #C[,c(2:4)] = C[,c(2:4)] / 1.08
 
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
  D <- read.csv('Targets/TS_02cAMP.txt', sep='\t');
  # Normalize data with constant 1.08
  #D[,c(2:4)] = D[,c(2:4)] / 1.08
 
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
  E <- read.csv('Targets/TS_05cAMP.txt', sep='\t');
  # Normalize data with constant 1.08
  #E[,c(2:4)] = E[,c(2:4)] / 1.08
  
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
  G <- read.csv('Targets/TS_1cAMP.txt', sep='\t');
  # Normalize data with constant 1.08
  #G[,c(2:4)] = G[,c(2:4)] / 1.08
  
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
  H <- read.csv('Targets/TS_2cAMP.txt', sep='\t')
  # Normalize data with constant 1.08
  #H[,c(2:4)] = H[,c(2:4)] / 1.08
  
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
  
  
  
  
  # #### Read txt fileS for Mutation Data ####
  # I <- read.csv('Targets/Mut1.txt', sep='\t') # 98A
  # # Normalize data with constant 1.08
  # #I[,c(2:4)] = I[,c(2:4)] / 1.08
  # 
  # # Normalize data and re-assign it to 'C'
  # First = ((I$WTCaN - 100) / ( 171.67 - 100)) * 0.2  ;
  # Second = ((I$WTalone - 100 ) / ( 171.67 - 100 )) * 0.2  ;
  # Third = ((I$NEACaN - 100 ) / ( 171.67 - 100 )) * 0.2  ;
  # normAKAR = cbind(I$Time, First, Second, Third);
  # colnames(normAKAR) = c("Time", "WTCaN", "WTalone", "NEACaN");
  # I = as.data.frame(normAKAR);
  # 
  # xtarget[[25]] <- I[,1]
  # ytarget[[25]] <- I[,2]; 
  # xtarget[[26]] <- I[,1]
  # ytarget[[26]] <- I[,3]; 
  # xtarget[[27]] <- I[,1]
  # ytarget[[27]] <- I[,4]; 
  # 
  # 
  # 
  # 
  # J <- read.csv('Targets/Mut2.txt', sep='\t') # 98E
  # # Normalize data with constant 1.08
  # #J[,c(2:3)] = J[,c(2:3)] / 1.08
  # # Normalize data and re-assign it to 'C'
  # First = ((J$EplusCaN - 100) / ( 171.67 - 100)) * 0.2  ;
  # Second = ((J$Ealone - 100 ) / ( 171.67 - 100 )) * 0.2  ;
  # normAKAR = cbind(J$Time, First, Second);
  # colnames(normAKAR) = c("Time", "EplusCaN", "Ealone");
  # J = as.data.frame(normAKAR);
  # 
  # xtarget[[28]] <- J[,1]
  # ytarget[[28]] <- J[,2]; 
  # xtarget[[29]] <- J[,1]
  # ytarget[[29]] <- J[,3]; 
  # 
  
  return(list(xtarget=xtarget, ytarget=ytarget))
  
}
