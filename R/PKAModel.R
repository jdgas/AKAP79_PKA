# R implementation of AKAP79-PKA model
# Copyright (C) 2021 Joao Antunes (joaodgantunes@gmail.com) and Olivia Eriksson (olivia@kth.se)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
PKAModel <- function(t, state, parameters) {

  # Species
  Rii = state[1];
  cAMP = state[2];
  RiiP = state[3]
  Rii_C = state[4];
  RiiP_cAMP = state[5];
  RiiP_C = state[6];
  RiiP_C_cAMP = state[7];
  C = state[8];
  Rii_cAMP = state[9];
  Rii_C_cAMP = state[10];
  CaN = state[12];
  RiiP_CaN = state[13];
  RiiP_cAMP_CaN = state[14];
  b_AKAP = state[20];
  AKAR4 = state[21];
  AKAR4_C = state[22];
  AKAR4p = state[23];
  
  # Parameters
  kf_Rii_C__RiiP_C <- parameters[1];
  kf_RiiP_CxcAMP__RiiP_C_cAMP <- parameters[2];
  kb_RiiP_CxcAMP__RiiP_C_cAMP <- parameters[3];
  kf_RiiP_cAMPxC__RiiP_C_cAMP <- parameters[4];
  kb_RiiP_cAMPxC__RiiP_C_cAMP <- parameters[5];
  kb_RiiPXcAMP__RiiP_cAMP <- parameters[6];
  kf_RiiPXcAMP__RiiP_cAMP <- parameters[7];
  kf_RiiPxC__RiiP_C <- parameters[8];
  kb_RiiPxC__RiiP_C <- parameters[9];
  kf_cAMPxRii__Rii_cAMP <- parameters[10];
  kb_cAMPxRii__Rii_cAMP <- parameters[11];
  kf_Rii_CxcAMP__Rii_C_cAMP <- parameters[12];
  kb_Rii_CxcAMP__Rii_C_cAMP <- parameters[13];
  kf_RiixC__Rii_C <- parameters[14];
  kf_Rii_cAMPxC__Rii_C_cAMP <- parameters[15];
  kb_Rii_cAMPxC__Rii_C_cAMP <- parameters[16];
  kf_Rii_C_cAMP__RiiP_C_cAMP <- parameters[1];
  kb_RiixC__Rii_C <- parameters[24];
  # AKAPoff params
  kf_RiiP_cAMP_CaN__CaNXRii_cAMP <- parameters[18];
  kf_RiiPxCaN__RiiP_CaN <- parameters[19];
  kb_RiiPxCaN__RiiP_CaN <- parameters[20];
  kf_RiiP_CaN__RiixCaN <- parameters[18];
  kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN <- parameters[19];
  kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN <- parameters[20];
  #AKAPon params
  AKAPon_1 <- parameters[21];
  AKAPon_2 <- parameters[22];
  AKAPon_3 <- parameters[23];
  #AKAR params
  kf_C_AKAR4 <- parameters[25];
  kb_C_AKAR4 <- parameters[26];
  kcat_AKARp <- parameters[27];
  # Km's OFF & ON
  kmOFF <- parameters[28];
  kmON <- parameters[29];
  KD_T <- parameters[30];
  
  
  # These KDs are changed depending on 'b_AKAP' value: if b_AKAP = 1 -> signal AKAPon; otherwise -> signal AKAPoff
  # Sidenote: b_AKAP is a specie (added by Jo?o in R) whose purpose is solely to indicate if simulation is to be run in AKAPon (1) or AKAPoff (0)
  kf_RiiP_cAMP_CaN__CaNXRii_cAMP = b_AKAP*AKAPon_1  +  (1 - b_AKAP)*kf_RiiP_cAMP_CaN__CaNXRii_cAMP;
  kf_RiiPxCaN__RiiP_CaN = b_AKAP*AKAPon_2  +  (1 - b_AKAP)*kf_RiiPxCaN__RiiP_CaN;
  kb_RiiPxCaN__RiiP_CaN = b_AKAP*AKAPon_3  +  (1 - b_AKAP)*kb_RiiPxCaN__RiiP_CaN;
  kf_RiiP_CaN__RiixCaN = b_AKAP*AKAPon_1  +  (1 - b_AKAP)*kf_RiiP_CaN__RiixCaN;
  kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP*AKAPon_2  +  (1 - b_AKAP)*kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN;
  kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP*AKAPon_3  +  (1 - b_AKAP)*kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN;
    
  
  # Force Km - much like thermodynamical constraints
  kf_RiiPxCaN__RiiP_CaN = b_AKAP *((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON )+
    (1 - b_AKAP)*(kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmOFF;
  kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP *((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON)  + 
    (1 - b_AKAP)* (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF;
  
  # Force Taylor's KD
  kb_RiiP_CxcAMP__RiiP_C_cAMP = kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T;
  
  #Fluxes
  reaction_51 = (kf_Rii_C__RiiP_C*Rii_C)
  reaction_14 = (kf_RiiPxC__RiiP_C*RiiP*C)-(kb_RiiPxC__RiiP_C*RiiP_C)
  reaction_12 = (kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP)-(kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP)
  reaction_43 = (kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP)-(kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP)
  reaction_23 = (kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C)-(kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP)
  reaction_78 = (kf_cAMPxRii__Rii_cAMP*cAMP*Rii)-(kb_cAMPxRii__Rii_cAMP*Rii_cAMP)
  reaction_56 = (kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP)-(kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP)
  reaction_76 = (kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C)-(kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP)
  reaction_62 = (kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP)
  reaction_58 = (kf_RiixC__Rii_C*Rii*C)-(kb_RiixC__Rii_C*Rii_C)
  reaction_44 = (kf_RiiPxCaN__RiiP_CaN*RiiP*CaN)-(kb_RiiPxCaN__RiiP_CaN*RiiP_CaN)
  reaction_33 = (kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP)-(kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN)
  reaction_48 = (kf_RiiP_CaN__RiixCaN*RiiP_CaN)
  reaction_37 = (kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN)
  reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
  reaction_2 = kcat_AKARp*AKAR4_C
  
  #Change concentration of species 
  vf_ <- vector(len = 23)
  vf_[1] = (-reaction_78 - reaction_58 + reaction_48)
  vf_[2] = (-reaction_12 - reaction_43 - reaction_78 - reaction_56)
  vf_[3] = (-reaction_14 - reaction_43 - reaction_44)
  vf_[4] = (-reaction_51 - reaction_56 + reaction_58)
  vf_[5] = (reaction_43 - reaction_23 - reaction_33)
  vf_[6] = (reaction_51 + reaction_14 - reaction_12)
  vf_[7] = (reaction_12 + reaction_23 + reaction_62)
  vf_[8] = (-reaction_14 - reaction_23 - reaction_76 - reaction_58)
  vf_[9] = (reaction_78 - reaction_76 + reaction_37)
  vf_[10] = (reaction_56 + reaction_76 - reaction_62)
  vf_[12] = (-reaction_44 - reaction_33 + reaction_48 + reaction_37)
  vf_[13] = (reaction_44 - reaction_48)
  vf_[14] = (reaction_33 - reaction_37)
  vf_[21] = -reaction_1;
  vf_[22] = reaction_1 - reaction_2;
  vf_[23] = reaction_2
  
  list( vf_ )
  
}

