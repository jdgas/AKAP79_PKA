# Uncertainty Quantification: Distance criterion
# Copyright (C) 2018 
# Alexandra Jauhiainen (alexandra.jauhiainen@gmail.com)
# Olivia Eriksson

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

getMaxScore  <- function(xtarget, ytarget, xx, yy){ # rename: getScore; original Alexandra: getMaxScore
  ms <- 1
  distance <- 1
  x <- xtarget
  y <- ytarget
  scores <- rep(NA, length(x))
  
  if (!(all(yy==0, na.rm=T) | any(yy<0, na.rm=T)| any(is.na(yy))| any(yy>120))){    
    
    # sumSqError
    for (j in 1:length(x)){
      #   sigma <- yN[j]*.05 + maxyN*.05
      sigma <- 1
      distance <- distance + ((yy[j]-ytarget[j])/sigma)^2
    }
    distance <- distance / length(ytarget)
  }
  
  #ms
  distance
}