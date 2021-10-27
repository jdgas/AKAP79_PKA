function draws=get_prior(ranges,nSamples, const_flag)
%  Functions for making simulated trajectories from posterior distribution
%  sample and reproducing figures in eLife 2021;10:e68164 
 
%  Copyright (C) 2021 Olivia Eriksson (olivia@kth.se) and Parul Tewatia (parul.tewatia@scilifelab.se)

%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License or
%  (at your option) any later version.
% 
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%  GNU General Public License for more details.

    draws=NaN(nSamples,length(ranges));
    for i=1:length(ranges)
        a=ranges(i,1);
        b=ranges(i,2);
        r = a + (b-a).*rand(nSamples,1);
        
        draws(:,i)=r;
       
    end
    
end