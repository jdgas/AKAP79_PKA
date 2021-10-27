function [score] = get_scores(Y_sims,Y_exp, sim_min,sim_max, exp_min, exp_max, nWT,nMut)
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
Y_exp_n=Y_exp;
Y_sims_n=Y_sims;

[nSamples,~]=size(Y_sims{1});

for i=1:(nWT+nMut)
    Y_exp_n{i}=(Y_exp{i}-exp_min)./(exp_max-exp_min);
    for j=1:nSamples
        Y_sims_n{i}(j,:)=(Y_sims{i}(j,:)-sim_min)./(sim_max-sim_min);
    end
end

score=NaN(nWT+nMut,nSamples);
 for i=1:(nWT+nMut)
     [~,nT]=size(Y_sims{i});
     for j=1:nSamples        
        score(i,j)=1/nT.*sum((Y_exp_n{i}-Y_sims_n{i}(j,:)).^2);
     end
 end

end

