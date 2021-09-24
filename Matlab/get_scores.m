function [score] = get_scores(Y_sims,Y_exp, sim_min,sim_max, exp_min, exp_max, nWT,nMut)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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
     %nT=length(L.sims.times{i});
     for j=1:nSamples        
        %score(i,j)=1/nT.*sqrt(sum((Y_exp_n{i}-Y_sims_n{i}(j,:)).^2));
        score(i,j)=1/nT.*sum((Y_exp_n{i}-Y_sims_n{i}(j,:)).^2);
     end
 end

end

