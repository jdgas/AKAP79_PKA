function [Y_exp, exp_name, exp_name2] = load_targets()
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
nWT=18;
nMut=6;
exp_data_dir='../Targets/';
exp_data_files={'0_cAMP.txt','0.2_cAMP.txt','1_cAMP.txt', '2_cAMP.txt', 'Calib_Curves.txt','mutant.mat'};

exp_name={'NoCNnoAKAP cAMP=0', 'CNonly cAMP=0', 'CNandAKAP cAMP=0', 'NoCNnoAKAP cAMP=0.2', 'CNonly cAMP=0.2', 'CNandAKAP cAMP=0.2','NoCNnoAKAP cAMP=1', 'CNonly cAMP=1', 'CNandAKAP cAMP=1','NoCNnoAKAP cAMP=2', 'CNonly cAMP=2', 'CNandAKAP cAMP=2', 'Calib C=0.4', 'Calib C=0.2', 'Calib C=0.1', 'Calib C=0.05', 'Calib C=0.025','Calib C=0.0','WT','WT79CaN','98A+79CN','98A alone','98E+79CN', '98E alone'};
exp_name2={' ??','??','??',' ??','??','??','PKA(II)', 'PKA(II)m+CN', {'PKA(II)+CN', '+ AKAP79_{c97}'},' ??','??','??','??','??',' ??','??','??','??','RII\alpha WT no CN','RII\alpha WT plus CN', 'RII\alpha S98A plus CN', 'RII\alpha S98A no CN', 'RII\alpha S98E plus CN','RII\alpha S98E no CN'};


data_idx=[1 2; 1 3; 1 4;...
          2 2; 2 3; 2 4;...
          3 2; 3 3; 3 4;...
          4 2; 4 3; 4 4;...
          5 2; 5 3; 5 4; 5 5; 5 6; 5 7;... 
          6 2; 6 3; 6 4; 6 5; 6 6; 6 7]; %file, column


Y_exp=cell(1,nWT+nMut);
nc=1.08;
 

for i=1:(nWT+nMut)
        filename=[exp_data_dir exp_data_files{data_idx(i,1)}];
        if (i<=12)
        data=table2array(readtable(filename));
        Y_exp{i}=data(:,data_idx(i,2))';
        elseif i>12 && i<=18
        data=table2array(readtable(filename));
        Y_exp{i}=data(:,data_idx(i,2))'./nc;
        elseif i>18
            T=load(filename);
            data=T.mutant{2:end,data_idx(i,2)};
            Y_exp{i}=data'./nc;
        end
end

end

