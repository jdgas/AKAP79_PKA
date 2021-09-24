function [Y_exp, exp_name, exp_name2] = load_targets()
nWT=18;
nMut=6;
exp_data_dir='Targets/';
exp_data_files={'TS_0cAMP.txt','TS_02cAMP.txt','TS_1cAMP.txt', 'TS_2cAMP.txt', 'Calib_Curves.txt','mutant.mat'};

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

