function out=simulations_from_posterior_sample(nIt,post_file, sims_file,const_flag)
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

%sims_file:name of file where to store simulation data.
%post_file: file with the posterior parameter distribution


out=0; 
AKARtot=0.2;
data=load(post_file);
S=sbioloadproject('../models/PKA_insee_05_06_20_AKAR');  %load model
modelobj = S.m1;
data.samples=add_constr_to_posterior(data.samples,const_flag);

[nSamples,~]=size(data.samples);
nExp=18;
nMut=6;
nExp_all=nExp+nMut;

output=cell(1,nExp_all);

%%%%%%%%%%%%%%%%SIMULATIONS SMALL MODEL - experiments%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
configsetObj = getconfigset(modelobj);
set(configsetObj, 'SolverType','ode15s');
%set(configsetObj.SolverOptions, 'AbsoluteTolerance', 1.0e-8);
%set(configsetObj.SolverOptions, 'RelativeTolerance', 1.0e-6);
set(configsetObj.SolverOptions, 'AbsoluteTolerance', 1.0e-10);
set(configsetObj.SolverOptions, 'RelativeTolerance', 1.0e-8);
set(configsetObj, 'MaximumWallClock', 300)
set(configsetObj.RuntimeOptions, 'StatesToLog', 'all');
nSpecs=length(configsetObj.RuntimeOptions.StatesToLog);

offParIdx=[1:16 1 18:20 18:20 24 31 25:27];
onParIdx=[1:16 1 21:23 21:23 24 31 25:27];

simtime = [605 605 605 605 605 605 605 605 605 605 605 605 1105 1105 1105 1105 1105 1105];
is_on=[false false true false false true false false true false false true false false false false false false];

species_changed=[6.30,0,0,0.63,0,6.93,0,AKARtot;...
                 6.30,0,0,0.63,0,6.93,1.5,AKARtot;...
                 6.30,0,0,0.63,0,6.93,1.5,AKARtot;...
                 6.30,0.2,0,0.63,0,6.93,0,AKARtot;...
                 6.30,0.2,0,0.63,0,6.93,1.5,AKARtot;...
                 6.30,0.2,0,0.63,0,6.93,1.5,AKARtot;...
                 6.30,1,0,0.63,0,6.93,0,AKARtot;...
                 6.30,1,0,0.63,0,6.93,1.5,AKARtot;...
                 6.30,1,0,0.63,0,6.93,1.5,AKARtot;...
                 6.30,2,0,0.63,0,6.93,0,AKARtot;...
                 6.30,2,0,0.63,0,6.93,1.5,AKARtot;...
                 6.30,2,0,0.63,0,6.93,1.5,AKARtot;...
                 0,0,0,0,0.4,0,0,AKARtot;
                 0,0,0,0,0.2,0,0,AKARtot;
                 0,0,0,0,0.1,0,0,AKARtot;
                 0,0,0,0,0.05,0,0,AKARtot;
                 0,0,0,0,0.025,0,0,AKARtot;
                 0,0,0,0,0,0,0,AKARtot];
             % Rii, cAMP, RiiP, Rii_C, C, total_Rii, CaN, AKAR

data.samples(:,31)=modelobj.parameters(25).Value*ones(nSamples,1);
AKAPoffPar= data.samples(:,offParIdx); %using AKAPoff parameters
AKAPonPar= data.samples(:,onParIdx);%using AKAPon parameters

modelobj.Rules(1).Active=1;%All Rii
modelobj.Rules(2).Active=1;%Total C
modelobj.Rules(3).Active=0;%cAMP
modelobj.Rules(4).Active=1;%Total_Rii
modelobj.Rules(5).Active=1;%AKAR ratio

%Needed in case simualtions are run in parallel
for i=1:nIt
    modelarray(i)=modelobj;
end

times=cell(1,nExp);
for k=1:nExp
    fprintf('Datafit exp k=%i\n',k);   
    times{k}=5:5:simtime(k);
    outs=NaN(nIt,length(times{k}),nSpecs);
    parfor i=1:nIt    
        configsetObj = getconfigset(modelarray(i));
        set(configsetObj, 'StopTime', times{k}(end));
        set(configsetObj.SolverOptions, 'OutputTimes', times{k});
        
        modelarray(i).species(1).InitialAmount=species_changed(k,1);%Rii
        modelarray(i).species(2).InitialAmount=species_changed(k,2);%cAMP
        modelarray(i).species(3).InitialAmount=species_changed(k,3);%RiiP
        modelarray(i).species(4).InitialAmount=species_changed(k,4);%Rii_C
        modelarray(i).species(8).InitialAmount=species_changed(k,5);%C
        modelarray(i).species(11).InitialAmount=species_changed(k,6);%total Rii
        modelarray(i).species(12).InitialAmount=species_changed(k,7);%CaN
        modelarray(i).species(17).InitialAmount=species_changed(k,8);%AKAR
        
        if is_on(k)
            for j=1:length(AKAPonPar(1,:))
                modelarray(i).parameters(j).Value=AKAPonPar(i,j);  
            end
        else
            for j=1:length(AKAPonPar(1,:))
                modelarray(i).parameters(j).Value=AKAPoffPar(i,j);  
            end
        end
        
        try
            [~, ya,~] = sbiosimulate(modelarray(i));
            outs(i,:,:)=ya;
        catch 
          fprintf('OBS:The simulation of the simple model iteration i=%i failed\n',i); 
        end
        
        
    end 
    output{k}=outs;
end

%%%%%%%%%%%%%%%%%%%%% SIMULATIONS EXTENDED MODEL - Mutations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S=sbioloadproject('../models/PKA_all_mutant_data_25_10'); %Have both Rii forms (alfa and beta).
modelobj=S.m1;

configsetObj = getconfigset(modelobj);
set(configsetObj, 'SolverType','ode15s');
% set(configsetObj.SolverOptions, 'AbsoluteTolerance', 1.0e-8);
% set(configsetObj.SolverOptions, 'RelativeTolerance', 1.0e-6);
 set(configsetObj.SolverOptions, 'AbsoluteTolerance', 1.0e-10);
 set(configsetObj.SolverOptions, 'RelativeTolerance', 1.0e-8);

set(configsetObj, 'MaximumWallClock', 300)
set(configsetObj.RuntimeOptions, 'StatesToLog', 'all');
nSpecs=length(configsetObj.RuntimeOptions.StatesToLog);

offParIdx=[1:17 18:20 18:20 24 4 5 18 19 20 19 20 18 15 16 14 24 12 13 1 7 6 2 3 8 9 10 11 26 25 27 17];
onParIdx=[1:17 21:23 21:23 24 4 5 21 22 23 22 23 21 15 16 14 24 12 13 1 7 6 2 3 8 9 10 11 26 25 27 17];
species_changed=NaN(nMut,length(modelobj.Species));
species_changed(1,:)=[0.9364,1,0,0.094,0,0,0,0,0,0,1.03,0,0,0,0.63,5.3638,0,0,0.536,0,0,0,0,0,0,5.9,6.93,AKARtot,0,0]; % WT AKAPoff
species_changed(2,:)=[0.9364,1,0,0.094,0,0,0,0,0,0,1.03,1.5,0,0,0.63,5.3638,0,0,0.536,0,0,0,0,0,0,5.9,6.93,AKARtot,0,0]; % WT+CAN 
species_changed(3,:)=[0.9364,1,0,0.094,0,0,0,0,0,0,1.03,1.5,0,0,0.63,5.3638,0,0,0.536,0,0,0,0,0,0,5.9,6.93,AKARtot,0,0]; % 98A+CaN
species_changed(4,:)=[0.9364,1,0,0.094,0,0,0,0,0,0,1.03,0,0,0,0.63,5.3638,0,0,0.536,0,0,0,0,0,0,5.9,6.93,AKARtot,0,0]; % 98A
species_changed(5,:)=[0.4044,1,0,0.6256,0,0,0,0.0035,0,0,1.03,1.5,0,0,0.63,0,0,0,0,5.8991,0.01,0,0,0,0,5.9,6.93,AKARtot,0,0]; %98E+CaN
species_changed(6,:)=[0.4044,1,0,0.6256,0,0,0,0.0035,0,0,1.03,0,0,0,0.63,0,0,0,0,5.8991,0.01,0,0,0,0,5.9,6.93,AKARtot,0,0]; %98E

is_on=[false true true false true false];
AKAPoffPar= data.samples(:,offParIdx); %all AKAPoff parameters
AKAPonPar= data.samples(:,onParIdx);%all AKAPon parameters


for i=1:nIt
    modelarray(i)=modelobj;
end


for k=1:nMut
    
    times{k+nExp}=5:5:605;
    ntp=length(times{k+nExp});
    outsMut=NaN(nIt,ntp,nSpecs);
    fprintf('Mutation k=%i\n',k);
    
    parfor i=1:nIt
    configsetObj = getconfigset(modelarray(i));
    set(configsetObj, 'StopTime', times{k+nExp}(end));
    set(configsetObj.SolverOptions, 'OutputTimes', times{k+nExp});
    if is_on(k)
        for j=1:length(AKAPonPar(1,:))
            modelarray(i).parameters(j).Value=AKAPonPar(i,j);  
        end  
    else
       for j=1:length(AKAPonPar(1,:))
        modelarray(i).parameters(j).Value=AKAPoffPar(i,j);  
       end
    end
    if k==3 || k==4
         modelarray(i).parameters(39).Value=0;
         modelarray(i).parameters(51).Value=0;
    end
    if k==5 || k==6
        modelarray(i).parameters(27).Value=0;
        modelarray(i).parameters(32).Value=0;
    end
    for  j=1:length(species_changed(k,:))
         modelarray(i).species(j).Value=species_changed(k,j);
    end
     
    try
        [~, yb,~] = sbiosimulate(modelarray(i));
        outsMut(i,:,:)=yb;
    catch 
         fprintf('OBS:The simulation of the extended model iteration i=%i failed\n',i); 
    end
        
     
    end
   output{nExp+k}=outsMut; 
end


sims.output=output;
sims.times=times;

save(sims_file,'sims', '-v7.3');
out=1;
end
    
