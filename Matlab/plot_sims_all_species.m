
function class=plot_sims_all_species(sims_file, draws_file)
%close all;

nWT=18;%Number of WT experiments
nMut=6;%Number of mutation experiments

fp=1.7167;%fully phosphorylated
bl=1.00;%base level

%%%% SET UP DATA  %%%%%%%%%%%

 [Y_exp, exp_name] = load_targets();

%%%%% Normalize between 0 and 1 ( 
exp_min=100; %zero phosphorylation
exp_max=171.6; %full phosphorylation

Y_exp_n=Y_exp;

for i=1:(nWT+nMut)
    Y_exp_01=(Y_exp{i}-exp_min)./(exp_max-exp_min); 
    Y_exp_n{i}=Y_exp_01.*(fp-bl)+bl;
    
end

%%%%%%%%%%%%%% Set up simulations for parameter estimation, normilize between 0 and 1 

%contains all simulation trajectories

L=load(sims_file);
[nIt,~,nPar]=size(L.sims.output{1});
AKAR4p_idx1=19; %AKAR4p idx in small model
AKAR4p_idx2=30; %AKAR4p idx in extended model

sim_min=0; %zero phosphorylation
sim_max=0.2; %full phosphorylation

Y_sims_n=cell(nWT+nMut);
for i=1:nWT+nMut
    if i<=nWT
        Y_sims_01=(L.sims.output{i}(:,:,AKAR4p_idx1)-sim_min)./(sim_max-sim_min);  
        
        Y_sims_n{i}=Y_sims_01.*(fp-bl)+bl;
       
         
    else
        Y_sims_01=(L.sims.output{i}(:,:,AKAR4p_idx2)-sim_min)./(sim_max-sim_min);
        Y_sims_n{i}=Y_sims_01.*(fp-bl)+bl;
       
    end
end

scores = get_scores(Y_sims_n,Y_exp_n, bl,fp, bl, fp, nWT,nMut);



scoresSum=sum(scores(1:24,:))./24;

[mm,mIdx]=min(scoresSum);

class=false(1,nIt);
for j=1:nIt
     class(j)=all(scores(19:24,j)<0.01);
end

disp(sum(class))

nth=1; 


plotIdx=[1:12 21:24];
plotLoc=[1:3 5:7 9:11 13:15 17:20];
figure()
for i=1:length(plotIdx)
     subplot(5,4,plotLoc(i))
     times=L.sims.times{plotIdx(i)};
     p1=plot(times,Y_sims_n{plotIdx(i)}(class(1:nth:end),:),'Color',[190 190 190]/255);
     title(exp_name{plotIdx(i)});
     ylabel('Emission ratio');
     xlabel('Time (s)')
     hold on
     p2=plot(times,Y_exp_n{plotIdx(i)},'o', 'MarkerFaceColor',[150 0 0]/255,'MarkerEdgeColor',[150 0 0]/255, 'MarkerSize',5);
     ylim([0.9 1.8])
end
legend([p1(1) p2],{ 'simulations', 'experiment'},'Location','southeast')
sgtitle('Modelling of responses at different cAMP levels')
h=gcf;
h.OuterPosition=[100 100 1000 3000];
saveas(gcf,'./figures/experiments_versus_simulations','png');


figure()
for i=1:length(plotIdx)
     subplot(5,4,plotLoc(i))
     times=L.sims.times{plotIdx(i)};
     p1=plot(times,Y_sims_n{plotIdx(i)}(class(1:nth:end),:),'Color',[190 190 190]/255);
     title(exp_name{plotIdx(i)});
     hold on
     p2=plot(times,Y_sims_n{plotIdx(i)}(mIdx,:),'-','Color',[150 0 0]/255, 'LineWidth',3);
     ylim([0.9 1.8])
     ylabel('Emission ratio');
     xlabel('Time (s)')
end
legend([p1(1) p2],{ 'simulations', 'best score'},'Location','southeast')
sgtitle('Modelling of responses at different cAMP levels')
h=gcf;
h.OuterPosition=[100 100 1000 3000];
saveas(gcf,'./figures/best_score_simulations','png');

end