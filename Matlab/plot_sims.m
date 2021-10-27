function out=plot_sims(post_file,sims_file)

close all;
SIM_MAX=0.2;
const_flag=0;
scale=1000;

nWT=18;%Number of WT experiments
nMut=6;%Number of mutation experiments

[Y_exp, exp_name, data_names2] = load_targets();
 

%%%%%%%%%%%%%% Set up simulations for parameter estimation 
L=load(sims_file);

[nIt,~,nPar]=size(L.sims.output{1});

for i=1:nWT+nMut
times{i}=L.sims.times{i};
end

AKAR4p_idx1=19; %AKAR4p idx in small model
AKAR4p_idx2=30; %AKAR4p idx in extended model

Y_sims=cell(nWT+nMut);
for i=1:nWT+nMut
    if i<=nWT
        Y_sims{i}=L.sims.output{i}(:,:,AKAR4p_idx1);
    else
        Y_sims{i}=L.sims.output{i}(:,:,AKAR4p_idx2);
    end
end

 

%%%%% Normalize between 0 and 1
exp_min=100; %zero phosphorylation
exp_max=171.6; %full phosphorylation
sim_min=0; %zero phosphorylation
sim_max=SIM_MAX; %full phosphorylation

Y_exp_n=Y_exp;
Y_sims_n=Y_sims;
for i=1:(nWT+nMut)
    Y_exp_n{i}=(Y_exp{i}-exp_min)./(exp_max-exp_min);
    for j=1:nIt
        Y_sims_n{i}(j,:)=(Y_sims{i}(j,:)-sim_min)./(sim_max-sim_min);
    end
end

 %%%%Calculate scores
 scores = get_scores(Y_sims_n,Y_exp_n, 0, 1, 0, 1, nWT,nMut);


%%%%%Classify
class_WT=false(1,nIt);
for j=1:nIt
     class_WT(j)=all(scores(1:18,j)<10); %The full posterior is kept (threshold from ABC used)
end


class_muts_good=false(1,nIt);
for j=1:nIt
     class_muts_good(j)=all(scores(19:24,j)<0.01);
end


class_muts_bad=and(~class_muts_good,class_WT);


%Form KD parameters
KdT=readtable('KD_idx_names');
KdNames=table2array(KdT(:,1));
Kdf_idx=table2array(KdT(:,2));
Kdb_idx=table2array(KdT(:,3));
nKd=length(KdNames);
KdNames(nKd+1)={'Km_OFF'};
KdNames(nKd+2)={'Km_ON'};


%Load ABC parameters and prior
data_post=load(post_file);
data_post.samples = add_constr_to_posterior(data_post.samples, const_flag);
ranges=get_ranges('pkaParms_restri2.txt', scale);
data_prior.samples=get_prior(ranges,length(data_post.samples), const_flag);
data_prior.samples =10.^data_prior.samples;
data_prior.samples = add_constr_to_posterior(data_prior.samples, const_flag);


Kd=NaN(nIt,nKd+2);
Kd_prior=NaN(nIt,nKd+2);

for i=1:nKd
    Kd(:,i)=data_post.samples(1:nIt,Kdb_idx(i))./data_post.samples(1:nIt,Kdf_idx(i));
    Kd_prior(:,i)=data_prior.samples(1:nIt,Kdb_idx(i))./data_prior.samples(1:nIt,Kdf_idx(i));
end

Kd(:,i+1)=data_post.samples(1:nIt,28);%Km_OFF
Kd(:,i+2)=data_post.samples(1:nIt,29);%Km_ON
Kd_prior(:,i+1)=data_prior.samples(1:nIt,28);%Km_OFF
Kd_prior(:,i+2)=data_prior.samples(1:nIt,29);%Km_ON

%Normilize for plotting (between 1 and 1.7)
fp=1.7167;%fully phosphorylated
bl=1.00;%base level
y_exp_n_o=Y_exp_n;
 for i=1:(nWT+nMut)
     Y_exp_n{i}=(Y_exp_n{i}.*(fp-bl))+bl;
      for j=1:nIt       
         Y_sims_n{i}(j,:)=(Y_sims_n{i}(j,:).*(fp-bl))+bl;
      end
  end
 


 %%%%%%%%%%%%%%%%%%%%% FIGURE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure()
  class=class_WT; 
  exp_indx=[7 8 9];
  clrs=[ 0 0 0; 65 105 225; 204 0 0]./255;
for i=1:3
    subplot(1,3,i)
    y_sims=Y_sims_n{exp_indx(i)};
    y_exp=Y_exp_n{exp_indx(i)};
    t=times{exp_indx(i)};
    p1=plot(t, y_sims(class,:),'linewidth',0.1,'Color',[160 160 160]./255);
    hold on
    p2=plot(t, y_exp, 'ko','MarkerFaceColor',clrs(i,:),'MarkerEdgeColor',clrs(i,:));
    ylim([1 1.8])
    xlim([0 605])
    title(data_names2{exp_indx(i)});
    ylabel('Emission ratio (Y/C)')
    xlabel('Time (seconds)')
    legend([p1(1) p2],{ 'simulations', 'experiment'},'Location','southeast')
 end
 h=gcf;
 h.OuterPosition=[100 100 1000 300];
 


 saveas(gcf,'./figures/parameter_estimation1','png');
 
 %%%%%%%%%%%%%%%%%%%%% FIGURE 1B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure()
  class=class_WT; 
  exp_indx=[10 11 12];
  clrs=[ 0 0 0; 65 105 225; 204 0 0]./255;
for i=1:3
    subplot(1,3,i)
    y_sims=Y_sims_n{exp_indx(i)};
    y_exp=Y_exp_n{exp_indx(i)};
    t=times{exp_indx(i)};
    p1=plot(t, y_sims(class,:),'linewidth',0.1,'Color',[160 160 160]./255);
    hold on
    p2=plot(t, y_exp, 'ko','MarkerFaceColor',clrs(i,:),'MarkerEdgeColor',clrs(i,:));
    ylim([1 1.8])
    xlim([0 605])
    title(data_names2{exp_indx(i)});
    ylabel('Emission ratio (Y/C)')
    xlabel('Time (seconds)')
    legend([p1(1) p2],{ 'simulations', 'experiment'},'Location','southeast')
 end
 h=gcf;
 h.OuterPosition=[100 100 1000 300];
 


 
 %%%%%% FIGURE 1C %%%%%%%%%%%%%%%%%%%%%%%%
 
 figure()
  class=class_WT; 
  
  clrs=[ 0 0 0; 65 105 225; 204 0 0]./255;
for i=1:3
    subplot(1,3,i)
    exp_indx=[7 8 9];
    y_sims=Y_sims_n{exp_indx(i)};
    y_exp=Y_exp_n{exp_indx(i)};
    t=times{exp_indx(i)};
    plot(t, mean(y_sims(class,:)), 'k--','LineWidth',3);
    hold on
    exp_indx=[ 4 5 6];
    y_sims=Y_sims_n{exp_indx(i)};
    y_exp=Y_exp_n{exp_indx(i)};
    t=times{exp_indx(i)};
    plot(t, mean(y_sims(class,:)), 'r--','LineWidth',3);
    ylim([1 1.8])
    xlim([0 605])
    title(data_names2{exp_indx(i)});
    ylabel('Emission ratio (Y/C)')
    xlabel('Time (seconds)')
    legend([p1(1) p2],{ 'simulations', 'experiment'},'Location','southeast')
 end
 h=gcf;
 h.OuterPosition=[100 100 1000 300];
 

 
 
 %%%%%% FIGURE 2 %%%%%%%%%%%%%%%%%%%%%%%%
 figure()
 class=class_muts_good;
 exp_indx=[ 13 16 18 14 15 17];
 exp_indx=exp_indx+6;
 clrs=[ 102 178 255;  255 0 0; 0 204 102;  0 0 153; 153 0 0; 51 102 0]./255;

for i=1:6
    subplot(2,3,i)
    y_sims=Y_sims_n{exp_indx(i)};    
    y_exp=Y_exp_n{exp_indx(i)};
    t=times{exp_indx(i)};
    p1=plot(t, y_sims(class,:),'linewidth',0.1,'Color',[160 160 160]./255);
    hold on
    p2=plot(t, y_exp, 'ko','MarkerFaceColor',clrs(i,:),'MarkerEdgeColor',clrs(i,:));
    ylim([1 1.8])
    xlim([0 605])
    title(data_names2(exp_indx(i)));
    ylabel('Emission ratio (Y/C)')
    xlabel('Time (seconds)')
    legend([p1(1) p2],{ 'simulations', 'experiment'},'Location','southeast')
 end
 h=gcf;
 h.OuterPosition=[100 100 1200 600];
 
 saveas(gcf,'./figures/mutations','png');

%%%%%%%%%%%%%%%%%%%%%% FIGURE 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()

mv=1;%0.005
sizep=1;%30

subplot(2,3,1)
Kd_idx=[1 2 3 4 5 6 7 8];
boxplot(log10(Kd_prior(:,Kd_idx)),KdNames(Kd_idx));
xtickangle(45);
ylim([-8 8])
title('Prior');
ylabel('log10(KD)')

subplot(2,3,2)
boxplot(log10(Kd(class_WT,Kd_idx)),KdNames(Kd_idx));
xtickangle(45);
ylim([-8 8])

title('Restricted by WT data');


subplot(2,3,3)
boxplot(log10(Kd(class_muts_good,Kd_idx)),KdNames(Kd_idx));
xtickangle(45);
ylim([-8 8])


title('Restricted by WT and mutation data');

 

mv=3/8;%0.005
sizep=0.5;%30

%Kd_idx=[4 5];
Kd_idx=[2 3];
%Kd_idx=[5 8];
for i=1:2
 subplot(2,3,3+i)
 histogram(log10(Kd(class_muts_bad,Kd_idx(i))),20, 'facecolor', 'r');
 hold on;
 histogram(log10(Kd(class_muts_good,Kd_idx(i))),20,'facecolor', 'b');
 xlabel(strcat('log10(',KdNames(Kd_idx(i)),')'), 'Interpreter','None');
 legend({ 'Far away'; 'Close to data'});
end
i=1;
j=2;
subplot(2,3,6)
scatter(log10(Kd(class_muts_bad,Kd_idx(i))),log10(Kd(class_muts_bad,Kd_idx(j))),'r','filled','MarkerFaceAlpha',mv,'MarkerEdgeAlpha',mv)
hold on;
scatter(log10(Kd(class_muts_good,Kd_idx(i))),log10(Kd(class_muts_good,Kd_idx(j))),'b','filled','MarkerFaceAlpha',mv,'MarkerEdgeAlpha',mv)
legend({'Far away';'Close to data';});

xlabel(strcat('log10(',KdNames(Kd_idx(i)),')'),'Interpreter','None');
ylabel(strcat('log10(',KdNames(Kd_idx(j)),')'),'Interpreter','None');


 h=gcf;
 h.OuterPosition=[100 100 1200 800];

 saveas(gcf,'./figures/KDs','png');



%%%%%%%%%%%%%%%%%%%%%% FIGURE 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

class_muts_bad_red=class_muts_bad(1:10:end); %reduce number of lines that are plotted (otherwise it takes to long time to plot)
class_muts_good_red=class_muts_good(1:1:end);

exp_idx=[7 8 9 21 22 23 24];
spec_idx=[19 3 1 4 6 5 7 9 10 8 2]; %relative to "names" below
AKAR4p_idx2=[NaN 30;3 20;1 16; 4 19; 6 21; 5 22;7 23;9 18; 10 17; NaN 8; NaN 2];
ylims=[2 7 7 1 1 1 0.5 1 1 1 2 10 10];


S=sbioloadproject('../models/PKA_insee_05_06_20_AKAR');
modelobj = S.m1;
[t, ya, names] = sbiosimulate(modelobj);

names={'Rii','cAMP','RiiP','Rii-C',{'RiiP-','cAMP'},'RiiP-C',{'RiiP-C-','cAMP'},...
     'C',{'Rii-','cAMP'},{'Rii-C-','cAMP'},'Total RII','CaN','RiiP-CaN',...
 {'RiiP-','cAMP_CaN'},'Total C','Free C','AKAR4','AKAR4-C',{'Emission', 'ratio (Y/C)'},'AKAR4 ratio'};


nWT=7;
nSpec=11;


figure();
cnt=1;
for j=1:nSpec
    for i=1:nWT 
        y_sims=L.sims.output{exp_idx(i)};
        t=times{exp_idx(i)};
        if i<4 %Experiment 1-3 (small model)
             y_sims_sub=y_sims(:,:,spec_idx(j));
             if j==1
                 y_sims_sub_n01=(y_sims_sub-sim_min)/(sim_max-sim_min);
                 y_sims_sub=y_sims_sub_n01*(fp-bl)+bl;
             end
        elseif i>=4 % Exp 15-18, extended model
            if ~isnan(AKAR4p_idx2(j,1)) %alfa and beta type, needs to be summed together 
                y_sims_sub=y_sims(:,:,AKAR4p_idx2(j,1))+y_sims(:,:,AKAR4p_idx2(j,2));
                 if j==1
                    y_sims_sub_n01=(y_sims_sub-sim_min)/(sim_max-sim_min);
                    y_sims_sub=y_sims_sub_n01*(fp-bl)+bl;
                 end
        elseif isnan(AKAR4p_idx2(j,1)) %not alfa and beta, single type
                y_sims_sub=y_sims(:,:,AKAR4p_idx2(j,2));
                 if j==1
                    y_sims_sub_n01=(y_sims_sub-sim_min)/(sim_max-sim_min);
                    y_sims_sub=y_sims_sub_n01*(fp-bl)+bl;
                 end
            end
        end
        
        subplot(nSpec,nWT,cnt);
        yb=y_sims_sub(class_muts_bad,:,:)';
        yg=y_sims_sub(class_muts_good,:,:)';
        plot(t, yb(:,1:70:end),'linewidth',0.1,'Color',[255 0 0]/255);
        hold on
        plot(t, yg(:,1:10:end) ,'linewidth',0.1,'Color',[0 0 255]/255);
        if j==1
            hold on
            plot(times{exp_idx(i)}, Y_exp_n{exp_idx(i)}, 'ko-','MarkerFaceColor',[0,0,0],'MarkerSize',3);
            title(data_names2{exp_idx(i)});
            
        end
        if i==1 
            ylabel(names{spec_idx(j)},'Interpreter','none')
        end
        if j<nSpec
            set(gca,'XTickLabel',[]);
        end
        if j==nSpec
           xlabel('time (s)');
        end
        sgtitle('Blue=close to mutation data, red=far away');
        ylim([0 ylims(j)]);
        if j==1
            ylim([1 1.8])
        end
        xlim([0 605])
        cnt=cnt+1;
    end   
end
h=gcf;
h.OuterPosition=[100 100 1000 3000];
saveas(gcf, './figures/species_posterior_classified','png');
end
