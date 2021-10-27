function param_sample = add_constr_to_posterior(param_sample, const_flag)
%param_sample is a nSamplesxnPars matrix with the posterior from ABC
%const_flag=0 Taylor Kd and Km constraints
%const_flag=1 Taylor Kd and Km and thermodynamic constraints

if const_flag>=0
%Add Km
param_sample(:,19)=(param_sample(:,20)+param_sample(:,18))./(param_sample(:,28)); %OFF parameters
param_sample(:,22)=(param_sample(:,23)+param_sample(:,21))./param_sample(:,29); %ON parameters
%Idx 28 and 29 is the Km value  

%Add Taylor (A paper reference Taylor):
param_sample(:,3)=param_sample(:,2).*param_sample(:,30);% KD12=0.7  (0.35 1.4)
end
%Add thermodynamic constr
if const_flag==1
      t_idx=[10, 11, 14, 12, 16, 24, 13, 15; %First thermodyn constr
            7, 6, 8, 2, 5, 9, 3, 4];% Kf34 %Second thermodyn constr
    for i=1:2
        param_sample(:,t_idx(i,1))=...
        param_sample(:,t_idx(i,2)).*param_sample(:,t_idx(i,3)).*param_sample(:,t_idx(i,4)).*param_sample(:,t_idx(i,5))./...
        (param_sample(:,t_idx(i,6)).*param_sample(:,t_idx(i,7)).*param_sample(:,t_idx(i,8)));
    end
end

end
