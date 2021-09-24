% post_file='../../ABCruns/Automatized_code_v2_large_free_k43_k87/DrawsNoThermoScale1000_7-13-19-22-8-14-20-23-9-15-21-24';
% sims_file=[post_file '_SIM'];%Where to save simulations
% %nIt=15351; %Must be less than total number of samples in post_file
% nIt=500;
% const_flag=0; %Should be=0, describes what constraints that are used
% 
% %SIMULATE MODEL (to run the full sammple nIt=15351, takes a few hours on a single core)
% simulations_from_posterior_sample(nIt,post_file, sims_file,const_flag);

%Reproduce article figures
plot_sims(post_file,sims_file);
plot_sims_all_species(sims_file, post_file);

