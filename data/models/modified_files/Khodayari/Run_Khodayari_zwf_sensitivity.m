clear all;
close all;
tic;

% 3 simulations - WT, zwf * 15, zwf KO
enzyme_id = 42; % ID of zwf
condition = 'aerobic_glucose';

sample_id = 'WT';
modification = 1;
fprintf('Working on %s\n',sample_id)
[vnet, conc] = Main_Module(enzyme_id,modification,condition);
Vnet = real(vnet); Conc = real(conc);
save(sprintf('khodayari_zwf_sens_%s.mat',sample_id),'Vnet','Conc');

sample_id = 'dzwf';
modification = 0;
fprintf('Working on %s\n',sample_id)
[vnet, conc] = Main_Module(enzyme_id,modification,condition);
Vnet = real(vnet); Conc = real(conc);
save(sprintf('khodayari_zwf_sens_%s.mat',sample_id),'Vnet','Conc');

sample_id = 'zwf(15)';
modification = 15;
fprintf('Working on %s\n',sample_id)
[vnet, conc] = Main_Module(enzyme_id,modification,condition);
Vnet = real(vnet); Conc = real(conc);
save(sprintf('khodayari_zwf_sens_%s.mat',sample_id),'Vnet','Conc');
