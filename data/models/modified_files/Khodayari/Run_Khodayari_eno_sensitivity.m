clear all;
close all;

%%% eno sensitivity
%% eno(0) eno*0.2, eno(50) eno*1.8, eno(200) eno*3.0, eno(500) eno*3.1

enzyme_id = 17; % ID of eno
condition = 'aerobic_glucose';

sample_id = 'WT';
modification = 1;
fprintf('Working on %s\n',sample_id)
[vnet, conc] = Main_Module(enzyme_id,modification,condition);
Vnet = real(vnet); Conc = real(conc);
save(sprintf('khodayari_eno_sens_%s.mat',sample_id),'Vnet','Conc');


sample_id = 'eno(0)';
modification = 0.2;
fprintf('Working on %s\n',sample_id)
[vnet, conc] = Main_Module(enzyme_id,modification,condition);
Vnet = real(vnet); Conc = real(conc);
save(sprintf('khodayari_eno_sens_%s.mat',sample_id),'Vnet','Conc');

sample_id = 'eno(50)';
modification = 1.8;
fprintf('Working on %s\n',sample_id)
[vnet, conc] = Main_Module(enzyme_id,modification,condition);
Vnet = real(vnet); Conc = real(conc);
save(sprintf('khodayari_eno_sens_%s.mat',sample_id),'Vnet','Conc');

sample_id = 'eno(200)';
modification = 3.0;
fprintf('Working on %s\n',sample_id)
[vnet, conc] = Main_Module(enzyme_id,modification,condition);
Vnet = real(vnet); Conc = real(conc);
save(sprintf('khodayari_eno_sens_%s.mat',sample_id),'Vnet','Conc');

sample_id = 'eno(500)';
modification = 3.1;
fprintf('Working on %s\n',sample_id)
[vnet, conc] = Main_Module(enzyme_id,modification,condition);
Vnet = real(vnet); Conc = real(conc);
save(sprintf('khodayari_eno_sens_%s.mat',sample_id),'Vnet','Conc');

return