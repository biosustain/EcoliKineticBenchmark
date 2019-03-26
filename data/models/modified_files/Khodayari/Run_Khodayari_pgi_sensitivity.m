clear all;
close all;

%%%
% 6 simulations for pgi - WT, pgi(0) dpgi,  0.2*pgi, pgi(20) 1.2*pgi, pgi(50) 2.4*pgi, pgi(100) 4.1*pgi

enzyme_id = 67; % ID of pgi
condition = 'aerobic_glucose';

sample_id = 'WT';
modification = 1;
fprintf('Working on %s\n',sample_id)
[vnet, conc] = Main_Module(enzyme_id,modification,condition);
Vnet = real(vnet); Conc = real(conc);
save(sprintf('khodayari_pgi_sens_%s.mat',sample_id),'Vnet','Conc');


sample_id = 'dpgi';
modification = 0;
fprintf('Working on %s\n',sample_id)
[vnet, conc] = Main_Module(enzyme_id,modification,condition);
Vnet = real(vnet); Conc = real(conc);
save(sprintf('khodayari_pgi_sens_%s.mat',sample_id),'Vnet','Conc');

sample_id = 'pgi(0)';
modification = 0.2;
fprintf('Working on %s\n',sample_id)
[vnet, conc] = Main_Module(enzyme_id,modification,condition);
Vnet = real(vnet); Conc = real(conc);
save(sprintf('khodayari_pgi_sens_%s.mat',sample_id),'Vnet','Conc');

sample_id = 'pgi(20)';
modification = 1.2;
fprintf('Working on %s\n',sample_id)
[vnet, conc] = Main_Module(enzyme_id,modification,condition);
Vnet = real(vnet); Conc = real(conc);
save(sprintf('khodayari_pgi_sens_%s.mat',sample_id),'Vnet','Conc');

sample_id = 'pgi(50)';
modification = 2.4;
fprintf('Working on %s\n',sample_id)
[vnet, conc] = Main_Module(enzyme_id,modification,condition);
Vnet = real(vnet); Conc = real(conc);
save(sprintf('khodayari_pgi_sens_%s.mat',sample_id),'Vnet','Conc');

sample_id = 'pgi(100)';
modification = 4.1;
fprintf('Working on %s\n',sample_id)
[vnet, conc] = Main_Module(enzyme_id,modification,condition);
Vnet = real(vnet); Conc = real(conc);
save(sprintf('khodayari_pgi_sens_%s.mat',sample_id),'Vnet','Conc');

return