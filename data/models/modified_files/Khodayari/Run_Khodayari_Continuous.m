clear all;
close all;
tic;

condition = 'aerobic_glucose';
modification = 0;

for i=1:1:18
switch i
    case 1
        sample_id = 'Delta_sucAB';
        enzyme_id = 8;
    
    case 2
        sample_id = 'Delta_fbaAB';
        enzyme_id = 37;

    case 3
        sample_id = 'Delta_zwf';
        enzyme_id = 42;

    case 4
        sample_id = 'Delta_pts';
        enzyme_id = 44;

    case 5
        sample_id = 'Delta_gnd';
        enzyme_id = 51;

    case 6
        sample_id = 'Delta_pfkAB';
        enzyme_id = 66;

    case 7
        sample_id = 'Delta_pgi';
        enzyme_id = 67;

    case 8
        sample_id = 'Delta_pgl';
        enzyme_id = 69;

    case 9
        sample_id = 'Delta_ppsA';
        enzyme_id = 74;

    case 10
        sample_id = 'Delta_pykF';
        enzyme_id = 76;
    
    case 11
        sample_id = 'Delta_rpe';
        enzyme_id = 78;    

    case 12
        sample_id = 'Delta_rpiAB';
        enzyme_id = 79;

    case 13
        sample_id = 'Delta_sdhCD';
        enzyme_id = 82;

    case 14
        sample_id = 'Delta_talAB';
        enzyme_id = 84;

    case 15
        sample_id = 'Delta_tktAB';
        enzyme_id = 85;
    
    case 16
        sample_id = 'Delta_pykA';
        enzyme_id = 88;

    case 17
        sample_id = 'Delta_fbp';
        enzyme_id = 38;

    case 18
        sample_id = 'Delta_tpi';
        enzyme_id = 87;
        
    case 19
        sample_id = 'Delta_eda';
        enzyme_id = 89;
        
    case 20
        sample_id = 'Delta_edd';
        enzyme_id = 90;
end

fprintf('Working on %s\n',sample_id)

[vnet, conc] = Main_Module(enzyme_id,modification,condition);
Vnet = real(vnet); Conc = real(conc);

save(sprintf('result_cont_%s.mat',sample_id),'Vnet','Conc');
toc;
end

return