%function Run_KurataModel_Continuous()
clear all;
close all;
tic;

global strain_no continuous_flg ;

options = odeset('RelTol',1e-10,'AbsTol',1e-10); %
n_FLUX = 155;
ode_file = @KurataModel_ODE;
flux_file = @KurataModel_Flux;

for i= 1:1:29
switch i
    case 1
    % RF06, WT(Oct),WT, 0.2h-1
    SampleID = 'RF06';
    strain_no = 1;
    continuous_flg = 1;
    case 2
    % GR02 WT, 0.4h-1
    SampleID = 'GR02';
    strain_no = 1;
    continuous_flg = 3;        
    case 3
    % GR03 WT, 0.5h-1
    SampleID = 'GR03';
    strain_no = 1;
    continuous_flg = 4;       
    case 4
    % GR04 WT, 0.7h-1
    SampleID = 'GR04';
    strain_no = 1;
    continuous_flg = 5;
    
    case 5
    SampleID = 'Delta_glk'; %Delta glk
    strain_no = 2;
    continuous_flg = 1;

    case 6
    SampleID = 'Delta_pgi'; %Delta pgi
    strain_no = 4;
    continuous_flg = 1;

    case 7
    SampleID = 'Delta_pfkA'; %Delta pfkA
    strain_no = 5;
    continuous_flg = 1;
    
    case 8
    SampleID = 'Delta_pfkB'; %Delta pfkB
    strain_no = 6;
    continuous_flg = 1;

    case 9
    SampleID = 'Delta_fbp'; %Delta fbp
    strain_no = 7;
    continuous_flg = 1;

    case 10
    SampleID = 'Delta_fbaB'; %Delta fbaB
    strain_no = 8;
    continuous_flg = 1;

    case 11
    SampleID = 'Delta_gpmA'; % Delta gpmA
    strain_no = 10;
    continuous_flg = 1;
    
    case 12
    SampleID = 'Delta_pykA'; % Delta pykA
    strain_no = 12;
    continuous_flg = 1;
    
    case 13
    SampleID = 'Delta_pykF'; % Delta pykF
    strain_no = 13;
    continuous_flg = 1;
    
    case 14
    SampleID = 'Delta_ppsA'; %Delta ppsA
    strain_no = 14;
    continuous_flg = 1;

    case 15
    SampleID = 'Delta_zwf'; %Delta zwf (G6pdh)
    strain_no = 15;
    continuous_flg = 1;

    case 16
    SampleID =  'Delta_pgl'; %Delta pgl
    strain_no = 16;
    continuous_flg = 1;

    case 17
    SampleID = 'Delta_gnd'; % Delta gnd (6Pgdh)
    strain_no = 17;
    continuous_flg = 1;

    case 18
    SampleID = 'Delta_rpe'; %Delta rpe 
    strain_no = 18;
    continuous_flg = 1;

    case 19
    SampleID = 'Delta_rpiA'; % Delta rpiA (R5pi)
    strain_no = 19;
    continuous_flg = 1;

    case 20
    SampleID = 'Delta_rpiB'; % Delta rpiB (R5pi)
    strain_no = 20;
    continuous_flg = 1;
    
    case 21
    SampleID = 'Delta_tktA'; % Delta tktA
    strain_no = 21;
    continuous_flg = 1;

    case 22
    SampleID = 'Delta_tktB'; %Delta tktB
    strain_no = 22;
    continuous_flg = 1;

    case 23
    SampleID = 'Delta_talA'; %Delta talA 
    strain_no = 23;
    continuous_flg = 1;
    
    case 24
    SampleID = 'Delta_talB'; %Delta talB 
    strain_no = 24;
    continuous_flg = 1;
    
    case 25
    % Delta ppc, D=0.2h-1
    SampleID = 'Delta_ppc'; %Delta ppc
    strain_no = 25;
    continuous_flg = 1;

    %% Denis additions
    case 26
    % Delta AKGDH
    SampleID = 'Delta_sucAC';
    strain_no = 26;
    continuous_flg = 1;
    
    case 27
    % Delta Pts
    SampleID = 'Delta_pts';
    strain_no = 27;
    continuous_flg = 1;
    
    case 28
    % Delta SdhC
    SampleID = 'Delta_sdhC';
    strain_no = 28;
    continuous_flg = 1;

    case 29
    % Delta tpi
    SampleID = 'Delta_tpi';
    strain_no = 29;
    continuous_flg = 1;
end

%% Simulation
span = -10:0.1:200;
y0 = getInitialCondition();
y0( 1) = 1; % X
y0( 2) = 1; % GLCex
y0(12) = 1; % ACEex
[ T, Y, FLUX ] = runSimulation(ode_file,flux_file,span,y0,options);
Y = real(Y); FLUX = real(FLUX);
%plotContinuousFlux(FLUX(end,:), SampleID);
save(sprintf('result_cont_%s.mat',SampleID),'T','Y','FLUX');
toc;

end

return

