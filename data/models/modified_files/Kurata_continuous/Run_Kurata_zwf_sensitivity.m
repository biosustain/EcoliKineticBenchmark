clear all;
close all;
tic;

global strain_no continuous_flg ;

options = odeset('RelTol',1e-10,'AbsTol',1e-10); %
n_FLUX = 155;
ode_file = @KurataModel_ODE;
flux_file = @KurataModel_Flux;

for i = 1:1:3

    switch i     
    case 1
    % GR04 WT, 0.7h-1
    SampleID = 'WT';
    strain_no = 1;
    continuous_flg = 5;
    
    % Delta zwf(G6Pdh)
    case 2
    SampleID = 'dzwf'; %Delta zwf (G6pdh)
    strain_no = 15;
    continuous_flg = 5;

    % zwf(15) overexpression
    case 3
    SampleID = 'zwf(15)';
    strain_no = 30;
    continuous_flg = 5;
    end

%% Simulation
span = -10:0.1:200;
y0 = getInitialCondition();
y0( 1) = 1; % X
y0( 2) = 1; % GLCex
y0(12) = 1; % ACEex
fprintf('Working on %s\n',SampleID);
[ T, Y, FLUX ] = runSimulation(ode_file,flux_file,span,y0,options);
Y = real(Y); FLUX = real(FLUX);
%plotContinuousFlux(FLUX(end,:), SampleID);
save(sprintf('kurata_zwf_sens_%s.mat',SampleID),'T','Y','FLUX');
toc;

end

return