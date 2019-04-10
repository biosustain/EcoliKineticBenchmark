clear all;
close all;
tic;

global strain_no continuous_flg ;

options = odeset('RelTol',1e-10,'AbsTol',1e-10); %
n_FLUX = 155;
ode_file = @KurataModel_ODE;
flux_file = @KurataModel_Flux;

% Dilutions are 0.2, 0.4, 0.6, 0.7

for i = 1:1:4
    switch i     
    case 1
    % 0.2h-1
    SampleID = 'dilution_02';
    strain_no = 1;
    continuous_flg = 1;
    
    % 0.4h-1
    case 2
    SampleID = 'dilution_04';
    strain_no = 1;
    continuous_flg = 3;

    % 0.6h-1
    case 3
    SampleID = 'dilution_06';
    strain_no = 1;
    continuous_flg = 7;

    % 0.7h-1
    case 4
    SampleID = 'dilution_07';
    strain_no = 1;
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
    save(sprintf('kurata_%s.mat',SampleID),'T','Y','FLUX');
    toc;
end

return