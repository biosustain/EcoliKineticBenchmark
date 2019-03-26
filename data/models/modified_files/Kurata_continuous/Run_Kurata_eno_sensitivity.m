clear all;
close all;
tic;

global strain_no continuous_flg ;

options = odeset('RelTol',1e-10,'AbsTol',1e-10); %
n_FLUX = 155;
ode_file = @KurataModel_ODE;
flux_file = @KurataModel_Flux;

for i = 1:1:5
    switch i     
    case 1
    % GR04 WT, 0.7h-1
    SampleID = 'WT';
    strain_no = 1;
    continuous_flg = 5;
    
    % eno(0)
    case 2
    SampleID = 'eno(0)';
    strain_no = 36;
    continuous_flg = 5;

    % eno(50)
    case 3
    SampleID = 'eno(50)';
    strain_no = 37;
    continuous_flg = 5;

    % eno(200)
    case 4
    SampleID = 'eno(200)';
    strain_no = 38;
    continuous_flg = 5;

    % eno(500)
    case 5
    SampleID = 'eno(500)';
    strain_no = 39;
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
    save(sprintf('kurata_eno_sens_%s.mat',SampleID),'T','Y','FLUX');
    toc;
end

return