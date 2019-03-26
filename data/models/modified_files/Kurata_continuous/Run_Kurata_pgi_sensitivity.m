clear all;
close all;
tic;

global strain_no continuous_flg ;

options = odeset('RelTol',1e-10,'AbsTol',1e-10); %
n_FLUX = 155;
ode_file = @KurataModel_ODE;
flux_file = @KurataModel_Flux;

for i = 1:1:6
    switch i     
    case 1
    % GR04 WT, 0.7h-1
    SampleID = 'WT';
    strain_no = 1;
    continuous_flg = 5;
        
    % Delta pgi
    case 2
    SampleID = 'dpgi'; %Delta zwf (G6pdh)
    strain_no = 4;
    continuous_flg = 5;
    
    % pgi(0)
    case 3
    SampleID = 'pgi(0)';
    strain_no = 31;
    continuous_flg = 5;
    
    % pgi(20)
    case 4
    SampleID = 'pgi(20)';
    strain_no = 32;
    continuous_flg = 5;

    % pgi(50)
    case 5
    SampleID = 'pgi(50)';
    strain_no = 33;
    continuous_flg = 5;

    % pgi(100)
    case 6
    SampleID = 'pgi(100)';
    strain_no = 34;
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
save(sprintf('kurata_pgi_sens_%s.mat',SampleID),'T','Y','FLUX');
toc;

end
return