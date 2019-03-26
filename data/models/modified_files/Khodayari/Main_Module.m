%------------------------------------------------------------
% This code integrates the ODEs of mass conservation for metabolites and
% enzymes (free and complexes), following decomposition of reactions into
% their elementary stpes (PMID:24928774).
%
%
% INPUTS:
% -------
%          x1: The id of the rection(s) that gets perturbed
%          x2: The level of the perturbation as follows:
%                       0: knock out
%             less than 1: down regulation (e.g., 0.5 means 2-fold
%             downregulation)
%          greater than 1: up regulation (e.g., 2 means 2-fold overexpression)
%        growth_condition: This represents the carbone substrate and growth
%        condition as follows:
%                       aerobic_glucose
%                       anaerobic_glucose
%                       aerobic_pyruvate
%                       anaerobic_acetate
%               
% 
% OUTPUTS:
% ------
%           Vnet_perturb: Rate of the reactions at different time intervals
%  concentration_perturb: Metabolite concentrations at different time intervals
%
%
% Ali Khodayari, Costas Maranas Lab @ Penn State
%------------------------------------------------------------

function [Vnet_perurb,concentration_perturb]=Main_Module(x1,x2,growth_condition)

% Loading model structure and parameters
load Data.mat
    
% Time interval for ODEs integration
t_interval = [0:1e2];% Change it accordingly

warning off

if nargin<3;
    
   fprintf('\n-------------------------------\n');
   fprintf('\nThree inputs are required to run this module, in form of:\n');
   fprintf('Main_Module(reaction ID, Perturbation level,growth_condition)\n');
   fprintf('Default: reaction ID=1, Perturbation level=1, growth_condition=aerobic_glucose\n')
   t_interval = [0:1e1];% 
   
   x1=1;
   x2=1;
   growth_condition='aerobic_glucose';
   
end

switch growth_condition
    
    case 'aerobic_glucose'
        
        [x1,x2]=parallel_automated(x1,x2,aerobic_enzyme_fraction);
    
    case 'anaerobic_glucose'
        
        x1=[x1 anaerobic_id' 33 448 449 450 452 453 454 455];
        x2=[x2 anaerobic_enzyme_fraction' 1e-8 zeros(1,8)];
                
    case 'aerobic_pyruvate'
        
        x1=[x1 25 44 446 447 448 449 450 451 453 454 455 pyr_id'];
        x2=[x2 1e-8 1e-8 zeros(1,9) pyr_enzyme_fraction'];
                
    case 'aerobic_acetate'
        
        x1=[x1 25 44 446 447 448 449 450 451 452 454 455 acetate_id'];
        x2=[x2 1e-8 1e-8 zeros(1,9) acetate_enzyme_fraction'];   
        
    otherwise
        
        disp('The growth_condition is not valid! ')
        
end

fprintf('\n-------------------------------\n');

       % Initial metabolite and enzyme ocncentrations
       y0 = perturbModel(network_data,metabolite_and_enzyme_conc,x1,x2);
           
       % Integration of ODEs
       [Vnet_perurb,concentration_perturb] = solve_ode(network_data,t_interval,y0,elementary_kinetic_parameter);
           
end
