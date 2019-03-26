function [value,isterminal,direction] = event_function(t_interval,y,CPU_t,kinetic_param,expanded_model)

% This function evaluates conservation of mass (for metabolites, enzymes and enzyme complexes) at each time interval and terminates ODE's
% integration once the model reaches to a steady-state solution (S.V~tol). 
% The function allows the integration to be continineud on CPU for 5000 seconds by which it will be stopped (PMID:24928774).
%
% INPUTS:
% -------
%     t_interval: The time interval for simulations
%              y: Vector of metabolite and enzyme concentrations
%          CPU_t: CPU time
%  kinetic_param: Vector of elementary kinetic parameters
%          S_f_b: Stoichiometric matrix for the elementary rxns
%
%
% OUTPUTS:
% ------
%          value: Value of the event function
%     isterminal: 1 if the integration is to terminate at a zero of this event function, otherwise, 0
%      direction: 0 if all zeros are to be located (the default), +1 if only zeros where the event 
%                 function is increasing, and -1 if only zeros where the event function is decreasing
%
%
% Ali Khodayari, Costas Maranas Lab @ Penn State

    % Elapsed CPU time
    c=cputime-CPU_t;
    
    direction = 0;
    Vnet_perturb_ss=get_Vnet(expanded_model,t_interval,y',kinetic_param);
    Error=max(abs(expanded_model.unexpanded_model.S*Vnet_perturb_ss(1:size(expanded_model.unexpanded_model.S,2))));
    tol=5;
    
    if (c<1000 && Error>tol)
        value=Error;
        isterminal = 0;
    else
        value=0;
        isterminal = 1;
    end
    
end