function dy = mass_balance_ode(t_interval,y,kinetic_param,S_f_b)
%
% This function computes the conservation of mass for all metabolites and
% enzymes in the form of free and complex (PMID:24928774).
%
%
% INPUTS:
% -------
%    t_interval: The time interval for simulations
%             y: Vector of initial concentrations of metabolites and
%             enzymes
% kinetic_param: Vector of elementary kinetic parameters
%         S_f_b: Stoichiometric matrix of decomposed reactions
%
%
% OUTPUTS:
% ---------
%            dy: Value of dy/dt
%
%
% Ali Khodayari, Costas Maranas Lab @ Penn State

dy = zeros(length(y),1);

% Reaction rates
v = zeros(length(kinetic_param),1);


% Compute the rate of reactions 
for j = 1:length(kinetic_param)

    v(j) = kinetic_param(j);

    for i=1:length(y) 
        if S_f_b(i,j) < 0
            v(j)=v(j)*(y(i)^abs(S_f_b(i,j))); 
        end
    end
end


% Loop over each metabolite
for i = 1:length(y)
    dy(i) = S_f_b(i,:)*v;
end
dy;
