function [y0] = perturbModel(expanded_model,initial_conc,x1,x2)

% This code perturbs the enzyme fraction(s) in the model (PMID:24928774).
%
%
% INPUTS:
% -------
%       expanded_model: model structure that containts the information of the
%       expanded and unexpanded model
%         initial_conc: The initial concentration of the metabolites and
%        enzymes after decomposition of the reactions into
%        their elementary steps
%                   x1: The id of the rection(s) that gets perturbed
%                   x2: The level of the perturbation as follows:
%                              0: knock out
%                    less than 1: down regulation
%                 greater than 1: up regulation (e.g., 2 means 2-fold overexpression)
%
%
% OUTPUTS:
% --------
%                  y0: The vector of metabolite and enzyme concetrations
%                  after imposing the effect of perturbation(s)
%
%
% Ali Khodayari, Costas Maranas Lab @ Penn State

% Set the initial concentrations
y0=initial_conc;
y0(1:length(expanded_model.unexpanded_model.metab))=1;
% Index of the perturbed enzymes
perturb_enz_indices = x1; 
for i=1:length(perturb_enz_indices)
    % Find the index of the free enzyme
    free_enz_index = expanded_model.metabIndexInfo{1,3} + perturb_enz_indices(i);

    % Find the indices of corresponding enzyme complex
    enzComplex_indices = find(expanded_model.enz_enzComplex(:,perturb_enz_indices(i)) ~= 0);

    %--- knockout ---
    if x2(i) == 0

       % Initial condition for this enzyme and its corresponding enzyme complexes are zero
       y0(free_enz_index) = 0;
       y0(enzComplex_indices) = 0;              

    %--- overexpression ---
    elseif x2(i) > 1

       % Sum of the enzyme fractions should be equal to the value of perturbation which
       % is greater than. This sum is currently one. So, if add the difference of the value
       % of perturbaiton and one to the initial value of the enzyme fraction for free enzyme
       % the sum will be equal to the value of perturbation 
       y0(free_enz_index) = y0(free_enz_index) + [x2(i) - 1];

    % if downregulaiton
    elseif x2(i) < 1
       % Sample e between zero and the perturbation value 
       % Total number of enzyme fractions associated with each enzyme
       enz_frac_num = 1 + length(enzComplex_indices); 


       % Draw one sample from the sample pool randomly
       isEqualToPerturbValue = 0;
       while ~isEqualToPerturbValue

       % Perform the sampling
       e_sample_pool=rand(1,enz_frac_num);

       % Normalize to the value of perturbation
       e_sample_pool = x2(i)*e_sample_pool./repmat(sum(e_sample_pool,2),1,enz_frac_num);
       
           if sum(e_sample_pool(1,:)) == x2(i)
               isEqualToPerturbValue = 1;
           end
       end  % end of while
       y0(free_enz_index) = e_sample_pool(1,1);
       y0(enzComplex_indices) = e_sample_pool(1,2:end);

    end % end of if expanded_model.Perturbations(perturb_enz_indices(ni),np) == 0 else ...
end   

clearvars -except y0
end
