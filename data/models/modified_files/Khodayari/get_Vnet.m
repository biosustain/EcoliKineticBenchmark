function [Vnet] = get_Vnet(expanded_model,t_interval,conc,kinetic_param)

% This function evaluates conservation of mass (for metabolites, enzymes and enzyme complexes) at each time 
% interval to calculate tolerance error for ODE's integration (PMID:24928774). 
%
% INPUTS:
% -------
%       expanded_model: model structure that containts the information of the
%           t_interval: The time interval for simulations
%                 conc: Vector of metabolite and enzyme concentrations
%        kinetic_param: Vector of elementary kinetic parameters
%
%
% OUTPUTS:
% ------
%                 Vnet: Rate of overall rxns in the unexpanded model at the given
%                       time point t_interval
%
%
% Ali Khodayari, Costas Maranas Lab @ Penn State

conc=conc';

S_f_b = expanded_model.S_f_b;
 
  %---------- Compute reaction rates ----------------
  % Total number of time points
  time_point_num=length(t_interval);
  
  v = zeros(length(kinetic_param),time_point_num);
  
  % First compute the rate of reactions 
  for j = 1:length(kinetic_param)
   for t=time_point_num
      v(j,t) = kinetic_param(j);
  
      for i=1:size(S_f_b,1)
          if S_f_b(i,j) < 0
              v(j,t)=v(j,t)*(conc(i,t)^abs(S_f_b(i,j)));
          end
      end
   end
  end
  
  %---------- Compute Vnet ----------------
  Vnet = zeros(size(expanded_model.unexp_rxn_info,1),time_point_num);
  
  for j=1:size(expanded_model.unexp_rxn_info,1)
    for t=time_point_num
       % If this is not an export reaction which is irreversible
       if expanded_model.unexp_rxn_info{j,6} ~= 3 
          % index of the first elementary rxn corresponding to j
          elem_index_first = expanded_model.unexp_rxn_info{j,2};
  
          % Corresponding indices of v_f and v_b. The first element is index of the forward
          % rxn whereas the second element is the index of the backward rxn
          v_indices=find(expanded_model.rxn_f_bInd_rxnInd==elem_index_first);
  
          % Vnet = v_f - v_b
          if expanded_model.rev(elem_index_first)==1
            Vnet(j,t) = v(v_indices(1),t) - v(v_indices(2),t);
          else
              if isempty(find(expanded_model.unexpanded_model.export==j))
                  Vnet(j,t) = v(v_indices(1),t);
              else
                  Vnet(j,t) = -v(v_indices(1),t);% for those irreversible rxns that already multiplied by -1 (transport rxn)
              end
          end
  
          % Check if v_f - v_b are equal for elemetnary steps
          for k=expanded_model.unexp_rxn_info{j,2}:expanded_model.unexp_rxn_info{j,3}
             v_k_indices=find(expanded_model.rxn_f_bInd_rxnInd==k);
             dev=[];
          if expanded_model.rev(elem_index_first)==1
            dev=abs(v(v_k_indices(1),t) - v(v_k_indices(2),t) - Vnet(j,t));
          else
            dev=abs(v(v_k_indices(1),t) - Vnet(j,t));
          end             
             
             if dev > 0.1
  %                 fprintf('(v_f - v_b) not equal for all elementary steps of reaction %s at time %d\n',expanded_model.unexp_rxn_info{j,1},t);
             end
          end
  
       % if this an export rxn which is irreversilbe
       else
          % Vnet = v
          elem_index = expanded_model.unexp_rxn_info{j,2};
          v_index=find(expanded_model.rxn_f_bInd_rxnInd==elem_index);
          Vnet(j,t) = v(v_index,t);
       end
    end   % end for t
  end
clearvars -except Vnet
Vnet=real(Vnet);
end
