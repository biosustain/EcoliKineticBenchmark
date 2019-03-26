function [Vnet,conc] = solve_ode(expanded_model,t_interval,y0,kinetic_param)

% This function solves the ODEs and finds the rate of elementary and
% overall rxns (PMID:24928774).
%
%
% INPUTS:
% -------
%       expanded_model: model structure that containts the information of the
%       expanded and unexpanded model
%           t_interval: The time interval for simulations
%                   y0: The initial concentration of the metabolites and
%                       enzymes after decomposition of the reactions into
%        kinetic_param: kinetic parameters of rxns
% 
%
% OUTPUTS:
% ------
%           conc: Normalized metabolite concentrations and enzyme fractions:
%                 Columns correspond to the time points 
%                 Rows correspond to metabolites 
%           Vnet: Rate of overall rxns in the unexpanded model at each time
%           point:
%                 Columns correspond to time points
%                 Row correspond to reactions
%
%
% Ali Khodayari, Costas Maranas Lab @ Penn State

S_f_b = expanded_model.S_f_b;

try

  % Solve system of ODEs
  a=cputime;
  eval(strcat('options=odeset(''Jacobian'',@(t,y)jacobFun_',expanded_model.model_name,...
      '(t,y,kinetic_param),''Events'',@(t,y)event_function(t,y,a,kinetic_param,expanded_model));')); 
  %eval(strcat('options=odeset(''Jacobian'',@(t,y)jacobFun_',expanded_model.model_name,'(t,y,kinetic_param));')); 
  
  % you may need to change the solver to ode23s for anaerobic condition
  [T,Y]=ode15s(@(t,y)mass_balance_ode(t,y,kinetic_param,S_f_b),t_interval,y0,options);

  %---------  Concentrations ---------------
  % Rows in conc correspond to metabolites and columns correspond to time points
  conc = Y'; 
  
  %---------- Compute reaction rates ----------------
  % Total number of time points
  time_point_num=length(T);
  
  v = zeros(length(kinetic_param),time_point_num);
  
  % First compute the rate of reactions 
  for j = 1:length(kinetic_param)
   for t=1:time_point_num
      v(j,t) = kinetic_param(j);
  
      for i=1:size(S_f_b,1)
          if S_f_b(i,j) < 0
              v(j,t)=v(j,t)*conc(i,t);
          end
      end
   end
  end
  
  %---------- Compute Vnet ----------------
  Vnet = zeros(size(expanded_model.unexp_rxn_info,1),time_point_num);
  
  for j=1:size(expanded_model.unexp_rxn_info,1)
    for t=1:time_point_num
       % If this is not an export reaction which is irreversible
       if expanded_model.unexp_rxn_info{j,6} ~= 3 
          % index of the first elementary rxn corresponding to j
          elem_index_first = expanded_model.unexp_rxn_info{j,2};
  
          % Corresponding indices of v_f and v_b. The first element is index of the forward
          % rxn whereas the second element is the index of the backward rxn
          v_indices=find(expanded_model.rxn_f_bInd_rxnInd==elem_index_first);
  
          % Vnet = v_f - v_b
          % added to take care of irreversible reactions, if any
          if expanded_model.rev(elem_index_first)==1
            Vnet(j,t) = v(v_indices(1),t) - v(v_indices(2),t);
          else
              if isempty(find(expanded_model.unexpanded_model.export,j))
                  Vnet(j,t) = v(v_indices(1),t);
              else
                  Vnet(j,t) = -v(v_indices(1),t);% for those irreversible rxns that already multiplied by -1
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
conc = Y(:,cell2mat(expanded_model.metabIndexInfo(1,2)):cell2mat(expanded_model.metabIndexInfo(1,3)))';

catch Me
  % If there is any errors in integration
  fprintf('\nError in integration\n');
  
  T=t_interval;
  conc = NaN(length(expanded_model.metab),length(T));
  v= NaN(length(expanded_model.rxn_f_b),length(T));
  Vnet = NaN(length(expanded_model.unexp_rxn_info(:,1)),length(T));

end   % end of try catch
 
clearvars -except Vnet conc
conc=real(conc);
Vnet=real(Vnet);
end
