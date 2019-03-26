function [x1,x2]=parallel_automated(id,level,aerobic_enzyme_fraction)

% This function assigns the estimated enzyme level changes for the training mutant condition (PMID:24928774).
%
% INPUTS:
% -------
%                        m: The id of the rection(s) that gets perturbed
%                    level: The level of the perturbation
%  aerobic_enzyme_fraction: the estimated enzyme levels under aerobic
%                           condition
%
% OUTPUTS:
% ------
%                        m: The id of the rection(s) that gets perturbed
%                           upon applying the enzyme level changes 
%                    level: The level of the perturbation upon applying the enzyme level changes
%
%
% Ali Khodayari, Costas Maranas Lab @ Penn State

    if length(id)==1 && length(level)==1
    
        if level<1e-7

                switch id
                    case 51
                        x1=[51 446:455];
                        x2=[1e-8,zeros(1,10)];
                        %Mutant='gnd';

                    case 67
                        x1=[67 446:455];
                        x2=[1e-8,zeros(1,10)];
                        %Mutant='pgi';

                    case 74
                        x1=[74 446:455];
                        x2=[1e-8,zeros(1,10)];
                        %Mutant='ppsA';

                    case 88
                        x1=[88 446:455];
                        x2=[1e-8,zeros(1,10)];
                        %Mutant='pykA';

                    case 76
                        x1=[76 446:455];
                        x2=[1e-8,zeros(1,10)];
                        %Mutant='pykF';

                    case 78
                        x1=[78 446:455];
                        x2=[1e-8,zeros(1,10)];
                        %Mutant='rpe';

                    case 42
                        x1=[42 446:455];
                        x2=[1e-8,zeros(1,10)];
                        %Mutant='zwf';  

                    case 66
                        x1=[66 446:455];
                        x2=[aerobic_enzyme_fraction(8),zeros(1,10)];
                        %Mutant='pfk';

                    case 37
                        x1=[37 446:455];
                        x2=[aerobic_enzyme_fraction(10),zeros(1,10)];
                        %Mutant='fbaB';  

                    case 70
                        x1=[70 446:455];
                        x2=[aerobic_enzyme_fraction(11),zeros(1,10)];
                        %Mutant='gpm';  

                    case 69
                        x1=[69 446:455];
                        x2=[aerobic_enzyme_fraction(13),zeros(1,10)];
                        %Mutant='pgl';  

                    case 79
                        x1=[79 446:455];
                        x2=[aerobic_enzyme_fraction(14),zeros(1,10)];
                        %Mutant='rpi';  

                    case 84
                        x1=[84 446:455];
                        x2=[aerobic_enzyme_fraction(16),zeros(1,10)];
                        %Mutant='tal';  

                    case 85
                        x1=[85 86 446:455];
                        x2=[aerobic_enzyme_fraction(18),aerobic_enzyme_fraction(19),zeros(1,10)];
                        %Mutant='tktA';  
                    case 86
                        x1=[85 86 446:455];
                        x2=[aerobic_enzyme_fraction(20),aerobic_enzyme_fraction(21),zeros(1,10)];
                        %Mutant='tktB';  

                    otherwise
                        x1=[id(:)' 446:455];
                        x2=[level(:)' zeros(1,10)];

                end

        else
            
            x1=[id(:)' 446:455];
            x2=[level(:)' zeros(1,10)];
            
        end
        
    else
        x1=[id(:)' 446:455];
        x2=[level(:); zeros(1,10)];

    end
    
end

