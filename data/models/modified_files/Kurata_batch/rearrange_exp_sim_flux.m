function [Flux_exp, Flux_sim] = rearrange_exp_sim_flux(T, FLUX, sampling_time)

global strain_no

n_flux=42;
i=1; % time index

[T_exp, F_exp] = ExpDataForBatchCulture_Flux();

    switch strain_no
        case 1
            j=1; %5h WT
        case 26
            j=1; %5h DpykA/pykF
        case 4
            j=1; %16h Dpgi
        case 25
            j=1; %Dppc  
    end
    % rearrangement of experimental fluxes 
    Flux_exp(i,1)  = F_exp(j,1); %vPts4
    Flux_exp(i,2)  = F_exp(j,2); %vE_Pgi
    Flux_exp(i,3)  = F_exp(j,3); %vE_Pfk-vE_Fbp
    Flux_exp(i,4)  = F_exp(j,4); %vE_Fba
    Flux_exp(i,5)  = F_exp(j,6); %vE_Pgk (Gapdh)
    Flux_exp(i,6)  = F_exp(j,7); %vE_Eno (Gapdh)
    Flux_exp(i,7)  = F_exp(j,8); %vE_Pyk-vE_Pps
    Flux_exp(i,8)  = F_exp(j,9); %vE_Pdh
    Flux_exp(i,9)  = F_exp(j,10); %vE_Ack
    Flux_exp(i,10) = F_exp(j,11); %vE_Pgk
    Flux_exp(i,11) = F_exp(j,12); %vE_Eno
    Flux_exp(i,12) = F_exp(j,14); %vE_R5pi
    Flux_exp(i,13) = F_exp(j,13); %vE_Rpe 
    Flux_exp(i,14) = F_exp(j,15); %vE_Tkt1
    Flux_exp(i,15) = F_exp(j,17); %vE_Tkt2
    Flux_exp(i,16) = F_exp(j,16); %vE_Tal
    Flux_exp(i,17) = F_exp(j,28); %vE_Edd
    Flux_exp(i,18) = F_exp(j,28); %vE_Eda
    
    Flux_exp(i,19) = F_exp(j,18); %vE_Cs
    Flux_exp(i,20) = F_exp(j,19); %vE_Icdh
    Flux_exp(i,21) = F_exp(j,20); %vE_akgdh
    Flux_exp(i,22) = F_exp(j,21); %vE_Sdh
    Flux_exp(i,23) = F_exp(j,22); %vE_Fum
    Flux_exp(i,24) = F_exp(j,23); %vE_Mdh
    Flux_exp(i,25) = F_exp(j,25); %vE_MaeB
    Flux_exp(i,26) = F_exp(j,24); %vE_Ppc-vE_Pck
    Flux_exp(i,27) = F_exp(j,26); %vE_Icl
    Flux_exp(i,28) = F_exp(j,27); %vE_Ms
    Flux_exp(i,29:1:n_flux) = 0;
    %Normalization
    Flux_exp(i,:) = Flux_exp(i,:)./F_exp(j,1).*100; % Total glucose uptake = 100

    % rearrangement of simulated fluxes
    i=1; k=101+10*sampling_time;
    Flux_sim(i,1) = FLUX(k,3); %vPts4
    Flux_sim(i,2) = FLUX(k,8); %vE_Pgi
    Flux_sim(i,3) = FLUX(k,9)-FLUX(k, 10); %vE_Pfk-vE_Fbp
    Flux_sim(i,4) = FLUX(k,11); %vE_Fba
    Flux_sim(i,5) = FLUX(k,12); %vE_Gapdh  Pgk
    Flux_sim(i,6) = FLUX(k,12); %vE_Gapdh  Eno
    Flux_sim(i,7) = FLUX(k,13)-FLUX(k,14); %vE_Pyk-vE_Pps
    Flux_sim(i,8) = FLUX(k,15); %vE_Pdh
    Flux_sim(i,9) = FLUX(k,17); %vE_Ack
    Flux_sim(i,10) = FLUX(k,34); %vE_G6pdh
    Flux_sim(i,11) = FLUX(k,38); %vE_6Pgdh
    Flux_sim(i,12) = FLUX(k,39); %vE_R5pi
    Flux_sim(i,13) = FLUX(k,40); %vE_Rpe 
    Flux_sim(i,14) = FLUX(k,41); %vE_Tkt1
    Flux_sim(i,15) = FLUX(k,42); %vE_Tkt2
    Flux_sim(i,16) = FLUX(k,43); %vE_Tal
    Flux_sim(i,17) = FLUX(k,36); %vE_Edd
    Flux_sim(i,18) = FLUX(k,37); %vE_Eda
    
    Flux_sim(i,19) = FLUX(k,21); %vE_Cs
    Flux_sim(i,20) = FLUX(k,22); %vE_Icdh
    Flux_sim(i,21) = FLUX(k,23); %vE_akgdh
    Flux_sim(i,22) = FLUX(k,24); %vE_Sdh
    Flux_sim(i,23) = FLUX(k,25); %vE_Fum
    Flux_sim(i,24) = FLUX(k,26); %vE_Mdh
    Flux_sim(i,25) = FLUX(k,27); %vE_MaeB
    Flux_sim(i,26) = FLUX(k,29)-FLUX(k,28); %vE_Ppc-vE_Pck
    Flux_sim(i,27) = FLUX(k,30); %vE_Icl
    Flux_sim(i,28) = FLUX(k,31); %vE_Ms

    Flux_sim(i,29)  = FLUX(k,123);%AcCoA
    Flux_sim(i,30)  = FLUX(k,129);%E4P
    Flux_sim(i,31)  = FLUX(k,119);%F6P
    Flux_sim(i,32)  = FLUX(k,120);%GAP
    Flux_sim(i,33)  = FLUX(k,118);%G6P
    Flux_sim(i,34)  = FLUX(k,127);%OAA
    Flux_sim(i,35)  = FLUX(k,121);%PEP
    Flux_sim(i,36)  = FLUX(k,122);%PYR
    Flux_sim(i,37)  = FLUX(k,128);%R5P
    Flux_sim(i,38)  = FLUX(k,124);%aKG

    Flux_sim(i,39)  = FLUX(k,13);%Pyk
    Flux_sim(i,40)  = FLUX(k,14);%Pps
    Flux_sim(i,41)  = FLUX(k,29);%Ppc
    Flux_sim(i,42)  = FLUX(k,28);%Pck

    % Normalization
    Flux_sim(i,:) = Flux_sim(i,:)./(FLUX(k,3) + FLUX(k,5)).*100; % Total glucose uptake = 100

return






