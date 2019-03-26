function plotContinuousFlux(Flux_sim, SampleID)

[~, Fidx] = setIndex();

Flux_exp = ExpDataForContinuousCulture(SampleID);

Flux_sim = Flux_sim / Flux_sim(Fidx.vPts4) * 100; %normalization by vPts4
Flux_exp = Flux_exp / Flux_exp(1) * 100;          %normalization by vPts4

%  1: Glucose + PEP -> G6P + PYR vs. vPts4
Flux_sim_temp(1) = Flux_sim(Fidx.vPts4); Flux_exp_temp(1) = Flux_exp(1);

%  2: G6P <-> F6P vs. vE_Pgi
Flux_sim_temp(2) = Flux_sim(Fidx.vE_Pgi); Flux_exp_temp(2) = Flux_exp(2);

%  3: F6P -> F1,6P vs. vE_Pfk - vE_Fbp
Flux_sim_temp(3) = Flux_sim(Fidx.vE_Pfk)-Flux_sim(Fidx.vE_Fbp); Flux_exp_temp(3) = Flux_exp(3);

%  4: F1,6P -> DHAP + G3P vs. vE_Fba
Flux_sim_temp(4) = Flux_sim(Fidx.vE_Fba); Flux_exp_temp(4) = Flux_exp(4);

%  6: G3P -> 3PG vs. vE_Pgk
Flux_sim_temp(5) = Flux_sim(Fidx.vE_Gapdh); Flux_exp_temp(5) = Flux_exp(6);

%  7: 3PG <-> PEP vs. vE_Eno
Flux_sim_temp(6) = Flux_sim(Fidx.vE_Gapdh); Flux_exp_temp(6) = Flux_exp(7);

%  8: PEP -> PYR vs. vE_Pyk - vE_Pps
Flux_sim_temp(7) = Flux_sim(Fidx.vE_Pyk)-Flux_sim(Fidx.vE_Pps); Flux_exp_temp(7) = Flux_exp(8);

%  9: PYR -> AcCoA + CO2 vs. vE_Pdh
Flux_sim_temp(8) = Flux_sim(Fidx.vE_Pdh); Flux_exp_temp(8) = Flux_exp(9);

% 29: AcCoA -> Acetate vs. vE_Ack - vE_Acs
Flux_sim_temp(9) = Flux_sim(Fidx.vE_Ack) - Flux_sim(Fidx.vE_Acs) ; Flux_exp_temp(9) = Flux_exp(29);

% 10: G6P -> 6PG vs. vE_G6pdh
Flux_sim_temp(10) = Flux_sim(Fidx.vE_G6pdh); Flux_exp_temp(10) = Flux_exp(10);

% 11: 6PG -> Ru5P + CO2 vs. vE_6Pgdh
Flux_sim_temp(11) = Flux_sim(Fidx.vE_6Pgdh); Flux_exp_temp(11) = Flux_exp(11);

% 13: RU5P -> R5P vs. vE_R5pi
Flux_sim_temp(12) = Flux_sim(Fidx.vE_R5pi); Flux_exp_temp(12) = Flux_exp(13);

% 12: RU5P -> X5P vs. vE_Rpe
Flux_sim_temp(13) = Flux_sim(Fidx.vE_Rpe); Flux_exp_temp(13) = Flux_exp(12);

% 14: R5P + X5P <-> S7P + G3P vs. vE_Tkt1
Flux_sim_temp(14) = Flux_sim(Fidx.vE_Tkt1); Flux_exp_temp(14) = Flux_exp(14);

% 16: X5P + E4P <-> F6P + G3P vs. vE_Tkt2
Flux_sim_temp(15) = Flux_sim(Fidx.vE_Tkt2); Flux_exp_temp(15) = Flux_exp(16);

% 15: S7P + G3P <-> E4P + F6P vs. vE_Tal
Flux_sim_temp(16) = Flux_sim(Fidx.vE_Tal); Flux_exp_temp(16) = Flux_exp(15);

% 28: 6-PG -> G3P + PYR vs. vE_Edd 
Flux_sim_temp(17) = Flux_sim(Fidx.vE_Edd); Flux_exp_temp(17) = Flux_exp(28);

% 28: 6-PG -> G3P + PYR vs. vE_Eda 
Flux_sim_temp(18) = Flux_sim(Fidx.vE_Eda); Flux_exp_temp(18) = Flux_exp(28);

% 17: AcCoA + OAA -> CIT vs. vE_Cs
Flux_sim_temp(19) = Flux_sim(Fidx.vE_Cs); Flux_exp_temp(19) = Flux_exp(17);

% 19: ICT -> 2-KG + CO2 vs. vE_Icdh
Flux_sim_temp(20) = Flux_sim(Fidx.vE_Icdh); Flux_exp_temp(20) = Flux_exp(19);

% 20: 2-KG -> SUC + CO2 vs. vE_akgdh
Flux_sim_temp(21) = Flux_sim(Fidx.vE_akgdh); Flux_exp_temp(21) = Flux_exp(20);

% 21: SUC -> FUM vs. vE_Sdh
Flux_sim_temp(22) = Flux_sim(Fidx.vE_Sdh); Flux_exp_temp(22) = Flux_exp(21);

% 22: FUM -> MAL vs. vE_Fum
Flux_sim_temp(23) = Flux_sim(Fidx.vE_Fum); Flux_exp_temp(23) = Flux_exp(22);

% 23: MAL <-> OAA vs. vE_Mdh
Flux_sim_temp(24) = Flux_sim(Fidx.vE_Mdh); Flux_exp_temp(24) = Flux_exp(23);

% 25: MAL -> PYR + CO2 vs. vE_MaeB
Flux_sim_temp(25) = Flux_sim(Fidx.vE_MaeB); Flux_exp_temp(25) = Flux_exp(25);

% 24: PEP + CO2 <-> OAA vs. vE_Ppc - vE_Pck
Flux_sim_temp(26) = Flux_sim(Fidx.vE_Ppc)-Flux_sim(Fidx.vE_Pck); Flux_exp_temp(26) = Flux_exp(24);

% 26: ICT -> Glyoxylate + SUC vs. vE_Icl
Flux_sim_temp(27) = Flux_sim(Fidx.vE_Icl); Flux_exp_temp(27) = Flux_exp(26);

% 27: Glyoxylate + AcCoA -> MAL vs. vE_Ms
Flux_sim_temp(28) = Flux_sim(Fidx.vE_Ms); Flux_exp_temp(28) = Flux_exp(27);

% 40: AcCoA -> (Cell synthesis) vs. vBM_AcCoA
Flux_sim_temp(29) = Flux_sim(Fidx.vBM_AcCoA); Flux_exp_temp(29) = Flux_exp(40);

% 35: E4P -> (Cell synthesis) vs. vBM_E4P
Flux_sim_temp(30) = Flux_sim(Fidx.vBM_E4P); Flux_exp_temp(30) = Flux_exp(35);

% 33: F6P -> (Cell synthesis) vs. vBM_F6P
Flux_sim_temp(31) = Flux_sim(Fidx.vBM_F6P); Flux_exp_temp(31) = Flux_exp(33);

% 36: G3P -> (Cell synthesis) + 37: 3PG -> (Cell synthesis) vs. vBM_GAP
Flux_sim_temp(32) = Flux_sim(Fidx.vBM_GAP); Flux_exp_temp(32) = Flux_exp(36)+Flux_exp(37);

% 32: G6P -> (Cell synthesis) vs. vBM_G6P
Flux_sim_temp(33) = Flux_sim(Fidx.vBM_G6P); Flux_exp_temp(33) = Flux_exp(32);

% 41: OAA -> (Cell synthesis) vs. vBM_OAA
Flux_sim_temp(34) = Flux_sim(Fidx.vBM_OAA); Flux_exp_temp(34) = Flux_exp(41);

% 38: PEP -> (Cell synthesis) vs. vBM_PEP
Flux_sim_temp(35) = Flux_sim(Fidx.vBM_PEP); Flux_exp_temp(35) = Flux_exp(38);

% 39: PYR -> (Cell synthesis) vs. vBM_PYR
Flux_sim_temp(36) = Flux_sim(Fidx.vBM_PYR); Flux_exp_temp(36) = Flux_exp(39);

% 34: R5P -> (Cell synthesis) vs. vBM_R5P
Flux_sim_temp(37) = Flux_sim(Fidx.vBM_R5P); Flux_exp_temp(37) = Flux_exp(34);

% 42: 2KG -> (Cell synthesis) vs. vBM_aKG
Flux_sim_temp(38) = Flux_sim(Fidx.vBM_aKG); Flux_exp_temp(38) = Flux_exp(42);


%%
scrsz = get(0,'ScreenSize');
h = figure('Position',[10 scrsz(4)*0.24 scrsz(3)*0.24 scrsz(4)*0.24]);

hold on;
plot([-1e+3 1e+3],[-1e+3 1e+3],'k-');
scatter(Flux_exp_temp, Flux_sim_temp);

xlabel('Measured Flux (-)','FontSize',10,'FontName','Arial');
ylabel('Simulated Flux (-)','FontSize',10,'FontName','Arial');
set(gca,'FontSize',10,'FontName','Arial');
title(SampleID);
xlim([-10 200])
ylim([-10 200]);
box on;
hold off;

return



