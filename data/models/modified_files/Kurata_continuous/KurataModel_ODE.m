function dydt = KurataModel_ODE(t,y)

for i=1:1:51
    if y(i)<0
        y(i)= 1e-12;
    end
end

Flux = KurataModel_Flux(t,y);

vgrowth        = Flux(1);
vPts1          = Flux(2);
vPts4          = Flux(3);
vPts4_medium   = Flux(4);
vNonpts        = Flux(5);
vNonpts_medium = Flux(6);
vE_Glk         = Flux(7);
vE_Pgi         = Flux(8);
vE_Pfk         = Flux(9);
vE_Fbp         = Flux(10);
vE_Fba         = Flux(11);
vE_Gapdh       = Flux(12);
vE_Pyk         = Flux(13);
vE_Pps         = Flux(14);
vE_Pdh         = Flux(15);
vE_Pta         = Flux(16);
vE_Ack         = Flux(17);
vE_Ack_medium  = Flux(18);
vE_Acs         = Flux(19);
vE_Acs_medium  = Flux(20);
vE_Cs          = Flux(21);
vE_Icdh        = Flux(22);
vE_akgdh       = Flux(23);
vE_Sdh         = Flux(24);
vE_Fum         = Flux(25);
vE_Mdh         = Flux(26);
vE_MaeB         = Flux(27);
vE_Pck         = Flux(28);
vE_Ppc         = Flux(29);
vE_Icl         = Flux(30);
vE_Ms          = Flux(31);
vE_AceKki      = Flux(32);
vE_AceKph      = Flux(33);
vE_G6pdh       = Flux(34);
vE_Pgl         = Flux(35);
vE_Edd         = Flux(36);
vE_Eda         = Flux(37);
vE_6Pgdh       = Flux(38);
vE_R5pi        = Flux(39);
vE_Rpe        = Flux(40);
vE_Tkt1        = Flux(41);
vE_Tkt2        = Flux(42);
vE_Tal         = Flux(43);
vE_Cya         = Flux(44);
vE_cAMPdegr    = Flux(45);
vG_glk         = Flux(46);
vG_pfkA        = Flux(47);
vG_fbp         = Flux(48);
vG_fbaA        = Flux(49);
vG_gapA        = Flux(50);
vG_pykF        = Flux(51);
vG_ppsA        = Flux(52);
vG_lpd         = Flux(53);
vG_acs         = Flux(54);
vG_gltA        = Flux(55);
vG_icdA        = Flux(56);
vG_sucAB       = Flux(57);
vG_sdhCDAB     = Flux(58);
vG_fumABC      = Flux(59);
vG_mdh         = Flux(60);
vG_maeB        = Flux(61);
vG_pckA        = Flux(62);
vG_ppc         = Flux(63);
vG_aceA        = Flux(64);
vG_aceB        = Flux(65);
vG_aceK        = Flux(66);
vD_X           = Flux(67);
vD_GLCfeed     = Flux(68);
vD_GLCex       = Flux(69);
vD_GLC         = Flux(70);
vD_G6P         = Flux(71);
vD_F6P         = Flux(72);
vD_FBP         = Flux(73);
vD_GAP         = Flux(74);
vD_PEP         = Flux(75);
vD_PYR         = Flux(76);
vD_AcCoA       = Flux(77);
vD_AcP         = Flux(78);
vD_ACEex       = Flux(79);
vD_ICIT        = Flux(80);
vD_aKG         = Flux(81);
vD_SUC         = Flux(82);
vD_FUM         = Flux(83);
vD_MAL         = Flux(84);
vD_OAA         = Flux(85);
vD_GOX         = Flux(86);
vD_6PGL        = Flux(87);
vD_6PG         = Flux(88);
vD_KDPG        = Flux(89);
vD_RU5P        = Flux(90);
vD_R5P         = Flux(91);
vD_X5P         = Flux(92);
vD_S7P         = Flux(93);
vD_E4P         = Flux(94);
vD_cAMP        = Flux(95);
vD_Glk         = Flux(96);
vD_Pfk         = Flux(97);
vD_Fbp         = Flux(98);
vD_Fba         = Flux(99);
vD_Gapdh       = Flux(100);
vD_Pyk         = Flux(101);
vD_Pps         = Flux(102);
vD_Pdh         = Flux(103);
vD_Acs         = Flux(104);
vD_Cs          = Flux(105);
vD_Icdh        = Flux(106);
vD_IcdhP       = Flux(107);
vD_akgdh       = Flux(108);
vD_Sdh         = Flux(109);
vD_Fum         = Flux(110);
vD_Mdh         = Flux(111);
vD_MaeB        = Flux(112);
vD_Pck         = Flux(113);
vD_Ppc         = Flux(114);
vD_Icl         = Flux(115);
vD_Ms          = Flux(116);
vD_AceK        = Flux(117);
vBM_G6P        = Flux(118);
vBM_F6P        = Flux(119);
vBM_GAP        = Flux(120);
vBM_PEP        = Flux(121);
vBM_PYR        = Flux(122);
vBM_AcCoA      = Flux(123);
vBM_aKG        = Flux(124);
vBM_SUC        = Flux(125);
vBM_FUM        = Flux(126);
vBM_OAA        = Flux(127);
vBM_R5P        = Flux(128);
vBM_E4P        = Flux(129);

%********** MODEL STATES
X_dot      = vgrowth - vD_X;
GLCex_dot  = vD_GLCfeed - vD_GLCex - vPts4_medium - vNonpts_medium;
GLC_dot    = vNonpts - vE_Glk - vD_GLC;
G6P_dot    = vPts4 + vE_Glk - vE_Pgi - vE_G6pdh - vD_G6P - vBM_G6P;
F6P_dot    = vE_Pgi + vE_Tkt2 + vE_Tal + vE_Fbp - vE_Pfk - vD_F6P - vBM_F6P;
FBP_dot    = vE_Pfk - vE_Fba - vE_Fbp - vD_FBP;
GAP_dot    = 2*vE_Fba + vE_Tkt1 + vE_Tkt2 + vE_Eda - vE_Tal - vE_Gapdh - vD_GAP - vBM_GAP;
PEP_dot    = vE_Gapdh + vE_Pck + vE_Pps - vE_Pyk - vPts1 - vE_Ppc - vD_PEP - vBM_PEP;
PYR_dot    = vE_Pyk + vPts1 + vE_MaeB + vE_Eda - vE_Pdh - vE_Pps - vD_PYR - vBM_PYR;
AcCoA_dot  = vE_Pdh + vE_Acs - vE_Cs - vE_Ms - vE_Pta - vD_AcCoA - vBM_AcCoA;
AcP_dot    = vE_Pta - vE_Ack - vD_AcP;
ACEex_dot  = vE_Ack_medium - vE_Acs_medium - vD_ACEex;
ICIT_dot   = vE_Cs - vE_Icdh - vE_Icl - vD_ICIT;
aKG_dot    = vE_Icdh - vE_akgdh - vD_aKG - vBM_aKG;
SUC_dot    = vE_akgdh + vE_Icl - vE_Sdh - vD_SUC - vBM_SUC;
FUM_dot    = vE_Sdh - vE_Fum - vD_FUM - vBM_FUM;
MAL_dot    = vE_Fum + vE_Ms - vE_Mdh - vE_MaeB - vD_MAL;
OAA_dot    = vE_Mdh + vE_Ppc - vE_Pck - vE_Cs - vD_OAA - vBM_OAA;
GOX_dot    = vE_Icl - vE_Ms - vD_GOX;
sixPGL_dot = vE_G6pdh - vE_Pgl - vD_6PGL;
sixPG_dot  = vE_Pgl - vE_6Pgdh - vE_Edd - vD_6PG;
KDPG_dot   = vE_Edd - vE_Eda - vD_KDPG;
RU5P_dot   = vE_6Pgdh - vE_Rpe - vE_R5pi - vD_RU5P;
R5P_dot    = vE_R5pi - vE_Tkt1 - vD_R5P - vBM_R5P;
X5P_dot    = vE_Rpe - vE_Tkt1 - vE_Tkt2 - vD_X5P;
S7P_dot    = vE_Tkt1 - vE_Tal - vD_S7P;
E4P_dot    = vE_Tal - vE_Tkt2 - vD_E4P - vBM_E4P;
cAMP_dot   = vE_Cya - vE_cAMPdegr - vD_cAMP;
EIIAP_dot  = vPts1 - vPts4;
Glk_dot    = vG_glk - vD_Glk;
Pfk_dot    = vG_pfkA - vD_Pfk;
Fbp_dot    = vG_fbp - vD_Fbp;
Fba_dot    = vG_fbaA - vD_Fba;
GAPDH_dot  = vG_gapA - vD_Gapdh;
Pyk_dot    = vG_pykF - vD_Pyk;
Pps_dot    = vG_ppsA - vD_Pps;
PDH_dot    = vG_lpd - vD_Pdh;
Acs_dot    = vG_acs - vD_Acs;
CS_dot     = vG_gltA - vD_Cs;
ICDH_dot   = vG_icdA + vE_AceKph - vE_AceKki - vD_Icdh;
ICDHP_dot  = vE_AceKki - vE_AceKph - vD_IcdhP;
aKGDH_dot  = vG_sucAB - vD_akgdh;
SDH_dot    = vG_sdhCDAB - vD_Sdh;
Fum_dot    = vG_fumABC - vD_Fum;
MDH_dot    = vG_mdh - vD_Mdh;
MaeB_dot    = vG_maeB - vD_MaeB;
Pck_dot    = vG_pckA - vD_Pck;
Ppc_dot    = vG_ppc - vD_Ppc;
Icl_dot    = vG_aceA - vD_Icl;
MS_dot     = vG_aceB - vD_Ms;
AceK_dot   = vG_aceK - vD_AceK;

if t < 0
    X_dot     = 0;
    GLCex_dot = 0;
    ACEex_dot = 0;
end

dydt = [ X_dot; GLCex_dot; GLC_dot; G6P_dot; F6P_dot; %  5
    FBP_dot; GAP_dot; PEP_dot; PYR_dot; AcCoA_dot;    % 10
    AcP_dot; ACEex_dot; ICIT_dot; aKG_dot; SUC_dot;   % 15
    FUM_dot; MAL_dot; OAA_dot; GOX_dot; sixPGL_dot;   % 20
    sixPG_dot; KDPG_dot; RU5P_dot; R5P_dot; X5P_dot;  % 25
    S7P_dot; E4P_dot; cAMP_dot; EIIAP_dot; Glk_dot;   % 30
    Pfk_dot; Fbp_dot; Fba_dot; GAPDH_dot; Pyk_dot;    % 35
    Pps_dot; PDH_dot; Acs_dot; CS_dot; ICDH_dot;      % 40
    ICDHP_dot; aKGDH_dot; SDH_dot; Fum_dot; MDH_dot;  % 45
    MaeB_dot; Pck_dot; Ppc_dot; Icl_dot; MS_dot;       % 50
    AceK_dot ];                                       % 51
