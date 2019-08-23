function Flux = KurataModel_Flux(t,y)

global strain_no continuous_flg ;

X      = y(1);
GLCex  = y(2);
GLC    = y(3);
G6P    = y(4);
F6P    = y(5);
FBP    = y(6);
GAP    = y(7);
PEP    = y(8);
PYR    = y(9);
AcCoA  = y(10);
AcP    = y(11);
ACEex  = y(12);
ICIT   = y(13);
aKG    = y(14);
SUC    = y(15);
FUM    = y(16);
MAL    = y(17);
OAA    = y(18);
GOX    = y(19);
sixPGL = y(20);
sixPG  = y(21);
KDPG   = y(22);
RU5P   = y(23);
R5P    = y(24);
X5P    = y(25);
S7P    = y(26);
E4P    = y(27);
cAMP   = y(28);
EIIAP  = y(29);
Glk    = y(30);
Pfk    = y(31);
Fbp    = y(32);
Fba    = y(33);
Gapdh  = y(34);
Pyk    = y(35);
Pps    = y(36);
Pdh    = y(37);
Acs    = y(38);
Cs     = y(39);
Icdh   = y(40);
IcdhP  = y(41);
akgdh  = y(42);
Sdh    = y(43);
Fum    = y(44);
Mdh    = y(45);
MaeB   = y(46);
Pck    = y(47);
Ppc    = y(48);
Icl    = y(49);
Ms     = y(50);
AceK   = y(51);

%********** MODEL PARAMETERS
% constant component 			
GLCfeed    	 =	2.22E+01	;
ATP        	 =	9.60E+00	;
ADP        	 =	5.60E-01	;
AMP        	 =	2.80E-01	;
NAD        	 =	2.60E+00	;
NADH       	 =	8.30E-02	;
NADP       	 =	2.10E-03	;
NADPH      	 =	1.20E-01	;
CoA        	 =	1.40E+00	;
Pi         	 =	1.00E+01	;
EIIAtotal  	 =	7.69E-02	;
Crptotal   	 =	1.15E-02	;
Cratotal   	 =	3.00E-04	;
PdhRtotal  	 =	6.66E-05	;
IclRtotal  	 =	8.30E-05	;
pH         	 =	7.50E+00	;
aceBAK_DNA 	 =	5.15E-07	;
% Kinetic parameter for PTS1 & PTS4			
kPts1       	 =	6.60E+04	;
kmPts1      	 =	2.77E+04	;
vPts4_max   	 =	3.94E+03	;
KPts_EIIA   	 =	2.11E-03	;
KPts_GLC    	 =	4.90E-03	;
% Kinetic parameter for Non-PTS			
vNonpts_max 	 =	4.00E+03	;
KNonpts_S   	 =	1.55E+00	;
KNonpts_I   	 =	1.00E-02	;
kGlk_cat    	 =	2.01E+07	;
KGlk_GLC_m  	 =	1.54E-01	;
KGlk_ATP_m  	 =	7.12E-01	;
KGlk_G6P_i  	 =	1.50E+01	;
% Kinetic parameter for Pgi			
vPgi_max        	 =	3.56E+06	;
KPgi_eq         	 =	1.44E+00	;
KPgi_G6P        	 =	2.46E+00	;
KPgi_F6P        	 =	3.42E-01	;
KPgi_F6P_6pginh 	 =	1.99E-01	;
KPgi_G6P_6pginh 	 =	1.81E-01	;
% Kinetic parameter for Pfk			
kPfk_cat   	 =	2.49E+11	;
KPfk_PEP   	 =	1.74E+00	;
KPfk_ADP_b 	 =	2.56E-01	;
KPfk_AMP_b 	 =	2.55E-02	;
KPfk_ADP_a 	 =	2.77E+02	;
KPfk_AMP_a 	 =	1.01E+01	;
KPfk_ATP_s 	 =	1.61E-01	;
KPfk_ADP_c 	 =	4.48E-01	;
KPfk_F6P_s 	 =	2.11E-02	;
LPfk       	 =	1.77E+06	;
nPfk       	 =	4.00E+00	;
% Kinetic parameter for Fbp			
kFbp_cat 	 =	7.86E+07    ;
KFbp_FBP 	 =	8.92E-03	;
KFbp_PEP 	 =	4.88E-01	;
LFbp     	 =	4.41E+06	;
nFbp     	 =	4.00E+00	;
% Kinetic parameter for Fba			
kFba_cat     	 =	6.95E+08	;
KFba_eq      	 =	3.72E-01	;
KFba_FBP     	 =	8.39E-02	;
KFba_GAP     	 =	1.54E-01	;
KFba_DHAP    	 =	8.80E-02	;
VFba_blf     	 =	1.54E+00	;
KFba_GAP_inh 	 =	6.00E-01	;
% Kinetic parameter for Gapdh			
kGapdh_cat  	 =	5.05E+09	;
KGapdh_eq   	 =	3.00E-01	;
KGapdh_GAP  	 =	1.52E-01	;
KGapdh_PGP  	 =	1.34E-01	;
KGapdh_NAD  	 =	4.50E-01	;
KGapdh_NADH 	 =	2.08E-02	;
% Kinetic parameter for Pyk			
kPyk_cat 	 =	8.16E+05	;
KPyk_PEP 	 =	3.10E-01	;
KPyk_FBP 	 =	2.53E-01	;
KPyk_AMP 	 =	2.55E-01	;
KPyk_ADP 	 =	2.10E-01	;
KPyk_ATP 	 =	2.02E+01	;
LPyk     	 =	9.97E+02	;
nPyk     	 =	4.00E+00	;
% Kinetic parameter for Pps			
kPps_cat 	 =	2.49E+05	;
KPps_PYR 	 =	7.13E-04	;
KPps_PEP 	 =	2.16E-04	;
LPps     	 =	1.04E-79	;
nPps     	 =	2.00E+00	;
% Kinetic parameter for Pdh			
kPdh_cat     	 =	4.83E+07	;
KPdh_i       	 =	6.83E+01	;
KPdh_PYR_m   	 =	1.00E+00	;
KPdh_NAD_m   	 =	4.01E-01	;
KPdh_NADH_m  	 =	4.71E-02	;
KPdh_CoA_m   	 =	4.82E-03	;
KPdh_AcCoA_m 	 =	8.00E-03	;
% Kinetic parameter for Pta			
vPta_max     	 =	5.36E+03	;
KPta_eq      	 =	2.84E-02	;
KPta_AcCoA_i 	 =	2.01E-01	;
KPta_CoA_i   	 =	8.10E-02	;
KPta_Pi_m    	 =	6.87E-01	;
KPta_Pi_i    	 =	2.11E+00	;
KPta_AcP_m   	 =	2.29E-01	;
KPta_AcP_i   	 =	3.20E-01	;
% Kinetic parameter for Ack			
vAck_max   	 =	1.95E+05	;
KAck_eq    	 =	2.34E+02	;
KAck_ADP_m 	 =	1.77E-01	;
KAck_AcP_m 	 =	4.63E-02	;
KAck_ACE_m 	 =	6.10E+00	;
KAck_ATP_m 	 =	8.86E-02	;
% Kinetic parameter for Acs			
kAcs_cat 	 =	1.29E+05	;
KAcs_ACE 	 =	2.40E-02	;
% Kinetic parameter for Cs			
kCs_cat       	 =	5.43E+06	;
KCs_aKG       	 =	1.87E-01	;
KCs_OAA_AcCoA 	 =	8.86E-05	;
KCs_AcCoA     	 =	2.98E-02	;
KCs_OAA       	 =	1.76E-03	;
% Kinetic parameter for Icdh			
%%%kIcdh_cat  	 =	8.58E+06 *4	;
kIcdh_cat  	 =	3.432E+07	;
KIcdh_ICIT 	 =	2.01E-04	;
KIcdh_PEP  	 =	5.48E-02	;
LIcdh      	 =	9.26E+01	;
nIcdh      	 =	2.00E+00	;
% Kinetic parameter for alphakgdh			
%%%kakgdh_cat    	 =	7.02E+08*0.001	;
kakgdh_cat    	 =	7.02E+05	;
Kakgdh_NAD_m  	 =	5.52E-02	;
Kakgdh_CoA_m  	 =	3.10E-03	;
Kakgdh_aKG_m  	 =	2.37E-01	;
Kakgdh_Z      	 =	4.17E+00	;
Kakgdh_SUC_I  	 =	2.07E+00	;
Kakgdh_NADH_I 	 =	1.79E-02	;
Kakgdh_aKG_I  	 =	1.00E+00	;
% Kinetic parameter for Sdh   			
kSdh1_cat  	 =	1.96E+06	;
kSdh2_cat  	 =	1.96E+06	;
KSdh_eq    	 =	2.97E+01	;
KSdh_SUC_m 	 =	1.69E-01	;
% Kinetic parameter for FUM 			
kFum1_cat  	 =	3.82E+06 	;
kFum2_cat  	 =	3.82E+06	;
KFum_eq    	 =	1.29E+01	;
KFum_FUM_m 	 =	9.65E-02	;
% Kinetic parameter for Mdh 			
%%%kMdh1_cat   	 =	4.60E+06 *5	;
kMdh1_cat   	 =	2.30E+07	;
%%%kMdh2_cat   	 =	4.60E+06 *5	;
kMdh2_cat   	 =	2.30E+07	;
KMdh_eq     	 =	7.50E-01	;
KMdh_NAD_m  	 =	5.22E-02	;
KMdh_NAD_I  	 =	1.14E-01	;
KMdh_NAD_II 	 =	9.45E-01	;
KMdh_MAL_m  	 =	1.77E-01	;
KMdh_MAL_I  	 =	3.46E-01	;
KMdh_OAA_m  	 =	2.08E-01	;
KMdh_OAA_I  	 =	3.44E-01	;
KMdh_OAA_II 	 =	7.32E-02	;
KMdh_NADH_m 	 =	3.13E-02	;
KMdh_NADH_I 	 =	2.48E-02	;
% Kinetic parameter for MaeB			
%%%kMaeB_cat   	 =	3.68E+05 *50	;
kMaeB_cat   	 =	1.84E+07 	;
KMaeB_MAL   	 =	1.46E-03	;
KMaeB_AcCoA 	 =	1.82E+00	;
KMaeB_cAMP  	 =	4.55E+00	;
LMaeB       	 =	2.39E+05	;
nMaeB       	 =	1.99E+00	;
% Kinetic parameter for Pck			
kPck_cat   	 =	3.45E+07	;
KPck_OAA   	 =	5.78E-01	;
KPck_ATP_i 	 =	4.01E-02	;
KPck_ADP_i 	 =	1.98E-02	;
KPck_PEP   	 =	7.01E-02	;
KPck_PEP_i 	 =	6.00E-02	;
KPck_OAA_I 	 =	3.47E-01	;
KPck_ATP_I 	 =	3.98E-02	;
% Kinetic parameter for Ppc			
%%%kPpc_cat 	 =	9.702E+04 * 5	;
kPpc_cat 	 =	4.851E+05	;
KPpc_PEP 	 =	2.72E-02	;
KPpc_FBP 	 =	1.89E-01	;
LPpc     	 =	5.65E+06	;
nPpc     	 =	5.00E+00	;
% Kinetic parameter for Icl			
kIcl_cat  	 =	9.42E+07	;
%%%KIcl_ICIT 	 =	1.51E-02 * 0.5	;
KIcl_ICIT 	 =	7.55E-03	;
KIcl_PEP  	 =	1.63E-02	;
KIcl_3PG  	 =	5.22E-01	;
KIcl_aKG  	 =	8.34E-01	;
LIcl      	 =	1.91E+05	;
nIcl      	 =	4.00E+00	;
% Kinetic parameter for Ms			
kMs_cat       	 =	1.70E+10	;
KMs_GOX_AcCoA 	 =	3.96E-02	;
KMs_AcCoA     	 =	4.62E-02	;
KMs_GOX       	 =	1.11E+00	;
% Kinetic parameter for AceK-ki & AceK-ph			
kAceKki_cat 	 =	1.00E+16	;
kAceKph_cat 	 =	6.00E+2	;
KAceK_Icdh  	 =	1.97E-01	;
KAceK_IcdhP 	 =	7.26E+00	;
KAceK_OAA  	 =	1.00E-02	;
LAceK      	 =	2.83E+08	;
nAceK      	 =	5.00E+00	;
KAceKph_OAA  =  1.00E-02;
LAceKph      =  7.00E+4  ;
nAceKph      =  5.00E+00	; 
% Kinetic parameter for G6Pdh			
vG6pdh_max           	 =	1.66E+04	;
%%%KG6pdh_G6P           	 =	7.12E-02 * 5	;
KG6pdh_G6P           	 =	3.56E-01 	;
KG6pdh_NADPH_g6pinh  	 =	1.93E-01	;
KG6pdh_NADP          	 =	3.82E-03	;
KG6pdh_NADPH_nadpinh 	 =	8.15E-02	;
% Kinetic parameter for Pgl			
vPgl_max    	 =	2.28E+04	;
KPgl_eq     	 =	4.27E+01	;
KPgl_6PGL_m 	 =	2.30E-02	;
KPgl_6PG_m  	 =	1.00E+01	;
KPgl_h1     	 =	4.37E-03	;
KPgl_h2     	 =	9.70E-06	;
% Kinetic parameter for Edd			
vEdd_max    	 =	5.15E+02	;
KEdd_eq     	 =	1.00E+03	;
KEdd_6PG_m  	 =	1.18E-01	;
KEdd_KDPG_m 	 =	2.02E+00	;
pH_Edd_m    	 =	7.53E+00	;
pK_Edd      	 =	8.74E+00	;
% Kinetic parameter for Eda			
vEda_max    	 =	6.67E+02	;
KEda_eq     	 =	5.02E-01	;
KEda_PYR_m  	 =	7.70E+00	;
KEda_KDPG_m 	 =	1.50E-01	;
KEda_GAP_m  	 =	1.18E+00	;
pH_Eda_m    	 =	1.04E+01	;
pK_Eda      	 =	3.70E+01	;
% Kinetic parameter for 6PGDH			
v6Pgdh_max       	 =	2.48E+05	;
K6Pgdh_6PG       	 =	1.01E-01	;
K6Pgdh_NADP      	 =	2.08E-02	;
K6Pgdh_NADPH_inh 	 =	4.21E-02	;
K6Pgdh_ATP_inh   	 =	3.01E+00	;
% Kinetic parameter for R5PI			
vR5pi_max 	 =	3.21E+04	;
KR5pi_eq  	 =	4.77E-01	;
% Kinetic parameter for Rpe			
vRpe_max 	 =	1.27E+04	;
KRpe_eq  	 =	1.41E+00	;
% Kinetic parameter for Tkt1			
vTkt1_max 	 =	8.87E+03	;
KTkt1_eq  	 =	1.20E+00	;
% Kinetic parameter for Tkt2			
vTkt2_max 	 =	3.80E+05	;
KTkt2_eq  	 =	9.97E+00	;
% Kinetic parameter for Tal			
vTal_max 	 =	7.17E+04	;
KTal_eq  	 =	1.05E+00	;
% Kinetic parameter for Cya			
vCya_max   	 =	9.45E+00	;
KCya_EIIAP 	 =	2.08E-03	;
% Kinetic parameter for cAMPdegr 			
vcAMPdegr_max  	 =	9.21E+00	;
KcAMPdegr_cAMP 	 =	4.83E-02	;
% Kinetic parameter for Glk			
Kglk_Cra         	 =	1.22E-08	;
vglk_Cra_unbound 	 =	1.50E-01	; 
vglk_Cra_bound   	 =	1.13E-03	; 
% Kinetic parameter for pfkA			
KpfkA_Cra         	 =	9.87E-09	;
vpfkA_Cra_unbound 	 =	3.75E-02	;
vpfkA_Cra_bound   	 =	1.03E-03	;
% Kinetic parameter for fbp			
Kfbp_Cra         	 =	3.75E-05	;
vfbp_Cra_unbound 	 =	0.00E+00	;
vfbp_Cra_bound   	 =	1.03E-03	;
% Kinetic parameter for fbaA			
KfbaA_Cra         	 =	3.26E-03	;
vfbaA_Cra_unbound 	 =	1.83E-03	;
vfbaA_Cra_bound   	 =	0.00E+00	;
KfbaA_Crp         	 =	9.48E-03	;
vfbaA_Crp_unbound 	 =	0.00E+00	;
vfbaA_Crp_bound   	 =	1.13E-03	;
% Kinetic parameter for gapA			
KgapA_Cra         	 =	1.61E-02	;
vgapA_Cra_unbound 	 =	1.60E-03	;
vgapA_Cra_bound   	 =	0.00E+00	;
KgapA_Crp         	 =	4.75E-02	;
vgapA_Crp_unbound 	 =	0.00E+00	;
vgapA_Crp_bound   	 =	2.18E-03	;
% Kinetic parameter for pykF			
KpykF_Cra         	 =	7.26E-05	;
vpykF_Cra_unbound 	 =	8.20E-03	;
vpykF_Cra_bound   	 =	3.12E-05	;
% Kinetic parameter for ppsA			
KppsA_Cra         	 =	4.88E-04	;
vppsA_Cra_unbound 	 =	0.00E+00	;
vppsA_Cra_bound   	 =	3.86E-03	;
% Kinetic parameter for lpd			
Klpd_PdhR         	 =	2.46E-05	;
vlpd_PdhR_unbound 	 =	9.29E-03	;
vlpd_PdhR_bound   	 =	5.45E-05	;
% Kinetic parameter for acs 			
Kacs_Crp         	 =	1.36E-03	;
vacs_Crp_unbound 	 =	0.00E+00	;
vacs_Crp_bound   	 =	2.62E-03	;
nacs             	 =	2.31E+00	;
% Kinetic parameter for gltA 			
KgltA_Crp         	 =	5.65E-02	;
vgltA_Crp_unbound 	 =	0.00E+00	;
vgltA_Crp_bound   	 =	2.05E-02	;
ngltA             	 =	1.07E+00	;
% Kinetic parameter for icdA 			
KicdA_Cra         	 =	2.92E-05	;
vicdA_Cra_unbound 	 =	5.13E-04	;
vicdA_Cra_bound   	 =	1.42E-03	;
% Kinetic parameter for sucAB 			
KsucAB_Crp         	 =	3.09E-01	;
vsucAB_Crp_unbound 	 =	0.00E+00	;
vsucAB_Crp_bound   	 =	1.74E-01	;
nsucAB             	 =	7.40E-01	;
% Kinetic parameter for sdhCDAB 			
KsdhCDAB_Crp         	 =	8.56E-02	;
vsdhCDAB_Crp_unbound 	 =	0.00E+00	;
vsdhCDAB_Crp_bound   	 =	1.77E-02	;
nsdhCDAB             	 =	7.40E-01	;
% Kinetic parameter for fumABC 			
KfumABC_Crp         	 =	1.49E-01	;
vfumABC_Crp_unbound 	 =	0.00E+00	;
vfumABC_Crp_bound   	 =	3.76E-02	;
nfumABC             	 =	7.40E-01	;
% Kinetic parameter for mdh			
Kmdh_Crp         	 =	2.14E-01	;
vmdh_Crp_unbound 	 =	0.00E+00	;
vmdh_Crp_bound   	 =	7.33E-02	;
% Kinetic parameter for maeB 			
SS_MaeB 	 =	7.50E-03 	;
% Kinetic parameter for pckA 	 =	0.00E+00	;
KpckA_Cra         	 =	1.52E-04	;
vpckA_Cra_unbound 	 =	0.00E+00	;
vpckA_Cra_bound   	 =	9.95E-03	;
% Kinetic parameter for ppc 			
SS_Ppc 	 =	7.40E-03 	;
% Kinetic parameter for aceA     			
KaceBAK_Cra         	 =	5.40E-04	;
vaceBAK_Cra_unbound 	 =	1.27E-05	;
vaceBAK_Cra_bound   	 =	3.70E-03	;
KaceBAK_Crp         	 =	5.26E-01	;
vaceBAK_Crp_unbound 	 =	2.51E-04	;
vaceBAK_Crp_bound   	 =	1.58E-06	;
KaceBAK_DNA         	 =	9.51E-07	;
KaceBAK_PYR         	 =	2.01E+00	;
KaceBAK_PYRprime    	 =	2.72E-03	;
KaceBAK_GOX         	 =	2.25E-03	;
kaceBAK_cat_IclR    	 =	3.20E-01	;
LaceBAK             	 =	4.12E+02	;
% Kinetic parameter for aceB   			
Factor_aceB         	 =	3.09E-01	;
% Kinetic parameter for aceK   	 =	0.00E+00	;
Factor_aceK         	 =	1.72E-02	;

% protein degradtion rate constant
kdegr 	 =	3.30E-01	;
% Kinetic parameter for Biomass 			
kBM_GLC_G6P   	 =	4.89E+01	;
kBM_GLC_F6P   	 =	1.27E+03	;
kBM_GLC_GAP   	 =	1.02E+02	;
kBM_GLC_PEP   	 =	9.57E+02	;
kBM_GLC_PYR   	 =	8.51E+02	;
%%%kBM_GLC_AcCoA 	 =	2.37E+03*0.7	;
kBM_GLC_AcCoA 	 =	1.659E+03	;
%%%kBM_GLC_aKG   	 =	3.00E+03*0.1	;
%%%kBM_GLC_SUC   	 =	1.90E+03*0.1	;
kBM_GLC_aKG   	 =	3.00E+02	;
kBM_GLC_SUC   	 =	1.90E+02	;
kBM_GLC_FUM   	 =	3.47E+03	;
kBM_GLC_OAA   	 =	2.31E+04	;
kBM_GLC_R5P   	 =	3.08E+02	;
kBM_GLC_E4P   	 =	1.51E+03	;
% Kinetic parameter for Crp-cAMP			
KCrpcAMP 	 =	4.20E-01	;
nCrpcAMP 	 =	1.00E+00	;
% Kinetic parameter for Cra-FBP			
KCraFBP 	 =	2.94E-02	;
nCraFBP 	 =	2.00E+00	;
% Kinetic parameter for PdhR-PYR			
KPdhRPYR 	 =	4.34E-02	;
nPdhRPYR 	 =	1.00E+00	;
%			
POratio         =	3.48E+00	;
POratio_prime 	=	1.49E+00	;
kATP 	=	1.26E-05	;
rho 	=	5.64E+02	;
D       =	0.00E+00	;

%%

% WT
if strain_no == 1
end

% Delta pykA/pykF (Pyk)
if strain_no == 26
    vpykF_Cra_unbound = 0; %  PykA/PykF
    vpykF_Cra_bound   = 0;
end

% Delta pykF
if strain_no == 13
    vpykF_Cra_unbound = 0.5*vpykF_Cra_unbound; %  PykA/PykF
    vpykF_Cra_bound   = 0.5*vpykF_Cra_bound;
end

% Delta pgi
if strain_no == 4
    vPgi_max = 0; 
end

% Delta ppc
if strain_no == 25
    SS_Ppc = 0;
end

% Delta glk
if strain_no == 2
    vglk_Cra_unbound = 0;
    vglk_Cra_bound   = 0;
end

% Delta pfkA
if strain_no == 5
    vpfkA_Cra_unbound = 0.5 * vpfkA_Cra_unbound; % PkA/PfkB
    vpfkA_Cra_bound   = 0.5 * vpfkA_Cra_bound;
end

% Delta pfkB
if strain_no == 6
    vpfkA_Cra_unbound = 0.5 * vpfkA_Cra_unbound; 
    vpfkA_Cra_bound   = 0.5 * vpfkA_Cra_bound;
end

% Delta fbp
if strain_no == 7
    vfbp_Cra_unbound = 0;
    vfbp_Cra_bound   = 0;
end

% Delta fbaB
if strain_no == 8
    vfbaA_Cra_unbound = 0.5* vfbaA_Cra_unbound; %FbaA/FbaB
    vfbaA_Cra_bound   = 0.5* vfbaA_Cra_bound;
    vfbaA_Crp_unbound = 0.5* vfbaA_Crp_unbound;
    vfbaA_Crp_bound   = 0.5* vfbaA_Crp_bound;
end

% Delta gpmA
if strain_no == 10
    vgapA_Cra_unbound =0;
    vgapA_Crp_bound=0;
end

%Delta pykA
if strain_no == 12
    vpykF_Cra_unbound = 0.5*vpykF_Cra_unbound; %  PykA/PykF
    vpykF_Cra_bound   = 0.5*vpykF_Cra_bound;
end

% Delta gpmA
if strain_no == 14
    vppsA_Cra_unbound = 0;
    vppsA_Cra_bound   = 0;
end

% Delta zwf(G6Pdh)
if strain_no == 15
    vG6pdh_max = 0;
end

% Delta pgl
if strain_no == 16
    vPgl_max = 0;
end

% Delta gnd (6PGDH)
if strain_no == 17
    v6Pgdh_max = 0 *v6Pgdh_max;
end

% Delta rpe 
if strain_no == 18
    vRpe_max = 0;
end

% Delta rpiA (R5PI)
if strain_no == 19
    vR5pi_max = 0.5* vR5pi_max; % RpiA/RpiB isozymes
end

% Delta rpiB (R5PI)
if strain_no == 20
    vR5pi_max = 0.5* vR5pi_max;
end

% Delta tktA
if strain_no == 21
    vTkt1_max = 0.5*vTkt1_max; % TktA/TktB isozymes
    vTkt2_max = 0.5*vTkt2_max;
end

% Delta tktB
if strain_no == 22
    vTkt1_max = 0.5*vTkt1_max; 
    vTkt2_max = 0.5*vTkt2_max;
end

% Delta talA
if strain_no == 23
    vTal_max = 0.5 * vTal_max ; % TalA/TalB isozymes
end

% Delta talB
if strain_no == 24
    vTal_max = 0.5 * vTal_max ; 
end

%%%%%%%%%%%%%%%%%%%
% additions from Denis
% Delta AKGDH
if strain_no == 26
    kakgdh_cat = 0;
end

% Delta Pts
if strain_no == 27
    vPts4_max = 0;
end

% Delta sdhC
if strain_no == 28
    kSdh1_cat = 0;
end

% Delta tpi
if strain_no == 29
    kGapdh_cat = 0;
end

% zwf(15) overexpression
if strain_no == 30
    vG6pdh_max = v6Pgdh_max * 15;
end

% pgi(0)
if strain_no == 31
    vPgi_max = vPgi_max * 0.2;
end

% pgi(20)
if strain_no == 32
    vPgi_max = vPgi_max * 1.2;
end

%pgi(50)
if strain_no == 33
    vPgi_max = vPgi_max * 2.4;
end

%pgi(100)
if strain_no == 34
    vPgi_max = vPgi_max * 4.1;
end

% eno(0)
if strain_no == 36
    kGapdh_cat = kGapdh_cat * 0.2;
end

% eno(50)
if strain_no == 37
    kGapdh_cat = kGapdh_cat * 1.8;
end

% eno(200)
if strain_no == 38
    kGapdh_cat = kGapdh_cat * 3.0;
end

% eno(500)
if strain_no == 39
    kGapdh_cat = kGapdh_cat * 3.1;
end

% delta eda
if strain_no == 40
    vEda_max = 0;
end

% delta edd
if strain_no == 41
    vEdd_max = 0;
end

%% Batch culture
if continuous_flg == 0
    GLCfeed = 0;
    D = 0;
end

%% Continuous culture
if continuous_flg ~= 0
    GLCfeed = 22.2;
end
if continuous_flg == 1
    D = 0.2;
end
if continuous_flg == 2
    D = 0.1;
end
if continuous_flg == 3
    D = 0.4;
end
if continuous_flg == 4
    D = 0.5;
end
if continuous_flg == 5
    D = 0.7;
end
if continuous_flg == 6
    D = 0.25;
end


%%
%********** MODEL VARIABLES (v-independent)
alphaGLC=1; alphaACE=0;

H        = power(10,-pH)*1e+3;
EIIA     = EIIAtotal - EIIAP;
CrpcAMP  = Crptotal*power(cAMP,nCrpcAMP)/(power(cAMP,nCrpcAMP)+power(KCrpcAMP,nCrpcAMP));
Crp      = Crptotal - CrpcAMP;
CraFBP   = Cratotal*power(FBP,nCraFBP)/(power(FBP,nCraFBP)+power(KCraFBP,nCraFBP));
Cra      = Cratotal - CraFBP;
PdhRPYR  = PdhRtotal*power(PYR,nPdhRPYR)/(power(PYR,nPdhRPYR)+power(KPdhRPYR,nPdhRPYR));
PdhR     = PdhRtotal - PdhRPYR;
icdh_b    = 6E-7; 

% growth-related constants
%vG_aceA growth-suppressed function
n_mu=4;
Kmu=0.60;

% growth inhibition by F6P
F6PWT=0.48; % critical
growth_sup = (F6P-F6PWT)/F6PWT;
if growth_sup > 0
    kATP = 1.26E-05 * (1- 0.8*growth_sup);
else   
    kATP = 1.26E-05;
end

%% Flux equation
vPts1          = kPts1*PEP*EIIA - kmPts1*PYR*EIIAP;
vPts4          = vPts4_max*EIIAP*GLCex/((KPts_EIIA+EIIAP)*(KPts_GLC+GLCex));
vPts4_medium   = vPts4*X/rho;
vNonpts        = vNonpts_max*GLCex/(KNonpts_S+(1+EIIA/KNonpts_I)*GLCex);
vNonpts_medium = vNonpts*X/rho;
vE_Glk         = Glk*kGlk_cat*(GLC/KGlk_GLC_m)*(ATP/(KGlk_ATP_m*(1+G6P/KGlk_G6P_i)))/(1+GLC/KGlk_GLC_m+ATP/(KGlk_ATP_m*(1+G6P/KGlk_G6P_i))+GLC*ATP/(KGlk_GLC_m*KGlk_ATP_m*(1+G6P/KGlk_G6P_i))+G6P/KGlk_G6P_i);
vE_Pgi         = vPgi_max*(G6P-F6P/KPgi_eq)/(KPgi_G6P*(1+F6P/(KPgi_F6P*(1+sixPG/KPgi_F6P_6pginh))+sixPG/KPgi_G6P_6pginh)+G6P);
 A             = 1 + PEP/KPfk_PEP + ADP/KPfk_ADP_b + AMP/KPfk_AMP_b;
 B             = 1 + ADP/KPfk_ADP_a + AMP/KPfk_AMP_a;
vE_Pfk         = Pfk*kPfk_cat*ATP*F6P/((ATP+KPfk_ATP_s*(1+ADP/KPfk_ADP_c))*(F6P+KPfk_F6P_s*A/B)*(1+LPfk/power(1+F6P*B/(KPfk_F6P_s*A),nPfk)));
vE_Fbp         = Fbp*kFbp_cat*FBP/KFbp_FBP*power(1+FBP/KFbp_FBP,nFbp-1)/(power(1+FBP/KFbp_FBP,nFbp)+LFbp/power(1+PEP/KFbp_PEP,nFbp));
vE_Fba         = Fba*kFba_cat*(FBP-power(GAP,2)/KFba_eq)/(KFba_FBP+FBP+KFba_GAP*GAP/(KFba_eq*VFba_blf)+KFba_DHAP*GAP/(KFba_eq*VFba_blf)+FBP*GAP/KFba_GAP_inh+power(GAP,2)/(KFba_eq*VFba_blf));
vE_Gapdh       = Gapdh*kGapdh_cat*(GAP*NAD-PEP*NADH/KGapdh_eq)/((KGapdh_GAP*(1+PEP/KGapdh_PGP)+GAP)*(KGapdh_NAD*(1+NADH/KGapdh_NADH)+NAD));
vE_Pyk         = Pyk*kPyk_cat*PEP*power(PEP/KPyk_PEP+1,nPyk-1)*ADP/(KPyk_PEP*(LPyk*power((1+ATP/KPyk_ATP)/(FBP/KPyk_FBP+AMP/KPyk_AMP+1),nPyk)+power(PEP/KPyk_PEP+1,nPyk))*(ADP+KPyk_ADP));
vE_Pps         = Pps*kPps_cat*PYR/KPps_PYR*power(1+PYR/KPps_PYR,nPps-1)/(power(1+PYR/KPps_PYR,nPps)+LPps*power(1+PEP/KPps_PEP,nPps));
vE_Pdh         = Pdh*kPdh_cat*(1/(1+KPdh_i*NADH/NAD))*(PYR/KPdh_PYR_m)*(NAD/KPdh_NAD_m)*(CoA/KPdh_CoA_m)/((1+PYR/KPdh_PYR_m)*(1+NAD/KPdh_NAD_m+NADH/KPdh_NADH_m)*(1+CoA/KPdh_CoA_m+AcCoA/KPdh_AcCoA_m));
vE_Pta         = vPta_max*(1/(KPta_AcCoA_i*KPta_Pi_m))*(AcCoA*Pi-AcP*CoA/KPta_eq)/(1+AcCoA/KPta_AcCoA_i+Pi/KPta_Pi_i+AcP/KPta_AcP_i+CoA/KPta_CoA_i+AcCoA*Pi/(KPta_AcCoA_i*KPta_Pi_m)+AcP*CoA/(KPta_AcP_m*KPta_CoA_i));
vE_Ack         = vAck_max*(1/(KAck_ADP_m*KAck_AcP_m))*(AcP*ADP-ACEex*ATP/KAck_eq)/((1+AcP/KAck_AcP_m+ACEex/KAck_ACE_m)*(1+ADP/KAck_ADP_m+ATP/KAck_ATP_m));
vE_Ack_medium  = vE_Ack*X/rho;
vE_Acs         = Acs*kAcs_cat*ACEex/(ACEex+KAcs_ACE);
vE_Acs_medium  = vE_Acs*X/rho;
vE_Cs          = Cs*kCs_cat*OAA*AcCoA/((1+aKG/KCs_aKG)*KCs_OAA_AcCoA*KCs_AcCoA+KCs_AcCoA*OAA+(1+aKG/KCs_aKG)*KCs_OAA*AcCoA+OAA*AcCoA);
vE_Icdh        = (Icdh + icdh_b)*kIcdh_cat*ICIT/KIcdh_ICIT*power(1+ICIT/KIcdh_ICIT,nIcdh-1)/(power(1+ICIT/KIcdh_ICIT,nIcdh)+LIcdh*power(1+PEP/KIcdh_PEP,nIcdh));
vE_akgdh       = akgdh*kakgdh_cat*aKG*CoA*NAD/(Kakgdh_NAD_m*aKG*CoA+Kakgdh_CoA_m*aKG*NAD+Kakgdh_aKG_m*CoA*NAD+aKG*CoA*NAD+Kakgdh_aKG_m*Kakgdh_Z*SUC*NADH/Kakgdh_SUC_I+Kakgdh_NAD_m*aKG*CoA*NADH/Kakgdh_NADH_I+Kakgdh_CoA_m*aKG*NAD*SUC/Kakgdh_SUC_I+Kakgdh_aKG_m*Kakgdh_Z*aKG*SUC*NADH/(Kakgdh_aKG_I*Kakgdh_SUC_I));
vSdh1_max      = Sdh*kSdh1_cat;
vSdh2_max      = Sdh*kSdh2_cat;
vE_Sdh         = vSdh1_max*vSdh2_max*(SUC-FUM/KSdh_eq)/(KSdh_SUC_m*vSdh2_max+vSdh2_max*SUC+vSdh1_max*FUM/KSdh_eq);
vFum1_max      = Fum*kFum1_cat;
vFum2_max      = Fum*kFum2_cat;
vE_Fum         = vFum1_max*vFum2_max*(FUM-MAL/KFum_eq)/(KFum_FUM_m*vFum2_max+vFum2_max*FUM+vFum1_max*MAL/KFum_eq);
vMdh1_max      = Mdh*kMdh1_cat;
vMdh2_max      = Mdh*kMdh2_cat;
vE_Mdh         = vMdh1_max*vMdh2_max*(NAD*MAL-NADH*OAA/KMdh_eq)/(KMdh_NAD_I*KMdh_MAL_m*vMdh2_max+KMdh_MAL_m*vMdh2_max*NAD+KMdh_NAD_m*vMdh2_max*MAL+ vMdh2_max*NAD*MAL+KMdh_OAA_m*vMdh1_max*NADH/KMdh_eq+KMdh_NADH_m*vMdh1_max*OAA/KMdh_eq+ vMdh1_max*NADH*OAA/KMdh_eq+vMdh1_max*KMdh_OAA_m*NAD*NADH/(KMdh_eq*KMdh_NAD_I)+vMdh2_max*KMdh_NAD_m*MAL*OAA/KMdh_OAA_I+ vMdh2_max*NAD*MAL*NADH/KMdh_NADH_I+vMdh1_max*MAL*NADH*OAA/(KMdh_eq*KMdh_MAL_I)+vMdh2_max*NAD*MAL*OAA/KMdh_OAA_II+ vMdh1_max*NAD*NADH*OAA/(KMdh_NAD_II*KMdh_eq)+KMdh_NAD_I*vMdh2_max*NAD*MAL*NADH*OAA/(KMdh_NAD_II*KMdh_OAA_m*KMdh_NADH_I));
vE_MaeB        = MaeB*kMaeB_cat*MAL/KMaeB_MAL*power(1+MAL/KMaeB_MAL,nMaeB-1)/(power(1+MAL/KMaeB_MAL,nMaeB)+LMaeB*power(1+AcCoA/KMaeB_AcCoA+cAMP/KMaeB_cAMP,nMaeB));
vE_Pck         = Pck*kPck_cat*OAA*ATP/ADP/(KPck_OAA*ATP/ADP+OAA*ATP/ADP+KPck_ATP_i*KPck_OAA/KPck_ADP_i+KPck_ATP_i*KPck_OAA/(KPck_PEP*KPck_ADP_i)*PEP+KPck_ATP_i*KPck_OAA/(KPck_PEP_i*KPck_ATP_I)*ATP/ADP*PEP+KPck_ATP_i*KPck_OAA/(KPck_ADP_i*KPck_OAA_I)*OAA);
vE_Ppc         = Ppc*kPpc_cat*PEP/KPpc_PEP*power(1+PEP/KPpc_PEP,nPpc-1)/(power(1+PEP/KPpc_PEP,nPpc)+LPpc/power(1+FBP/KPpc_FBP,nPpc));
vE_Icl         = Icl*kIcl_cat*ICIT/KIcl_ICIT*power(1+ICIT/KIcl_ICIT,nIcl-1)/(power(1+ICIT/KIcl_ICIT,nIcl)+LIcl*power(1+PEP/KIcl_PEP+GAP/KIcl_3PG+aKG/KIcl_aKG,nIcl));
vE_Ms          = Ms*kMs_cat*GOX*AcCoA/(KMs_GOX_AcCoA*KMs_AcCoA+KMs_AcCoA*GOX+KMs_GOX*AcCoA+GOX*AcCoA);
vE_AceKki      = AceK*kAceKki_cat*Icdh/KAceK_Icdh*power(1+Icdh/KAceK_Icdh,nAceK-1)...
                /(power(1+Icdh/KAceK_Icdh,nAceK)+LAceK*power(1+OAA/KAceK_OAA,nAceK));
vE_AceKph      = AceK*kAceKph_cat*IcdhP/KAceK_IcdhP*power(1+IcdhP/KAceK_IcdhP,nAceKph-1)...
                /(power(1+IcdhP/KAceK_IcdhP,nAceKph)+LAceKph/power(1+OAA/KAceKph_OAA,nAceKph));
vE_G6pdh       = vG6pdh_max*G6P*NADP/((G6P+KG6pdh_G6P)*(1+NADPH/KG6pdh_NADPH_g6pinh)*(KG6pdh_NADP*(1+NADPH/KG6pdh_NADPH_nadpinh)+NADP));
vE_Pgl         = vPgl_max*(sixPGL-sixPG/KPgl_eq)/((1+H/KPgl_h1+KPgl_h2/H)*(KPgl_6PGL_m+sixPGL+KPgl_6PGL_m/KPgl_6PG_m*sixPG));
 QEdd_pH       = 1+2*power(10,pH_Edd_m-pK_Edd)/(1+power(10,pH-pK_Edd)+power(10,2*pH_Edd_m-pH-pK_Edd));
vE_Edd         = vEdd_max*QEdd_pH*(sixPG-KDPG/KEdd_eq)/(KEdd_6PG_m+sixPG+KEdd_6PG_m*KDPG/KEdd_KDPG_m);
 QEda_pH       = 1+2*power(10,pH_Eda_m-pK_Eda)/(1+power(10,pH-pK_Eda)+power(10,2*pH_Eda_m-pH-pK_Eda));
vE_Eda         = vEda_max*QEda_pH*(KDPG-GAP*PYR/KEda_eq)/(KEda_KDPG_m+KDPG+KEda_KDPG_m*(PYR/KEda_PYR_m+GAP/KEda_GAP_m+PYR*GAP/(KEda_PYR_m*KEda_GAP_m)));
vE_6Pgdh       = v6Pgdh_max*sixPG*NADP/((sixPG+K6Pgdh_6PG)*(NADP+K6Pgdh_NADP*(1+NADPH/K6Pgdh_NADPH_inh)*(1+ATP/K6Pgdh_ATP_inh)));
vE_R5pi        = vR5pi_max*(RU5P-R5P/KR5pi_eq);
vE_Rpe        = vRpe_max*(RU5P-X5P/KRpe_eq);
vE_Tkt1        = vTkt1_max*(R5P*X5P-S7P*GAP/KTkt1_eq);
vE_Tkt2        = vTkt2_max*(X5P*E4P-F6P*GAP/KTkt2_eq);
vE_Tal         = vTal_max*(GAP*S7P-E4P*F6P/KTal_eq);
vE_Cya         = vCya_max*EIIAP/(EIIAP+KCya_EIIAP);
vE_cAMPdegr    = vcAMPdegr_max*cAMP/(cAMP+KcAMPdegr_cAMP);

%********** MODEL VARIABLES (v-dependent)
OP_NADH  = (vE_Gapdh+vE_Pdh+vE_akgdh+vE_Mdh)*POratio;
OP_FADH2 = vE_Sdh*POratio_prime;
vATP     = OP_NADH + OP_FADH2 - vE_Glk - vE_Pfk + vE_Gapdh + vE_Pyk - vE_Pps + vE_Ack - vE_Acs + vE_akgdh - vE_Pck - vE_AceKki - vE_Cya;
mu       = kATP*vATP;

%********** MODEL REACTIONS (mu-dependent)
vgrowth    = mu*X;
vG_glk     = mu*((1-Cra/(Cra+Kglk_Cra))*vglk_Cra_unbound+Cra/(Cra+Kglk_Cra)*vglk_Cra_bound);
vG_pfkA    = mu*((1-Cra/(Cra+KpfkA_Cra))*vpfkA_Cra_unbound+Cra/(Cra+KpfkA_Cra)*vpfkA_Cra_bound);
vG_fbp     = mu*((1-Cra/(Cra+Kfbp_Cra))*vfbp_Cra_unbound+Cra/(Cra+Kfbp_Cra)*vfbp_Cra_bound);
vG_fbaA    = mu*((1-Cra/(Cra+KfbaA_Cra))*vfbaA_Cra_unbound+Cra/(Cra+KfbaA_Cra)*vfbaA_Cra_bound+(1-CrpcAMP/(CrpcAMP+KfbaA_Crp))*vfbaA_Crp_unbound+CrpcAMP/(CrpcAMP+KfbaA_Crp)*vfbaA_Crp_bound);
vG_gapA    = mu*((1-Cra/(Cra+KgapA_Cra))*vgapA_Cra_unbound+Cra/(Cra+KgapA_Cra)*vgapA_Cra_bound+(1-CrpcAMP/(CrpcAMP+KgapA_Crp))*vgapA_Crp_unbound+CrpcAMP/(CrpcAMP+KgapA_Crp)*vgapA_Crp_bound);
vG_pykF    = mu*((1-Cra/(Cra+KpykF_Cra))*vpykF_Cra_unbound+Cra/(Cra+KpykF_Cra)*vpykF_Cra_bound);
vG_ppsA    = mu*((1-Cra/(Cra+KppsA_Cra))*vppsA_Cra_unbound+Cra/(Cra+KppsA_Cra)*vppsA_Cra_bound);
vG_lpd     = mu*((1-PdhR/(PdhR+Klpd_PdhR))*vlpd_PdhR_unbound+PdhR/(PdhR+Klpd_PdhR)*vlpd_PdhR_bound);
vG_acs     = mu*((1-power(CrpcAMP,nacs)/(power(CrpcAMP,nacs)+power(Kacs_Crp,nacs)))*vacs_Crp_unbound+power(CrpcAMP,nacs)/(power(CrpcAMP,nacs)+power(Kacs_Crp,nacs))*vacs_Crp_bound);
vG_gltA    = mu*((1-power(CrpcAMP,ngltA)/(power(CrpcAMP,ngltA)+power(KgltA_Crp,ngltA)))*vgltA_Crp_unbound+power(CrpcAMP,ngltA)/(power(CrpcAMP,ngltA)+power(KgltA_Crp,ngltA))*vgltA_Crp_bound);
vG_icdA    = mu*((1-Cra/(Cra+KicdA_Cra))*vicdA_Cra_unbound+Cra/(Cra+KicdA_Cra)*vicdA_Cra_bound);
vG_sucAB   = mu*((1-power(CrpcAMP,nsucAB)/(power(CrpcAMP,nsucAB)+power(KsucAB_Crp,nsucAB)))*vsucAB_Crp_unbound+power(CrpcAMP,nsucAB)/(power(CrpcAMP,nsucAB)+power(KsucAB_Crp,nsucAB))*vsucAB_Crp_bound);
vG_sdhCDAB = mu*((1-power(CrpcAMP,nsdhCDAB)/(power(CrpcAMP,nsdhCDAB)+power(KsdhCDAB_Crp,nsdhCDAB)))*vsdhCDAB_Crp_unbound+power(CrpcAMP,nsdhCDAB)/(power(CrpcAMP,nsdhCDAB)+power(KsdhCDAB_Crp,nsdhCDAB))*vsdhCDAB_Crp_bound);
vG_fumABC  = mu*((1-power(CrpcAMP,nfumABC)/(power(CrpcAMP,nfumABC)+power(KfumABC_Crp,nfumABC)))*vfumABC_Crp_unbound+power(CrpcAMP,nfumABC)/(power(CrpcAMP,nfumABC)+power(KfumABC_Crp,nfumABC))*vfumABC_Crp_bound);
vG_mdh     = mu*((1-CrpcAMP/(CrpcAMP+Kmdh_Crp))*vmdh_Crp_unbound+CrpcAMP/(CrpcAMP+Kmdh_Crp)*vmdh_Crp_bound);
vG_maeB    = mu*SS_MaeB;
vG_pckA    = mu*((1-Cra/(Cra+KpckA_Cra))*vpckA_Cra_unbound+Cra/(Cra+KpckA_Cra)*vpckA_Cra_bound);
vG_ppc     = mu*SS_Ppc;
vG_aceA    = mu*power(Kmu,n_mu)/(power(mu,n_mu) + power(Kmu,n_mu)) *((1-Cra/(Cra+KaceBAK_Cra))*vaceBAK_Cra_unbound + Cra/(Cra+KaceBAK_Cra)*vaceBAK_Cra_bound + (1-CrpcAMP/(CrpcAMP+KaceBAK_Crp))*vaceBAK_Crp_unbound + CrpcAMP/(CrpcAMP+KaceBAK_Crp)*vaceBAK_Crp_bound + ( 1-aceBAK_DNA/KaceBAK_DNA*(1+PYR/KaceBAK_PYRprime)/(1+1/LaceBAK*(GOX/KaceBAK_GOX)*(1+GOX/KaceBAK_GOX)+aceBAK_DNA/KaceBAK_DNA+PYR/KaceBAK_PYR+aceBAK_DNA/KaceBAK_DNA*PYR/KaceBAK_PYRprime))*kaceBAK_cat_IclR*IclRtotal);
vG_aceB    = Factor_aceB*vG_aceA;
vG_aceK    = Factor_aceK*vG_aceA;
vD_X       = D*X;
vD_GLCfeed = D*GLCfeed;
vD_GLCex   = D*GLCex;
vD_GLC     = mu*GLC;
vD_G6P     = mu*G6P;
vD_F6P     = mu*F6P;
vD_FBP     = mu*FBP;
vD_GAP     = mu*GAP;
vD_PEP     = mu*PEP;
vD_PYR     = mu*PYR;
vD_AcCoA   = mu*AcCoA;
vD_AcP     = mu*AcP;
vD_ACEex   = D*ACEex;
vD_ICIT    = mu*ICIT;
vD_aKG     = mu*aKG;
vD_SUC     = mu*SUC;
vD_FUM     = mu*FUM;
vD_MAL     = mu*MAL;
vD_OAA     = mu*OAA;
vD_GOX     = mu*GOX;
vD_6PGL    = mu*sixPGL;
vD_6PG     = mu*sixPG;
vD_KDPG    = mu*KDPG;
vD_RU5P    = mu*RU5P;
vD_R5P     = mu*R5P;
vD_X5P     = mu*X5P;
vD_S7P     = mu*S7P;
vD_E4P     = mu*E4P;
vD_cAMP    = mu*cAMP;
vD_Glk     = (mu+kdegr)*Glk;
vD_Pfk     = (mu+kdegr)*Pfk;
vD_Fbp     = (mu+kdegr)*Fbp;
vD_Fba     = (mu+kdegr)*Fba;
vD_Gapdh   = (mu+kdegr)*Gapdh;
vD_Pyk     = (mu+kdegr)*Pyk;
vD_Pps     = (mu+kdegr)*Pps;
vD_Pdh     = (mu+kdegr)*Pdh;
vD_Acs     = (mu+kdegr)*Acs;
vD_Cs      = (mu+kdegr)*Cs;
vD_Icdh    = (mu+kdegr)*Icdh;
vD_IcdhP   = (mu+kdegr)*IcdhP;
vD_akgdh   = (mu+kdegr)*akgdh;
vD_Sdh     = (mu+kdegr)*Sdh;
vD_Fum     = (mu+kdegr)*Fum;
vD_Mdh     = (mu+kdegr)*Mdh;
vD_MaeB    = (mu+kdegr)*MaeB;
vD_Pck     = (mu+kdegr)*Pck;
vD_Ppc     = (mu+kdegr)*Ppc;
vD_Icl     = (mu+kdegr)*Icl;
vD_Ms      = (mu+kdegr)*Ms;
vD_AceK    = (mu+kdegr)*AceK;

vBM_G6P    = kBM_GLC_G6P*G6P;
vBM_F6P    = kBM_GLC_F6P*F6P;
vBM_GAP    = kBM_GLC_GAP*GAP;
vBM_PEP    = kBM_GLC_PEP*PEP;
vBM_PYR    = kBM_GLC_PYR*PYR;
vBM_AcCoA  = kBM_GLC_AcCoA*AcCoA;
vBM_aKG    = kBM_GLC_aKG*aKG;
vBM_SUC    = kBM_GLC_SUC*SUC;
vBM_FUM    = kBM_GLC_FUM*FUM;
vBM_OAA    = kBM_GLC_OAA*OAA;
vBM_R5P    = kBM_GLC_R5P*R5P;
vBM_E4P    = kBM_GLC_E4P*E4P;

Flux = [ vgrowth vPts1 vPts4 vPts4_medium vNonpts ...     %   5
    vNonpts_medium vE_Glk vE_Pgi vE_Pfk vE_Fbp ...        %  10
    vE_Fba vE_Gapdh vE_Pyk vE_Pps vE_Pdh ...              %  15
    vE_Pta vE_Ack vE_Ack_medium vE_Acs vE_Acs_medium ...  %  20
    vE_Cs vE_Icdh vE_akgdh vE_Sdh vE_Fum ...              %  25
    vE_Mdh vE_MaeB vE_Pck vE_Ppc vE_Icl ...                %  30
    vE_Ms vE_AceKki vE_AceKph vE_G6pdh vE_Pgl ...         %  35
    vE_Edd vE_Eda vE_6Pgdh vE_R5pi vE_Rpe ...            %  40
    vE_Tkt1 vE_Tkt2 vE_Tal vE_Cya vE_cAMPdegr ...         %  45
    vG_glk vG_pfkA vG_fbp vG_fbaA vG_gapA ...             %  50
    vG_pykF vG_ppsA vG_lpd vG_acs vG_gltA ...             %  55
    vG_icdA vG_sucAB vG_sdhCDAB vG_fumABC vG_mdh ...      %  60
    vG_maeB vG_pckA vG_ppc vG_aceA vG_aceB ...            %  65
    vG_aceK vD_X vD_GLCfeed vD_GLCex vD_GLC ...           %  70
    vD_G6P vD_F6P vD_FBP vD_GAP vD_PEP ...                %  75
    vD_PYR vD_AcCoA vD_AcP vD_ACEex vD_ICIT ...           %  80
    vD_aKG vD_SUC vD_FUM vD_MAL vD_OAA ...                %  85
    vD_GOX vD_6PGL vD_6PG vD_KDPG vD_RU5P ...             %  90
    vD_R5P vD_X5P vD_S7P vD_E4P vD_cAMP ...               %  95
    vD_Glk vD_Pfk vD_Fbp vD_Fba vD_Gapdh ...              % 100
    vD_Pyk vD_Pps vD_Pdh vD_Acs vD_Cs ...                 % 105
    vD_Icdh vD_IcdhP vD_akgdh vD_Sdh vD_Fum ...           % 110
    vD_Mdh vD_MaeB vD_Pck vD_Ppc vD_Icl ...                % 115
    vD_Ms vD_AceK vBM_G6P vBM_F6P vBM_GAP ...             % 120
    vBM_PEP vBM_PYR vBM_AcCoA vBM_aKG vBM_SUC ...         % 125
    vBM_FUM vBM_OAA vBM_R5P vBM_E4P alphaGLC ...          % 130
    alphaACE H EIIA CrpcAMP Crp ...                       % 135
    CraFBP Cra PdhRPYR PdhR OP_NADH ...                   % 140
    OP_FADH2 vATP mu A B ...                              % 145
    vSdh1_max vSdh2_max vFum1_max vFum2_max vMdh1_max ... % 150
    vMdh2_max QEdd_pH QEda_pH SS_MaeB SS_Ppc ];            % 155

return
