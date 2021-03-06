function [Yidx, Fidx] = setIndex()

% Index for variable;
Yidx.X      = 1;
Yidx.GLCex  = 2;
Yidx.GLC    = 3;
Yidx.G6P    = 4;
Yidx.F6P    = 5;
Yidx.FBP    = 6;
Yidx.GAP    = 7;
Yidx.PEP    = 8;
Yidx.PYR    = 9;
Yidx.AcCoA  = 10;
Yidx.AcP    = 11;
Yidx.ACEex  = 12;
Yidx.ICIT   = 13;
Yidx.aKG    = 14;
Yidx.SUC    = 15;
Yidx.FUM    = 16;
Yidx.MAL    = 17;
Yidx.OAA    = 18;
Yidx.GOX    = 19;
Yidx.sixPGL = 20;
Yidx.sixPG  = 21;
Yidx.KDPG   = 22;
Yidx.RU5P   = 23;
Yidx.R5P    = 24;
Yidx.X5P    = 25;
Yidx.S7P    = 26;
Yidx.E4P    = 27;
Yidx.cAMP   = 28;
Yidx.EIIAP  = 29;
Yidx.Glk    = 30;
Yidx.Pfk    = 31;
Yidx.Fbp    = 32;
Yidx.Fba    = 33;
Yidx.Gapdh  = 34;
Yidx.Pyk    = 35;
Yidx.Pps    = 36;
Yidx.Pdh    = 37;
Yidx.Acs    = 38;
Yidx.Cs     = 39;
Yidx.Icdh   = 40;
Yidx.IcdhP  = 41;
Yidx.akgdh  = 42;
Yidx.Sdh    = 43;
Yidx.Fum    = 44;
Yidx.Mdh    = 45;
Yidx.MaeB   = 46;
Yidx.Pck    = 47;
Yidx.Ppc    = 48;
Yidx.Icl    = 49;
Yidx.Ms     = 50;
Yidx.AceK   = 51;

% Index for flux
Fidx.vgrowth        = 1;
Fidx.vPts1          = 2;
Fidx.vPts4          = 3;
Fidx.vPts4_medium   = 4;
Fidx.vNonpts        = 5;
Fidx.vNonpts_medium = 6;
Fidx.vE_Glk         = 7;
Fidx.vE_Pgi         = 8;
Fidx.vE_Pfk         = 9;
Fidx.vE_Fbp         = 10;
Fidx.vE_Fba         = 11;
Fidx.vE_Gapdh       = 12;
Fidx.vE_Pyk         = 13;
Fidx.vE_Pps         = 14;
Fidx.vE_Pdh         = 15;
Fidx.vE_Pta         = 16;
Fidx.vE_Ack         = 17;
Fidx.vE_Ack_medium  = 18;
Fidx.vE_Acs         = 19;
Fidx.vE_Acs_medium  = 20;
Fidx.vE_Cs          = 21;
Fidx.vE_Icdh        = 22;
Fidx.vE_akgdh       = 23;
Fidx.vE_Sdh         = 24;
Fidx.vE_Fum         = 25;
Fidx.vE_Mdh         = 26;
Fidx.vE_MaeB        = 27;
Fidx.vE_Pck         = 28;
Fidx.vE_Ppc         = 29;
Fidx.vE_Icl         = 30;
Fidx.vE_Ms          = 31;
Fidx.vE_AceKki      = 32;
Fidx.vE_AceKph      = 33;
Fidx.vE_G6pdh       = 34;
Fidx.vE_Pgl         = 35;
Fidx.vE_Edd         = 36;
Fidx.vE_Eda         = 37;
Fidx.vE_6Pgdh       = 38;
Fidx.vE_R5pi        = 39;
Fidx.vE_Rpe         = 40;
Fidx.vE_Tkt1        = 41;
Fidx.vE_Tkt2        = 42;
Fidx.vE_Tal         = 43;
Fidx.vE_Cya         = 44;
Fidx.vE_cAMPdegr    = 45;
Fidx.vG_glk         = 46;
Fidx.vG_pfkA        = 47;
Fidx.vG_fbp         = 48;
Fidx.vG_fbaA        = 49;
Fidx.vG_gapA        = 50;
Fidx.vG_pykF        = 51;
Fidx.vG_ppsA        = 52;
Fidx.vG_lpd         = 53;
Fidx.vG_acs         = 54;
Fidx.vG_gltA        = 55;
Fidx.vG_icdA        = 56;
Fidx.vG_sucAB       = 57;
Fidx.vG_sdhCDAB     = 58;
Fidx.vG_fumABC      = 59;
Fidx.vG_mdh         = 60;
Fidx.vG_maeB        = 61;
Fidx.vG_pckA        = 62;
Fidx.vG_ppc         = 63;
Fidx.vG_aceA        = 64;
Fidx.vG_aceB        = 65;
Fidx.vG_aceK        = 66;
Fidx.vD_X           = 67;
Fidx.vD_GLCfeed     = 68;
Fidx.vD_GLCex       = 69;
Fidx.vD_GLC         = 70;
Fidx.vD_G6P         = 71;
Fidx.vD_F6P         = 72;
Fidx.vD_FBP         = 73;
Fidx.vD_GAP         = 74;
Fidx.vD_PEP         = 75;
Fidx.vD_PYR         = 76;
Fidx.vD_AcCoA       = 77;
Fidx.vD_AcP         = 78;
Fidx.vD_ACEex       = 79;
Fidx.vD_ICIT        = 80;
Fidx.vD_aKG         = 81;
Fidx.vD_SUC         = 82;
Fidx.vD_FUM         = 83;
Fidx.vD_MAL         = 84;
Fidx.vD_OAA         = 85;
Fidx.vD_GOX         = 86;
Fidx.vD_6PGL        = 87;
Fidx.vD_6PG         = 88;
Fidx.vD_KDPG        = 89;
Fidx.vD_RU5P        = 90;
Fidx.vD_R5P         = 91;
Fidx.vD_X5P         = 92;
Fidx.vD_S7P         = 93;
Fidx.vD_E4P         = 94;
Fidx.vD_cAMP        = 95;
Fidx.vD_Glk         = 96;
Fidx.vD_Pfk         = 97;
Fidx.vD_Fbp         = 98;
Fidx.vD_Fba         = 99;
Fidx.vD_Gapdh       = 100;
Fidx.vD_Pyk         = 101;
Fidx.vD_Pps         = 102;
Fidx.vD_Pdh         = 103;
Fidx.vD_Acs         = 104;
Fidx.vD_Cs          = 105;
Fidx.vD_Icdh        = 106;
Fidx.vD_IcdhP       = 107;
Fidx.vD_akgdh       = 108;
Fidx.vD_Sdh         = 109;
Fidx.vD_Fum         = 110;
Fidx.vD_Mdh         = 111;
Fidx.vD_MaeB         = 112;
Fidx.vD_Pck         = 113;
Fidx.vD_Ppc         = 114;
Fidx.vD_Icl         = 115;
Fidx.vD_Ms          = 116;
Fidx.vD_AceK        = 117;
Fidx.vBM_G6P        = 118;
Fidx.vBM_F6P        = 119;
Fidx.vBM_GAP        = 120;
Fidx.vBM_PEP        = 121;
Fidx.vBM_PYR        = 122;
Fidx.vBM_AcCoA      = 123;
Fidx.vBM_aKG        = 124;
Fidx.vBM_SUC        = 125;
Fidx.vBM_FUM        = 126;
Fidx.vBM_OAA        = 127;
Fidx.vBM_R5P        = 128;
Fidx.vBM_E4P        = 129;
Fidx.alphaGLC       = 130;
Fidx.alphaACE       = 131;
Fidx.H              = 132;
Fidx.EIIA           = 133;
Fidx.CrpcAMP        = 134;
Fidx.Crp            = 135;
Fidx.CraFBP         = 136;
Fidx.Cra            = 137;
Fidx.PdhRPYR        = 138;
Fidx.PdhR           = 139;
Fidx.OP_NADH        = 140;
Fidx.OP_FADH2       = 141;
Fidx.vATP           = 142;
Fidx.mu             = 143;
Fidx.A              = 144;
Fidx.B              = 145;
Fidx.vSdh1_max      = 146;
Fidx.vSdh2_max      = 147;
Fidx.vFum1_max      = 148;
Fidx.vFum2_max      = 149;
Fidx.vMdh1_max      = 150;
Fidx.vMdh2_max      = 151;
Fidx.QEdd_pH        = 152;
Fidx.QEda_pH        = 153;
Fidx.SS_MaeB         = 154;
Fidx.SS_Ppc         = 155;
