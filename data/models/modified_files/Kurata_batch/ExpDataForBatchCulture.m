function [ T_exp, X_exp, GLCex_exp, ACEex_exp ] = ExpDataForBatchCulture()

global strain_no

switch strain_no 
    % Wild Type (Toya 2010)
    case 1
        % T: h, X: g/L, GLCex: g/L, ACEex: g/L
        T_exp( 1,1) =  0.0; X_exp( 1,1) = 0.012; GLCex_exp( 1,1) = 4.042; ACEex_exp( 1,1) = 0.000;
        T_exp( 2,1) =  1.0; X_exp( 2,1) = 0.015; GLCex_exp( 2,1) = 4.033; ACEex_exp( 2,1) = 0.000;
        T_exp( 3,1) =  2.0; X_exp( 3,1) = 0.033; GLCex_exp( 3,1) = 3.964; ACEex_exp( 3,1) = 0.020;
        T_exp( 4,1) =  3.0; X_exp( 4,1) = 0.067; GLCex_exp( 4,1) = 3.843; ACEex_exp( 4,1) = 0.046;
        T_exp( 5,1) =  4.0; X_exp( 5,1) = 0.128; GLCex_exp( 5,1) = 3.696; ACEex_exp( 5,1) = 0.070;
        T_exp( 6,1) =  5.0; X_exp( 6,1) = 0.303; GLCex_exp( 6,1) = 3.351; ACEex_exp( 6,1) = 0.136;
        T_exp( 7,1) =  6.0; X_exp( 7,1) = 0.642; GLCex_exp( 7,1) = 2.440; ACEex_exp( 7,1) = 0.230;
        T_exp( 8,1) =  7.0; X_exp( 8,1) = 1.299; GLCex_exp( 8,1) = 0.622; ACEex_exp( 8,1) = 0.344;
        T_exp( 9,1) =  8.0; X_exp( 9,1) = 1.563; GLCex_exp( 9,1) = 0.013; ACEex_exp( 9,1) = 0.312;
        T_exp(10,1) =  8.5; X_exp(10,1) = 1.548; GLCex_exp(10,1) = 0.002; ACEex_exp(10,1) = 0.187;
        T_exp(11,1) =  9.0; X_exp(11,1) = 1.650; GLCex_exp(11,1) = 0.001; ACEex_exp(11,1) = 0.003;
        T_exp(12,1) = 10.0; X_exp(12,1) = 1.680; GLCex_exp(12,1) = 0.003; ACEex_exp(12,1) = 0.003;
        
    % Delta pykA/pykF (Toya 2010)
    case 26
        % T: h, X: g/L, GLCex: g/L, ACEex: g/L
        T_exp( 1,1) = 0.0; X_exp( 1,1) = 0.011; GLCex_exp( 1,1) = 3.86; ACEex_exp( 1,1) = 0.00;
        T_exp( 2,1) = 1.0; X_exp( 2,1) = 0.018; GLCex_exp( 2,1) = 3.83; ACEex_exp( 2,1) = 0.00;
        T_exp( 3,1) = 2.0; X_exp( 3,1) = 0.039; GLCex_exp( 3,1) = 3.93; ACEex_exp( 3,1) = 0.00;
        T_exp( 4,1) = 3.0; X_exp( 4,1) = 0.068; GLCex_exp( 4,1) = 3.67; ACEex_exp( 4,1) = 0.00;
        T_exp( 5,1) = 4.0; X_exp( 5,1) = 0.135; GLCex_exp( 5,1) = 3.63; ACEex_exp( 5,1) = 0.04;
        T_exp( 6,1) = 5.0; X_exp( 6,1) = 0.265; GLCex_exp( 6,1) = 3.13; ACEex_exp( 6,1) = 0.08;
        T_exp( 7,1) = 6.0; X_exp( 7,1) = 0.537; GLCex_exp( 7,1) = 2.42; ACEex_exp( 7,1) = 0.16;
        T_exp( 8,1) = 7.0; X_exp( 8,1) = 0.987; GLCex_exp( 8,1) = 1.11; ACEex_exp( 8,1) = 0.28;
        T_exp( 9,1) = 8.0; X_exp( 9,1) = 1.500; GLCex_exp( 9,1) = 0.00; ACEex_exp( 9,1) = 0.38;
        T_exp(10,1) = 9.0; X_exp(10,1) = 1.527; GLCex_exp(10,1) = 0.00; ACEex_exp(10,1) = 0.30;
        
    % Delta pgi (Toya 2010)
    case 4
        % T: h, X: g/L, GLCex: g/L, ACEex: g/L
        T_exp( 1,1) =  0; X_exp( 1,1) = 0.011; GLCex_exp( 1,1) = 4.08; ACEex_exp( 1,1) = 0.00;
        T_exp( 2,1) =  2; X_exp( 2,1) = 0.018; GLCex_exp( 2,1) = 4.09; ACEex_exp( 2,1) = 0.00;
        T_exp( 3,1) =  4; X_exp( 3,1) = 0.025; GLCex_exp( 3,1) = 4.16; ACEex_exp( 3,1) = 0.00;
        T_exp( 4,1) =  6; X_exp( 4,1) = 0.037; GLCex_exp( 4,1) = 4.14; ACEex_exp( 4,1) = 0.00;
        T_exp( 5,1) =  8; X_exp( 5,1) = 0.057; GLCex_exp( 5,1) = 4.16; ACEex_exp( 5,1) = 0.00;
        T_exp( 6,1) = 10; X_exp( 6,1) = 0.087; GLCex_exp( 6,1) = 4.04; ACEex_exp( 6,1) = 0.00;
        T_exp( 7,1) = 12; X_exp( 7,1) = 0.128; GLCex_exp( 7,1) = 4.02; ACEex_exp( 7,1) = 0.00;
        T_exp( 8,1) = 14; X_exp( 8,1) = 0.193; GLCex_exp( 8,1) = 3.79; ACEex_exp( 8,1) = 0.00;
        T_exp( 9,1) = 15; X_exp( 9,1) = 0.233; GLCex_exp( 9,1) = 3.70; ACEex_exp( 9,1) = 0.00;
        T_exp(10,1) = 16; X_exp(10,1) = 0.292; GLCex_exp(10,1) = 3.58; ACEex_exp(10,1) = 0.00;
        T_exp(11,1) = 17; X_exp(11,1) = 0.360; GLCex_exp(11,1) = 3.45; ACEex_exp(11,1) = 0.00;
        T_exp(12,1) = 18; X_exp(12,1) = 0.440; GLCex_exp(12,1) = 3.16; ACEex_exp(12,1) = 0.00;
        T_exp(13,1) = 19; X_exp(13,1) = 0.545; GLCex_exp(13,1) = 2.87; ACEex_exp(13,1) = 0.00;
        T_exp(14,1) = 20; X_exp(14,1) = 0.669; GLCex_exp(14,1) = 2.45; ACEex_exp(14,1) = 0.00;
        T_exp(15,1) = 21; X_exp(15,1) = 0.861; GLCex_exp(15,1) = 2.02; ACEex_exp(15,1) = 0.00;
        T_exp(16,1) = 22; X_exp(16,1) = 1.062; GLCex_exp(16,1) = 1.47; ACEex_exp(16,1) = 0.00;
        T_exp(17,1) = 23; X_exp(17,1) = 1.284; GLCex_exp(17,1) = 0.74; ACEex_exp(17,1) = 0.00;
        T_exp(18,1) = 25; X_exp(18,1) = 1.527; GLCex_exp(18,1) = 0.01; ACEex_exp(18,1) = 0.00;
        
    % Delta ppc (Kadir 2010)
    case 25
        % T: h, X: g/L, GLCex: g/L, ACEex: g/L
        T_exp( 1,1) =  0;     X_exp( 1,1) = 0.0358; GLCex_exp( 1,1) = 4.74;  ACEex_exp( 1,1) = 0;
        T_exp( 2,1) =  6;     X_exp( 2,1) = 0.268;  GLCex_exp( 2,1) = 4.71;  ACEex_exp( 2,1) = 0.0039;
        T_exp( 3,1) =  6.5;   X_exp( 3,1) = 0.3456; GLCex_exp( 3,1) = 4.7;   ACEex_exp( 3,1) = 0.0041;
        T_exp( 4,1) =  7;     X_exp( 4,1) = 0.44;   GLCex_exp( 4,1) = 4.62;  ACEex_exp( 4,1) = 0.0056;
        T_exp( 5,1) =  7.5;   X_exp( 5,1) = 0.6;    GLCex_exp( 5,1) = 4.49;  ACEex_exp( 5,1) = 0.0068;
        T_exp( 6,1) =  8;     X_exp( 6,1) = 0.765;  GLCex_exp( 6,1) = 4.34;  ACEex_exp( 6,1) = 0.0062;
        T_exp( 7,1) =  8.33;  X_exp( 7,1) = 0.871;  GLCex_exp( 7,1) = 4.27;  ACEex_exp( 7,1) = 0.0057;
        T_exp( 8,1) =  8.67;  X_exp( 8,1) = 1.01;   GLCex_exp( 8,1) = 4.16;  ACEex_exp( 8,1) = 0.005;
        T_exp( 9,1) =  8.88;  X_exp( 9,1) = 1.1;    GLCex_exp( 9,1) = 4.09;  ACEex_exp( 9,1) = 0.004;
        T_exp(10,1) =  9.33;  X_exp(10,1) = 1.33;   GLCex_exp(10,1) = 3.87;  ACEex_exp(10,1) = 0.0051;
        T_exp(11,1) =  9.67;  X_exp(11,1) = 1.508;  GLCex_exp(11,1) = 3.59;  ACEex_exp(11,1) = 0.0051;
        T_exp(12,1) = 10.32;  X_exp(12,1) = 1.856;  GLCex_exp(12,1) = 2.97;  ACEex_exp(12,1) = 0.003;
        T_exp(13,1) = 10.67;  X_exp(13,1) = 2.01;   GLCex_exp(13,1) = 2.6;   ACEex_exp(13,1) = 0.0028;
        T_exp(14,1) = 11;     X_exp(14,1) = 2.11;   GLCex_exp(14,1) = 2.31;  ACEex_exp(14,1) = 0.0062;
        T_exp(15,1) = 11.15;  X_exp(15,1) = 2.16;   GLCex_exp(15,1) = 2.02;  ACEex_exp(15,1) = 0.005;
        T_exp(16,1) = 11.33;  X_exp(16,1) = 2.23;   GLCex_exp(16,1) = 1.59;  ACEex_exp(16,1) = 0.004;
        T_exp(17,1) = 11.5;   X_exp(17,1) = 2.26;   GLCex_exp(17,1) = 1.275; ACEex_exp(17,1) = 0.003;
        T_exp(18,1) = 11.66;  X_exp(18,1) = 2.3;    GLCex_exp(18,1) = 1.055; ACEex_exp(18,1) = 0.0039;
        T_exp(19,1) = 12;     X_exp(19,1) = 2.37;   GLCex_exp(19,1) = 0.34;  ACEex_exp(19,1) = 0.003;
        T_exp(20,1) = 12.133; X_exp(20,1) = 2.406;  GLCex_exp(20,1) = 0;     ACEex_exp(20,1) = 0.0025;
        T_exp(21,1) = 12.275; X_exp(21,1) = 2.4058; GLCex_exp(21,1) = 0;     ACEex_exp(21,1) = 0.0032;
        T_exp(22,1) = 12.67;  X_exp(22,1) = 2.4;    GLCex_exp(22,1) = 0;     ACEex_exp(22,1) = 0.0027;
        T_exp(23,1) = 12.883; X_exp(23,1) = 2.399;  GLCex_exp(23,1) = 0;     ACEex_exp(23,1) = 0.0012;

    otherwise
        fprintf('Unexpected Strain Name!\n');
end

% Unit conversion: g/L -> mM (mmol/Lmedium)
GLCex_exp = GLCex_exp / 180 * 1e+3;
ACEex_exp = ACEex_exp /  60 * 1e+3;

return

 	 	 		 	