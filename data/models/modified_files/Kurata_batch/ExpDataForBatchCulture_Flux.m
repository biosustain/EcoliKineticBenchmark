function [ T_exp, F_exp] = ExpDataForBatchCulture_Flux()
global strain_no

switch strain_no
    %{
    % Wild Type (Toya 2010, Tables S-IV and S-V) Note that After 8 h, the direction of Ack reaction is reversed.
    case 1
        T_exp(1)    =   5; T_exp(2)    =   6; T_exp(3)    =   7; T_exp(4)    =      8; T_exp(5)    =    8.5; T_exp(6)    =      9;
        F_exp(1, 1) = 100; F_exp(2, 1) = 100; F_exp(3, 1) = 100; F_exp(4, 1) = -1e+10; F_exp(5, 1) = -1e+10; F_exp(6, 1) = -1e+10;
        F_exp(1, 2) =  55; F_exp(2, 2) =  60; F_exp(3, 2) =  66; F_exp(4, 2) = -1e+10; F_exp(5, 2) = -1e+10; F_exp(6, 2) = -1e+10;
        F_exp(1, 3) =  77; F_exp(2, 3) =  80; F_exp(3, 3) =  82; F_exp(4, 3) = -1e+10; F_exp(5, 3) = -1e+10; F_exp(6, 3) = -1e+10;
        F_exp(1, 4) =  77; F_exp(2, 4) =  80; F_exp(3, 4) =  82; F_exp(4, 4) = -1e+10; F_exp(5, 4) = -1e+10; F_exp(6, 4) = -1e+10;
        F_exp(1, 5) =  77; F_exp(2, 5) =  80; F_exp(3, 5) =  82; F_exp(4, 5) = -1e+10; F_exp(5, 5) = -1e+10; F_exp(6, 5) = -1e+10;
        F_exp(1, 6) = 163; F_exp(2, 6) = 167; F_exp(3, 6) = 169; F_exp(4, 6) = -1e+10; F_exp(5, 6) = -1e+10; F_exp(6, 6) = -1e+10;
        F_exp(1, 7) = 152; F_exp(2, 7) = 157; F_exp(3, 7) = 160; F_exp(4, 7) = -1e+10; F_exp(5, 7) = -1e+10; F_exp(6, 7) = -1e+10;
        F_exp(1, 8) =  25; F_exp(2, 8) =  18; F_exp(3, 8) =  13; F_exp(4, 8) = -1e+10; F_exp(5, 8) = -1e+10; F_exp(6, 8) = -1e+10;
        F_exp(1, 9) = 103; F_exp(2, 9) = 112; F_exp(3, 9) = 116; F_exp(4, 9) = -1e+10; F_exp(5, 9) = -1e+10; F_exp(6, 9) = -1e+10;
        F_exp(1,10) =  37; F_exp(2,10) =  23; F_exp(3,10) =  20; F_exp(4,10) =    100; F_exp(5,10) =    100; F_exp(6,10) =    100;
        F_exp(1,11) =  44; F_exp(2,11) =  39; F_exp(3,11) =  33; F_exp(4,11) = -1e+10; F_exp(5,11) = -1e+10; F_exp(6,11) = -1e+10;
        F_exp(1,12) =  44; F_exp(2,12) =  39; F_exp(3,12) =  33; F_exp(4,12) = -1e+10; F_exp(5,12) = -1e+10; F_exp(6,12) = -1e+10;
        F_exp(1,13) =  23; F_exp(2,13) =  20; F_exp(3,13) =  17; F_exp(4,13) = -1e+10; F_exp(5,13) = -1e+10; F_exp(6,13) = -1e+10;
        F_exp(1,14) =  20; F_exp(2,14) =  18; F_exp(3,14) =  16; F_exp(4,14) = -1e+10; F_exp(5,14) = -1e+10; F_exp(6,14) = -1e+10;
        F_exp(1,15) =  13; F_exp(2,15) =  12; F_exp(3,15) =  10; F_exp(4,15) = -1e+10; F_exp(5,15) = -1e+10; F_exp(6,15) = -1e+10;
        F_exp(1,16) =  13; F_exp(2,16) =  12; F_exp(3,16) =  10; F_exp(4,16) = -1e+10; F_exp(5,16) = -1e+10; F_exp(6,16) = -1e+10;
        F_exp(1,17) =  10; F_exp(2,17) =   9; F_exp(3,17) =   7; F_exp(4,17) = -1e+10; F_exp(5,17) = -1e+10; F_exp(6,17) = -1e+10;
        F_exp(1,18) =  36; F_exp(2,18) =  62; F_exp(3,18) =  69; F_exp(4,18) =     45; F_exp(5,18) =     56; F_exp(6,18) =     64;
        F_exp(1,19) =  36; F_exp(2,19) =  62; F_exp(3,19) =  69; F_exp(4,19) =      4; F_exp(5,19) =     18; F_exp(6,19) =     34;
        F_exp(1,20) =  27; F_exp(2,20) =  54; F_exp(3,20) =  61; F_exp(4,20) =      0; F_exp(5,20) =     16; F_exp(6,20) =     32;
        F_exp(1,21) =  27; F_exp(2,21) =  54; F_exp(3,21) =  61; F_exp(4,21) =     41; F_exp(5,21) =     54; F_exp(6,21) =     63;
        F_exp(1,22) =  27; F_exp(2,22) =  54; F_exp(3,22) =  61; F_exp(4,22) =     41; F_exp(5,22) =     54; F_exp(6,22) =     63;
        F_exp(1,23) =  27; F_exp(2,23) =  40; F_exp(3,23) =  39; F_exp(4,23) =     82; F_exp(5,23) =     91; F_exp(6,23) =     93;
        F_exp(1,24) =  23; F_exp(2,24) =  36; F_exp(3,24) =  43; F_exp(4,24) =     31; F_exp(5,24) =     32; F_exp(6,24) =     26;
        F_exp(1,25) =   0; F_exp(2,25) =  14; F_exp(3,25) =  23; F_exp(4,25) = -1e+10; F_exp(5,25) = -1e+10; F_exp(6,25) = -1e+10;
        F_exp(1,26) =   0; F_exp(2,26) =   0; F_exp(3,26) =   0; F_exp(4,26) =     31; F_exp(5,26) =     32; F_exp(6,26) =     26;
        F_exp(1,27) =   0; F_exp(2,27) =   0; F_exp(3,27) =   0; F_exp(4,27) =     29; F_exp(5,27) =     31; F_exp(6,27) =     26;
        F_exp(1,28) =   0; F_exp(2,28) =   0; F_exp(3,28) =   0; F_exp(4,28) = -1e+10; F_exp(5,28) = -1e+10; F_exp(6,28) = -1e+10;
        
    % Delta pykA/pykF (Toya 2010, Table S-VI)
    case 13
        T_exp(1)    =   5; T_exp(2)   =   6; T_exp(3)  =   7;
        F_exp(1, 1) = 100; F_exp(2, 1) = 100; F_exp(3, 1) = 100;
        F_exp(1, 2) =  57; F_exp(2, 2) =  46; F_exp(3, 2) =  51;
        F_exp(1, 3) =  79; F_exp(2, 3) =  75; F_exp(3, 3) =  77;
        F_exp(1, 4) =  79; F_exp(2, 4) =  75; F_exp(3, 4) =  77;
        F_exp(1, 5) =  79; F_exp(2, 5) =  75; F_exp(3, 5) =  77;
        F_exp(1, 6) = 168; F_exp(2, 6) = 163; F_exp(3, 6) = 165;
        F_exp(1, 7) = 158; F_exp(2, 7) = 154; F_exp(3, 7) = 156;
        F_exp(1, 8) =   0; F_exp(2, 8) =   0; F_exp(3, 8) =   0;
        F_exp(1, 9) = 117; F_exp(2, 9) = 110; F_exp(3, 9) = 113;
        F_exp(1,10) =  32; F_exp(2,10) =  30; F_exp(3,10) =  28;
        F_exp(1,11) =  42; F_exp(2,11) =  52; F_exp(3,11) =  47;
        F_exp(1,12) =  42; F_exp(2,12) =  52; F_exp(3,12) =  47;
        F_exp(1,13) =  23; F_exp(2,13) =  30; F_exp(3,13) =  26;
        F_exp(1,14) =  19; F_exp(2,14) =  23; F_exp(3,14) =  21;
        F_exp(1,15) =  13; F_exp(2,15) =  16; F_exp(3,15) =  15;
        F_exp(1,16) =  13; F_exp(2,16) =  16; F_exp(3,16) =  15;
        F_exp(1,17) =  10; F_exp(2,17) =  13; F_exp(3,17) =  12;
        F_exp(1,18) =  60; F_exp(2,18) =  54; F_exp(3,18) =  59;
        F_exp(1,19) =  60; F_exp(2,19) =  54; F_exp(3,19) =  59;
        F_exp(1,20) =  52; F_exp(2,20) =  46; F_exp(3,20) =  51;
        F_exp(1,21) =  52; F_exp(2,21) =  46; F_exp(3,21) =  51;
        F_exp(1,22) =  52; F_exp(2,22) =  46; F_exp(3,22) =  51;
        F_exp(1,23) =  17; F_exp(2,23) =  16; F_exp(3,23) =  19;
        F_exp(1,24) =  55; F_exp(2,24) =  50; F_exp(3,24) =  52;
        F_exp(1,25) =  36; F_exp(2,25) =  30; F_exp(3,25) =  32;
        F_exp(1,26) =   0; F_exp(2,26) =   0; F_exp(3,26) =   0;
        F_exp(1,27) =   0; F_exp(2,27) =   0; F_exp(3,27) =   0;
        F_exp(1,28) =   0; F_exp(2,28) =   0; F_exp(3,28) =   0;
        
    % Delta pgi (Toya 2010, Table S-VII)
    case 4
        T_exp(1)    =  16; T_exp(2)    =  21; T_exp(3)    =  23;
        F_exp(1, 1) = 100; F_exp(2, 1) = 100; F_exp(3, 1) = 100;
        F_exp(1, 2) =   0; F_exp(2, 2) =   0; F_exp(3, 2) =   0;
        F_exp(1, 3) =  54; F_exp(2, 3) =  53; F_exp(3, 3) =  54;
        F_exp(1, 4) =  54; F_exp(2, 4) =  53; F_exp(3, 4) =  54;
        F_exp(1, 5) =  54; F_exp(2, 5) =  53; F_exp(3, 5) =  54;
        F_exp(1, 6) = 139; F_exp(2, 6) = 139; F_exp(3, 6) = 143;
        F_exp(1, 7) = 126; F_exp(2, 7) = 129; F_exp(3, 7) = 134;
        F_exp(1, 8) =   8; F_exp(2, 8) =  25; F_exp(3, 8) =   0;
        F_exp(1, 9) = 103; F_exp(2, 9) = 134; F_exp(3, 9) = 117;
        F_exp(1,10) =   0; F_exp(2,10) =   0; F_exp(3,10) =   0;
        F_exp(1,11) =  98; F_exp(2,11) =  98; F_exp(3,11) =  99;
        F_exp(1,12) =  93; F_exp(2,12) =  88; F_exp(3,12) =  88;
        F_exp(1,13) =  55; F_exp(2,13) =  53; F_exp(3,13) =  54;
        F_exp(1,14) =  38; F_exp(2,14) =  35; F_exp(3,14) =  34;
        F_exp(1,15) =  29; F_exp(2,15) =  28; F_exp(3,15) =  28;
        F_exp(1,16) =  29; F_exp(2,16) =  28; F_exp(3,16) =  28;
        F_exp(1,17) =  26; F_exp(2,17) =  25; F_exp(3,17) =  26;
        F_exp(1,18) =  41; F_exp(2,18) =  64; F_exp(3,18) =  82;
        F_exp(1,19) =  13; F_exp(2,19) =  23; F_exp(3,19) =  70;
        F_exp(1,20) =   3; F_exp(2,20) =  14; F_exp(3,20) =  63;
        F_exp(1,21) =  31; F_exp(2,21) =  56; F_exp(3,21) =  75;
        F_exp(1,22) =  31; F_exp(2,22) =  56; F_exp(3,22) =  75;
        F_exp(1,23) =  43; F_exp(2,23) =  78; F_exp(3,23) =  62;
        F_exp(1,24) =  14; F_exp(2,24) =   0; F_exp(3,24) =  31;
        F_exp(1,25) =  16; F_exp(2,25) =  20; F_exp(3,25) =  25;
        F_exp(1,26) =  28; F_exp(2,26) =  42; F_exp(3,26) =  12;
        F_exp(1,27) =  28; F_exp(2,27) =  42; F_exp(3,27) =  12;
        F_exp(1,28) =   5; F_exp(2,28) =  10; F_exp(3,28) =  10;
        
    otherwise
        fprintf('Unexpected Strain Name!\n');
    %}    
    % Wild Type (Toya 2010, Tables S-IV and S-V) Note that After 8 h, the direction of Ack reaction is reversed.
    case 1
        T_exp(1)    =    5; T_exp(2)    =    6; T_exp(3)    =    7; T_exp(4)    =      8; T_exp(5)    =    8.5; T_exp(6)    =      9;
        F_exp(1, 1) = 11.7; F_exp(2, 1) = 11.1; F_exp(3, 1) =  7.7; F_exp(4, 1) = -1e+10; F_exp(5, 1) = -1e+10; F_exp(6, 1) = -1e+10;
        F_exp(1, 2) =  6.4; F_exp(2, 2) =  6.6; F_exp(3, 2) =  5.1; F_exp(4, 2) = -1e+10; F_exp(5, 2) = -1e+10; F_exp(6, 2) = -1e+10;
        F_exp(1, 3) =    9; F_exp(2, 3) =  8.8; F_exp(3, 3) =  6.3; F_exp(4, 3) = -1e+10; F_exp(5, 3) = -1e+10; F_exp(6, 3) = -1e+10;
        F_exp(1, 4) =    9; F_exp(2, 4) =  8.8; F_exp(3, 4) =  6.3; F_exp(4, 4) = -1e+10; F_exp(5, 4) = -1e+10; F_exp(6, 4) = -1e+10;
        F_exp(1, 5) =    9; F_exp(2, 5) =  8.8; F_exp(3, 5) =  6.3; F_exp(4, 5) = -1e+10; F_exp(5, 5) = -1e+10; F_exp(6, 5) = -1e+10;
        F_exp(1, 6) = 19.1; F_exp(2, 6) = 18.5; F_exp(3, 6) =   13; F_exp(4, 6) = -1e+10; F_exp(5, 6) = -1e+10; F_exp(6, 6) = -1e+10;
        F_exp(1, 7) = 17.8; F_exp(2, 7) = 17.4; F_exp(3, 7) = 12.3; F_exp(4, 7) = -1e+10; F_exp(5, 7) = -1e+10; F_exp(6, 7) = -1e+10;
        F_exp(1, 8) =  2.9; F_exp(2, 8) =    2; F_exp(3, 8) =    1; F_exp(4, 8) = -1e+10; F_exp(5, 8) = -1e+10; F_exp(6, 8) = -1e+10;
        F_exp(1, 9) =   12; F_exp(2, 9) = 12.5; F_exp(3, 9) =  8.9; F_exp(4, 9) = -1e+10; F_exp(5, 9) = -1e+10; F_exp(6, 9) = -1e+10;
        F_exp(1,10) =  4.3; F_exp(2,10) =  2.6; F_exp(3,10) =  1.5; F_exp(4,10) =    2.9; F_exp(5,10) =    3.2; F_exp(6,10) =    3.6;
        F_exp(1,11) =  5.1; F_exp(2,11) =  4.3; F_exp(3,11) =  2.5; F_exp(4,11) = -1e+10; F_exp(5,11) = -1e+10; F_exp(6,11) = -1e+10;
        F_exp(1,12) =  5.1; F_exp(2,12) =  4.3; F_exp(3,12) =  2.5; F_exp(4,12) = -1e+10; F_exp(5,12) = -1e+10; F_exp(6,12) = -1e+10;
        F_exp(1,13) =  2.7; F_exp(2,13) =  2.3; F_exp(3,13) =  1.3; F_exp(4,13) = -1e+10; F_exp(5,13) = -1e+10; F_exp(6,13) = -1e+10;
        F_exp(1,14) =  2.4; F_exp(2,14) =    2; F_exp(3,14) =  1.2; F_exp(4,14) = -1e+10; F_exp(5,14) = -1e+10; F_exp(6,14) = -1e+10;
        F_exp(1,15) =  1.6; F_exp(2,15) =  1.3; F_exp(3,15) =  0.7; F_exp(4,15) = -1e+10; F_exp(5,15) = -1e+10; F_exp(6,15) = -1e+10;
        F_exp(1,16) =  1.6; F_exp(2,16) =  1.3; F_exp(3,16) =  0.7; F_exp(4,16) = -1e+10; F_exp(5,16) = -1e+10; F_exp(6,16) = -1e+10;
        F_exp(1,17) =  1.2; F_exp(2,17) =    1; F_exp(3,17) =  0.5; F_exp(4,17) = -1e+10; F_exp(5,17) = -1e+10; F_exp(6,17) = -1e+10;
        F_exp(1,18) =  4.2; F_exp(2,18) =  6.9; F_exp(3,18) =  5.3; F_exp(4,18) =    1.3; F_exp(5,18) =    1.8; F_exp(6,18) =    2.3;
        F_exp(1,19) =  4.2; F_exp(2,19) =  6.9; F_exp(3,19) =  5.3; F_exp(4,19) =    0.1; F_exp(5,19) =    0.6; F_exp(6,19) =    1.2;
        F_exp(1,20) =  3.1; F_exp(2,20) =    6; F_exp(3,20) =  4.7; F_exp(4,20) =      0; F_exp(5,20) =    1.5; F_exp(6,20) =    1.2;
        F_exp(1,21) =  3.1; F_exp(2,21) =    6; F_exp(3,21) =  4.7; F_exp(4,21) =    1.2; F_exp(5,21) =    1.7; F_exp(6,21) =    2.3;
        F_exp(1,22) =  3.1; F_exp(2,22) =    6; F_exp(3,22) =  4.7; F_exp(4,22) =    1.2; F_exp(5,22) =    1.7; F_exp(6,22) =    2.3;
        F_exp(1,23) =  3.1; F_exp(2,23) =  4.4; F_exp(3,23) =    3; F_exp(4,23) =    2.4; F_exp(5,23) =    2.9; F_exp(6,23) =    3.4;
        F_exp(1,24) =  2.7; F_exp(2,24) =  3.9; F_exp(3,24) =  3.3; F_exp(4,24) =    0.9; F_exp(5,24) =      1; F_exp(6,24) =      1;
        F_exp(1,25) =    0; F_exp(2,25) =  1.6; F_exp(3,25) =  1.8; F_exp(4,25) = -1e+10; F_exp(5,25) = -1e+10; F_exp(6,25) = -1e+10;
        F_exp(1,26) =    0; F_exp(2,26) =    0; F_exp(3,26) =    0; F_exp(4,26) =    1.2; F_exp(5,26) =    1.2; F_exp(6,26) =    1.1;
        F_exp(1,27) =    0; F_exp(2,27) =    0; F_exp(3,27) =    0; F_exp(4,27) =    1.2; F_exp(5,27) =    1.2; F_exp(6,27) =    1.1;
        F_exp(1,28) =    0; F_exp(2,28) =    0; F_exp(3,28) =    0; F_exp(4,28) = -1e+10; F_exp(5,28) = -1e+10; F_exp(6,28) = -1e+10;
        
    % Delta pykA/pykF (Toya 2010, Table S-VI)
    case 26
        T_exp(1)    =    5; T_exp(2)   =     6; T_exp(3)  =      7;
        F_exp(1, 1) = 12.3; F_exp(2, 1) = 10.2; F_exp(3, 1) =  7.2;
        F_exp(1, 2) =    7; F_exp(2, 2) =  4.7; F_exp(3, 2) =  3.7;
        F_exp(1, 3) =  9.7; F_exp(2, 3) =  7.7; F_exp(3, 3) =  5.6;
        F_exp(1, 4) =  9.7; F_exp(2, 4) =  7.7; F_exp(3, 4) =  5.6;
        F_exp(1, 5) =  9.7; F_exp(2, 5) =  7.7; F_exp(3, 5) =  5.6;
        F_exp(1, 6) = 20.6; F_exp(2, 6) = 16.6; F_exp(3, 6) = 11.9;
        F_exp(1, 7) = 19.5; F_exp(2, 7) = 15.6; F_exp(3, 7) = 11.2;
        F_exp(1, 8) =    0; F_exp(2, 8) =    0; F_exp(3, 8) =    0;
        F_exp(1, 9) = 14.4; F_exp(2, 9) = 11.2; F_exp(3, 9) =  8.1;
        F_exp(1,10) =    4; F_exp(2,10) =  3.1; F_exp(3,10) =  2.1;
        F_exp(1,11) =  5.1; F_exp(2,11) =  5.3; F_exp(3,11) =  3.4;
        F_exp(1,12) =  5.1; F_exp(2,12) =  5.3; F_exp(3,12) =  3.4;
        F_exp(1,13) =  2.8; F_exp(2,13) =    3; F_exp(3,13) =  1.9;
        F_exp(1,14) =  2.3; F_exp(2,14) =  2.3; F_exp(3,14) =  1.5;
        F_exp(1,15) =  1.6; F_exp(2,15) =  1.6; F_exp(3,15) =  1.1;
        F_exp(1,16) =  1.6; F_exp(2,16) =  1.6; F_exp(3,16) =  1.1;
        F_exp(1,17) =  1.2; F_exp(2,17) =  1.4; F_exp(3,17) =  0.9;
        F_exp(1,18) =  7.3; F_exp(2,18) =  5.5; F_exp(3,18) =  4.2;
        F_exp(1,19) =  7.3; F_exp(2,19) =  5.5; F_exp(3,19) =  4.2;
        F_exp(1,20) =  6.4; F_exp(2,20) =  4.7; F_exp(3,20) =  3.7;
        F_exp(1,21) =  6.4; F_exp(2,21) =  4.7; F_exp(3,21) =  3.7;
        F_exp(1,22) =  6.4; F_exp(2,22) =  4.7; F_exp(3,22) =  3.7;
        F_exp(1,23) =    2; F_exp(2,23) =  1.6; F_exp(3,23) =  1.4;
        F_exp(1,24) =  6.8; F_exp(2,24) =  5.1; F_exp(3,24) =  3.8;
        F_exp(1,25) =  4.4; F_exp(2,25) =    3; F_exp(3,25) =  2.3;
        F_exp(1,26) =    0; F_exp(2,26) =    0; F_exp(3,26) =    0;
        F_exp(1,27) =    0; F_exp(2,27) =    0; F_exp(3,27) =    0;
        F_exp(1,28) =    0; F_exp(2,28) =    0; F_exp(3,28) =    0;
        
    % Delta pgi (Toya 2010, Table S-VII)
    case 4
        T_exp(1)    =  16; T_exp(2)    =  21; T_exp(3)    =  23;
        F_exp(1, 1) = 2.5; F_exp(2, 1) = 3.3; F_exp(3, 1) =   3;
        F_exp(1, 2) =   0; F_exp(2, 2) =   0; F_exp(3, 2) =   0;
        F_exp(1, 3) = 1.4; F_exp(2, 3) = 1.7; F_exp(3, 3) = 1.6;
        F_exp(1, 4) = 1.4; F_exp(2, 4) = 1.7; F_exp(3, 4) = 1.6;
        F_exp(1, 5) = 1.4; F_exp(2, 5) = 1.7; F_exp(3, 5) = 1.6;
        F_exp(1, 6) = 3.5; F_exp(2, 6) = 4.6; F_exp(3, 6) = 4.3;
        F_exp(1, 7) = 3.2; F_exp(2, 7) = 4.3; F_exp(3, 7) =   4;
        F_exp(1, 8) = 0.2; F_exp(2, 8) = 0.8; F_exp(3, 8) =   0;
        F_exp(1, 9) = 2.6; F_exp(2, 9) = 4.4; F_exp(3, 9) = 3.5;
        F_exp(1,10) =   0; F_exp(2,10) =   0; F_exp(3,10) =   0;
        F_exp(1,11) = 2.5; F_exp(2,11) = 3.2; F_exp(3,11) = 2.9;
        F_exp(1,12) = 2.4; F_exp(2,12) = 2.9; F_exp(3,12) = 2.6;
        F_exp(1,13) = 1.4; F_exp(2,13) = 1.8; F_exp(3,13) = 1.6;
        F_exp(1,14) =   1; F_exp(2,14) = 1.2; F_exp(3,14) =   1;
        F_exp(1,15) = 0.8; F_exp(2,15) = 0.9; F_exp(3,15) = 0.8;
        F_exp(1,16) = 0.8; F_exp(2,16) = 0.9; F_exp(3,16) = 0.8;
        F_exp(1,17) = 0.7; F_exp(2,17) = 0.8; F_exp(3,17) = 0.8;
        F_exp(1,18) =   1; F_exp(2,18) = 2.1; F_exp(3,18) = 2.4;
        F_exp(1,19) = 0.3; F_exp(2,19) = 0.8; F_exp(3,19) = 2.1;
        F_exp(1,20) = 0.1; F_exp(2,20) = 0.5; F_exp(3,20) = 1.9;
        F_exp(1,21) = 0.8; F_exp(2,21) = 1.8; F_exp(3,21) = 2.2;
        F_exp(1,22) = 0.8; F_exp(2,22) = 1.8; F_exp(3,22) = 2.2;
        F_exp(1,23) = 1.1; F_exp(2,23) = 2.6; F_exp(3,23) = 1.9;
        F_exp(1,24) = 0.4; F_exp(2,24) =   0; F_exp(3,24) = 0.9;
        F_exp(1,25) = 0.4; F_exp(2,25) = 0.6; F_exp(3,25) = 0.7;
        F_exp(1,26) = 0.7; F_exp(2,26) = 1.4; F_exp(3,26) = 0.4;
        F_exp(1,27) = 0.7; F_exp(2,27) = 1.4; F_exp(3,27) = 0.4;
        F_exp(1,28) = 0.1; F_exp(2,28) = 0.3; F_exp(3,28) = 0.3;
        
    otherwise
        fprintf('no experimental data\n');
        T_exp(1:1:3)=0;
        F_exp(1:1:3,1:1:28)=0;

end
return

%%
%{
    % T_exp: Time (hour)
    % F_exp(:, 1): vPts4
    % F_exp(:, 2): vE_Pgi
    % F_exp(:, 3): vE_Pfk-vE_Fbp
    % F_exp(:, 4): vE_Fba
    % F_exp(:, 5): vE_Tpi (not used)
    % F_exp(:, 6): vE_Pgk
    % F_exp(:, 7): vE_Eno
    % F_exp(:, 8): vE_Pyk-vE_Pps
    % F_exp(:, 9): vE_Pdh
    % F_exp(:,10): vE_Ack
    % F_exp(:,11): vE_Pgk
    % F_exp(:,12): vE_Eno
    % F_exp(:,13): vE_Rpe
    % F_exp(:,14): vE_R5pi
    % F_exp(:,15): vE_Tkt1
    % F_exp(:,16): vE_Tal
    % F_exp(:,17): vE_Tkt2
    % F_exp(:,18): vE_Cs
    % F_exp(:,19): vE_Icdh
    % F_exp(:,20): vE_akgdh
    % F_exp(:,21): vE_Sdh
    % F_exp(:,22): vE_Fum
    % F_exp(:,23): vE_Mdh
    % F_exp(:,24): vE_Ppc-vE_Pck
    % F_exp(:,25): vE_MaeB
    % F_exp(:,26): vE_Icl
    % F_exp(:,27): vE_Ms
    % F_exp(:,28): vE_Eda
%}
