% Special program to compute the decay width for exclusively transitions of
% p1 hybrid states to s quarkonia states

% We work in GeV!!

%Valor r0
setr0(3.964)
%massa
load("dades.mat","m_c","m_b")
setm_q(m_c) %Interesting state 1++ of charm (4140)

%Change this depending on the transition

ComputationIE = TransitionsP1toS('HQp1tS_full');      

[I_Cos_cell, I_Sin_cell, I_Si_cell, DeltaE, M] = ComputationIE(2, 2, 0.28, 0.467, 50);

file_name = 'GammaVSmass_HQcharm_1p1-2s_full_bad.txt';

% Value of pion mass 140MeV
m = 0.14;
% Value of the global constant that divides everithing in GeV. Fpi=92MeV
FracFpi = 1/(pi^3 * 0.092^4 * 0.187) * 4^2;
%0.092 és f_pi i 0.187GeV^2 is the string tension
%4^2 per compensar que Ce=eta/4 i tenim eta^2

%Obtain the I from the cell
I_Cos_R3 = I_Cos_cell{1, 1}; 
I_Sin_R2 = I_Sin_cell{1, 1};
I_Sin_R3 = I_Sin_cell{1, 2};
I_Sin_R4 = I_Sin_cell{1, 3};
I_Si_R = I_Si_cell{1, 1};
I_Si_R2 = I_Si_cell{1, 2};
I_Si_R3 = I_Si_cell{1, 3};
I_Si_R4 = I_Si_cell{1, 4};


% Filter values where DeltaE > M
valid_idx = (DeltaE > M);
M_valid = M(valid_idx);

I_Cos_R3 = I_Cos_R3(valid_idx); 
I_Sin_R2 = I_Sin_R2(valid_idx);
I_Sin_R3 = I_Sin_R3(valid_idx);
I_Sin_R4 = I_Sin_R4(valid_idx);
I_Si_R = I_Si_R(valid_idx);
I_Si_R2 = I_Si_R2(valid_idx);
I_Si_R3 = I_Si_R3(valid_idx);
I_Si_R4 = I_Si_R4(valid_idx);


num_valid = length(M_valid);
DW_values = zeros(1, num_valid);


for i = 1:num_valid
    dwMin1 = P1toSdwM1(I_Cos_R3(i), I_Sin_R2(i), I_Sin_R3(i), I_Sin_R4(i), I_Si_R(i), ...
    I_Si_R3(i), I_Si_R4(i), m, M_valid(i), DeltaE);
    
    dwMin0a = P1toSdwM0_a(I_Cos_R3(i), I_Sin_R3(i), I_Sin_R4(i), I_Si_R(i), I_Si_R2(i),...
    I_Si_R3(i), I_Si_R4(i), m, M_valid(i), DeltaE);

    dwMin0b = P1toSdwM0_b(I_Cos_R3(i), I_Sin_R3(i), I_Sin_R4(i), I_Si_R(i), I_Si_R2(i),...
    I_Si_R3(i), I_Si_R4(i), m, M_valid(i), DeltaE);

    dwMin0c = P1toSdwM0_c(I_Cos_R3(i), I_Sin_R3(i), I_Sin_R4(i), I_Si_R(i), I_Si_R2(i),...
    I_Si_R3(i), I_Si_R4(i), m, M_valid(i), DeltaE);

    dwMin0 = (dwMin0a + dwMin0b + dwMin0c);
   
   % Total decay width is the average over initial states Min and addition
   % to final Min=0. Whe have Min=+1,-1,0 and Jin=1
    DW_values(i) = ( dwMin1*2 + dwMin0 )/3;
end

% Now with each part of the gamma separated by parameters we define the
% parameters computed with mathematica:
% Each row is a set: [CE]  -> corresponds to [Ce]
param_sets = [
    0.0534711;
    0.115396;
    -0.115396;
    -0.0534711
];

%Define a decaywidth matrix where each row is a gamma computed with the
%diferent 4 parameters. If weverything goes well 4 of the results should be
%the same 2 by 2
Gamma = zeros(4, num_valid);  % Rows = parameter sets, columns = data points

for i = 1:num_valid
    DW = DW_values(i);
    
    for N = 1:4
        CE = param_sets(N);
        
        Gamma(N, i) = FracFpi * (CE^2 * DW);
    end
end

%To cross check we compute also the decay width integrated
% Compute M^2
M2 = M.^2;
M2_valid = M2(valid_idx);
% Perform numerical integration for each term
GammaFinal = zeros(1,4);
for N = 1:4
    GammaFinal(N) = trapz(M2_valid, Gamma(N, :));
end
% Display results
fprintf('Gamma (GeV): %.15f\n', GammaFinal(1));
fprintf('Gamma (GeV): %.15f\n', GammaFinal(2));
fprintf('Gamma (GeV): %.15f\n', GammaFinal(3));
fprintf('Gamma (GeV): %.15f\n', GammaFinal(4));

% Perform numerical integration for each term
DW_integral = trapz(M2_valid, DW_values);

GammaIntegral = zeros(1,4);
for N=1:4
    CE = param_sets(N);

    GammaIntegral(N) = FracFpi * (CE^2 * DW_integral);
end



% Define file path and name
% The following code is to create .txt files where the variables are stored
% The basename is the one mentioned above
folder_path = '/Users/sandra/Documents/Doctorat/Projectes PhD/Transicions a 2 pions/Lower order Lagrangian/Hibrids/DecayWidth';

full_path = fullfile(folder_path, file_name);

% Open the file for writing
fileID = fopen(full_path, 'w');

% Check if the file was opened successfully
if fileID == -1
    error('Could not open file for writing.');
end

% Write header
fprintf(fileID, 'Decay Width results for the different constant sets:\n');
fprintf(fileID, 'Decay width 1 (GeV)       Decay width 2 (GeV)       Decay width 3 (GeV)       Decay width 4 (GeV)\n');
fprintf(fileID,    '%.15f %.15f %.15f %.15f\n', GammaFinal(1), GammaFinal(2), GammaFinal(3), GammaFinal(4));
fprintf(fileID, '----------------------------------------------------------------------------------------\n');
fprintf(fileID, 'Mass (GeV)     Gamma1 (GeV)       Gamma2 (GeV)       Gamma3 (GeV)       Gamma4 (GeV)\n');
fprintf(fileID, '----------------------------------------------------------------------------------------\n');

for i = 1:num_valid
    fprintf(fileID,    '%.4f %.15f %.15f %.15f %.15f\n', M_valid(i), Gamma(1,i), Gamma(2,i), Gamma(3,i), Gamma(4,i));
end

% Close the file
fclose(fileID);

fprintf('Results saved to %s\n', full_path);




function termM1 = P1toSdwM1(I_Cos_R3, I_Sin_R2, I_Sin_R3, I_Sin_R4, I_Si_R, ...
    I_Si_R3, I_Si_R4, m, M, DeltaE)

    s=sqrt(DeltaE^2 - M^2);

    constM1 = ( M^2 * sqrt(1 - (4 * m^2) / (M^2))^5 ) / (1280 * s^7 );

    termCos3 = 32 * DeltaE * s * I_Cos_R3;
    termSin3 = - 4 * pi * s^2 * I_Sin_R3;
    termSin4 = 8 * DeltaE * (-8 + pi^2) * I_Sin_R4;
    termSin2 = 8 * DeltaE * s^2 * I_Sin_R2;
    termSi4 = 2 * pi^3 * DeltaE * I_Si_R4;
    termSi1 = - s^4 * I_Si_R;
    termSi3 = - pi^2 * s^2 * I_Si_R3;

    termM1 = constM1 * (termCos3 + termSin3 + termSin4 + termSin2 + termSi4 + termSi1 + termSi3)^2;

end

function termM0_a = P1toSdwM0_a(I_Cos_R3, I_Sin_R3, I_Sin_R4, I_Si_R, I_Si_R2,...
    I_Si_R3, I_Si_R4, m, M, DeltaE)

    s=sqrt(DeltaE^2 - M^2);

    constM0_a = (3/15) * sqrt(1 - (4 * m^2)/(M^2))^5 * (8*M^4 + 4*M^2*DeltaE^2 + 3*DeltaE^4) ...
        / (256 * s^11 );

    termA1 = (-32*M^2 + 4*M^2*pi^2 - 64*DeltaE^2 + 8*pi^2*DeltaE^2) * I_Sin_R4;
    termA2 = pi^3 * (M^2 + 2*DeltaE^2) * I_Si_R4;
    termA3 = pi * s^2 * (-M^2 - 2*DeltaE^2) * I_Si_R2;
    termA4 = 16 * s * (2*DeltaE^2 + M^2) * I_Cos_R3;
    termA5 = 8 * pi * DeltaE * (M^2 - DeltaE^2) * I_Sin_R3;
    termA6 = 2 * pi^2 * DeltaE * (M^2 - DeltaE^2) * I_Si_R3;
    termA7 = 2 * DeltaE * s^4 * I_Si_R;

    termM0_a = constM0_a * (termA1+ termA2 + termA3 + termA4 + termA5 + termA6 + termA7)^2;

end

function termM0_b = P1toSdwM0_b(I_Cos_R3, I_Sin_R3, I_Sin_R4, I_Si_R, I_Si_R2,...
    I_Si_R3, I_Si_R4, m, M, DeltaE)

    s=sqrt(DeltaE^2 - M^2);

    constM0_b = 3 * sqrt(1 - (4 * m^2)/(M^2)) / (256 * s^11 );

    termB1 = ( 8*M^4*(-8+pi^2) + 384*m^2*DeltaE^2 + 32*M^2*DeltaE^2 - 48*m^2*pi^2*DeltaE^2 - ...
        4*M^2*pi^2*DeltaE^2 - 64*DeltaE^4 + 8*pi^2*DeltaE^4 ) * I_Sin_R4;
    termB2 = pi^3 * ( 2*M^4 -12*m^2*DeltaE^2 - M^2*DeltaE^2 + 2*DeltaE^4 ) * I_Si_R4;
    termB3 = pi * s^2 * ( -2*M^4 +12*m^2*DeltaE^2 + M^2*DeltaE^2 - 2*DeltaE^4) * I_Si_R2;
    termB4 = 16 * s * ( 2*M^4 - (12*m^2+M^2)*DeltaE^2 + 2*DeltaE^4 ) * I_Cos_R3;
    termB5 = 8 * pi * DeltaE * ( -4*m^2*M^2 + 4*m^2*DeltaE^2 + M^2*DeltaE^2 - DeltaE^4 ) * I_Sin_R3;
    termB6 = 2 * pi^2 * DeltaE * ( -4*m^2*M^2 + 4*m^2*DeltaE^2 + M^2*DeltaE^2 - DeltaE^4 ) * I_Si_R3;
    termB7 = 2 * DeltaE * s^2 * ( 4*m^2*M^2 - 4*m^2*DeltaE^2 - M^2*DeltaE^2 + DeltaE^4 ) * I_Si_R;

    termM0_b = constM0_b * (termB1+ termB2 + termB3 + termB4 + termB5 + termB6 + termB7)^2;

end

function termM0_c = P1toSdwM0_c(I_Cos_R3, I_Sin_R3, I_Sin_R4, I_Si_R, I_Si_R2,...
    I_Si_R3, I_Si_R4, m, M, DeltaE)

    s=sqrt(DeltaE^2 - M^2);

    constM0_c = 2 * sqrt(1 - (4 * m^2)/(M^2))^3 * (2*M^2 + DeltaE^2) ...
        / (256 * s^11 );

    termA1 = (-32*M^2 + 4*M^2*pi^2 - 64*DeltaE^2 + 8*pi^2*DeltaE^2) * I_Sin_R4;
    termA2 = pi^3 * (M^2 + 2*DeltaE^2) * I_Si_R4;
    termA3 = pi * s^2 * (-M^2 - 2*DeltaE^2) * I_Si_R2;
    termA4 = 16 * s * (2*DeltaE^2 + M^2) * I_Cos_R3;
    termA5 = 8 * pi * DeltaE * (M^2 - DeltaE^2) * I_Sin_R3;
    termA6 = 2 * pi^2 * DeltaE * (M^2 - DeltaE^2) * I_Si_R3;
    termA7 = 2 * DeltaE * s^4 * I_Si_R;
    termB1 = ( 8*M^4*(-8+pi^2) + 384*m^2*DeltaE^2 + 32*M^2*DeltaE^2 - 48*m^2*pi^2*DeltaE^2 - ...
        4*M^2*pi^2*DeltaE^2 - 64*DeltaE^4 + 8*pi^2*DeltaE^4 ) * I_Sin_R4;
    termB2 = pi^3 * ( 2*M^4 -12*m^2*DeltaE^2 - M^2*DeltaE^2 + 2*DeltaE^4 ) * I_Si_R4;
    termB3 = pi * s^2 * ( -2*M^4 +12*m^2*DeltaE^2 + M^2*DeltaE^2 - 2*DeltaE^4) * I_Si_R2;
    termB4 = 16 * s * ( 2*M^4 - (12*m^2+M^2)*DeltaE^2 + 2*DeltaE^4 ) * I_Cos_R3;
    termB5 = 8 * pi * DeltaE * ( -4*m^2*M^2 + 4*m^2*DeltaE^2 + M^2*DeltaE^2 - DeltaE^4 ) * I_Sin_R3;
    termB6 = 2 * pi^2 * DeltaE * ( -4*m^2*M^2 + 4*m^2*DeltaE^2 + M^2*DeltaE^2 - DeltaE^4 ) * I_Si_R3;
    termB7 = 2 * DeltaE * s^2 * ( 4*m^2*M^2 - 4*m^2*DeltaE^2 - M^2*DeltaE^2 + DeltaE^4 ) * I_Si_R;

    termM0_c = constM0_c * (termA1+ termA2 + termA3 + termA4 + termA5 + termA6 + termA7) * ...
        (-termB1 - termB2 - termB3 - termB4 - termB5 - termB6 - termB7);

end



%VALOR DE LA MASSA
%m=1.4702; charm
%m=4.8802; bottom
function x0=m_q
global v0
x0=v0;
end

function setm_q(val0)
global v0
v0=val0;
end

% VALOR DEL SPIN
function x7=spin
global v7
x7=v7;
end 

function setspin(val7)
global v7
v7 = val7;
end

% VALOR DEL r0
function x3=r0
global v3
x3=v3;
end 

function setr0(val3)
global v3
v3 = val3;
end
