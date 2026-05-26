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

%%

file_name = 'GammaVSmass_HQcharm_1p1-2s_full56.txt';

% Value of pion mass 140MeV
m = 0.14;
% Value of the global constant that divides everithing in GeV. Fpi=92MeV
FracFpi = 1/(pi^3 * 0.092^4 * 0.187) ;
%0.092 és f_pi i 0.187GeV^2 is the string tension

%Obtain the I from the cell
I_Cos_R3 = I_Cos_cell{1, 1}; 
I_Sin_R2 = I_Sin_cell{1, 1};
I_Sin_R3 = I_Sin_cell{1, 2};
I_Sin_R4 = I_Sin_cell{1, 3};
I_Si_R = I_Si_cell{1, 1};
I_Si_R2 = I_Si_cell{1, 2};
I_Si_R3 = I_Si_cell{1, 3};
I_Si_R4 = I_Si_cell{1, 4};
I_Si_R0 = I_Si_cell{1, 5};


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
I_Si_R0 = I_Si_R0(valid_idx);


num_valid = length(M_valid);
DW_values = zeros(4, num_valid);

% We compute the DW_values for each set of lamb, eta, frac:
% Each row is a set: [lambda, eta, lambda'/2B0] = [lamb, eta, frac]
param_sets = [  0.0346988,  0.0928852,  -0.44125;  % Red set
               -0.162547,   0.28906,    -0.439649;   % Green set
                0.162547,  -0.28906,     0.439649;   % Green set
               -0.0346988, -0.0928852,   0.44125 ]; % Red set

%Define a decaywidth matrix where each row is a gamma computed with the
%diferent 4 parameters. If weverything goes well 4 of the results should be
%the same 2 by 2
Gamma = zeros(4, num_valid);  % Rows = parameter sets, columns = data points


for i = 1:num_valid
    for N = 1:4
        lamb = param_sets(N,1);
        eta = param_sets(N,2);
        frac = param_sets(N,3);

        %multiplicat per eta^2
        dwMin1 = P1toSdwM1(I_Cos_R3(i), I_Sin_R2(i), I_Sin_R3(i), I_Sin_R4(i), I_Si_R(i), ...
        I_Si_R2(i), I_Si_R3(i), I_Si_R4(i), m, M_valid(i), DeltaE, eta);
    
        %multiplicat per eta^2
        dwMin0a = P1toSdwM0_a(I_Cos_R3(i), I_Sin_R2(i), I_Sin_R3(i), I_Sin_R4(i), I_Si_R0(i), I_Si_R(i), I_Si_R2(i),...
        I_Si_R3(i), I_Si_R4(i), m, M_valid(i), DeltaE, eta);
    
        %
        dwMin0b = P1toSdwM0_b(I_Cos_R3(i), I_Sin_R2(i), I_Sin_R3(i), I_Sin_R4(i), I_Si_R0(i), I_Si_R(i), I_Si_R2(i),...
        I_Si_R3(i), I_Si_R4(i), m, M_valid(i), DeltaE, lamb, eta, frac);

        dwMin0c = P1toSdwM0_c(I_Cos_R3(i), I_Sin_R2(i), I_Sin_R3(i), I_Sin_R4(i), I_Si_R0(i), I_Si_R(i), I_Si_R2(i),...
        I_Si_R3(i), I_Si_R4(i), m, M_valid(i), DeltaE, lamb, eta, frac);

        dwMin0 = (dwMin0a + dwMin0b + dwMin0c);
   
        % Total decay width is the average over initial states Min and addition
        % to final Min=0. Whe have Min=+1,-1,0 and Jin=1
        DW_values(N, i) = ( dwMin1*2 + dwMin0 )/3;

        %Total gamma with constants
        Gamma(N, i) = FracFpi * DW_values(N, i);
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



%Corregit afegint tant g5 com g6
function termM1 = P1toSdwM1(I_Cos_R3, I_Sin_R2, I_Sin_R3, I_Sin_R4, I_Si_R, ...
    I_Si_R2, I_Si_R3, I_Si_R4, m, M, DeltaE, eta)

    s=sqrt(DeltaE^2 - M^2);

    constM1 = ( M^2 * sqrt(1 - (4 * m^2) / (M^2))^5 * eta^2) / (5120 * s^7 );

    termCos3 = 48 * DeltaE * s * I_Cos_R3;
    termSin3 = - 8 * pi * s^2 * I_Sin_R3;
    termSin4 = 4 * DeltaE * (-24 + 3*pi^2) * I_Sin_R4;
    termSin2 = 16 * DeltaE * s^2 * I_Sin_R2;
    termSi4 = 3 * pi^3 * DeltaE * I_Si_R4;
    termSi1 = - 2*s^4 * I_Si_R;
    termSi3 = - 2*pi^2 * s^2 * I_Si_R3;
    termSi2 = pi * DeltaE * s^2 * I_Si_R2;

    termM1 = constM1 * (termCos3 + termSin3 + termSin4 + termSin2 + termSi4 + termSi1 + termSi3 + termSi2)^2;

end

%Corregit afegint tant g5 com g6
function termM0_a = P1toSdwM0_a(I_Cos_R3, I_Sin_R2, I_Sin_R3, I_Sin_R4, I_Si_R0, I_Si_R, I_Si_R2,...
    I_Si_R3, I_Si_R4, m, M, DeltaE, eta)

    s=sqrt(DeltaE^2 - M^2);

    constM0_a = sqrt(1 - (4 * m^2)/(M^2))^5 * (8*M^4 + 4*M^2*DeltaE^2 + 3*DeltaE^4) * eta^2 ...
        / (20480 * pi^2 * s^11 );
 

    termA1 = 48 * pi * s * (2*DeltaE^2 + M^2) * I_Cos_R3;
    termA2 = 2*pi*(-8+pi^2)*(12*DeltaE^2 + 6*M^2) * I_Sin_R4;
    termA3 = -32 * pi^2 * s^2 * DeltaE * I_Sin_R3;
    termA4 = 2 * pi * s^2 * (3*s^2 + M^2 + DeltaE^2) * I_Sin_R2;
    termA5 = 3*pi^4 * (M^2 + 2*DeltaE^2) * I_Si_R4;
    termA6 = -8 * pi^3 * DeltaE * s^2 * I_Si_R3;
    termA7 = -4 * pi^2 * s^2 * (M^2 + DeltaE^2) * I_Si_R2;
    termA8 = 8 * DeltaE * s^4 * pi * I_Si_R;
    termA9 = - s^4 * (-M^2 + 2*DeltaE^2) * I_Si_R0;

    termM0_a = constM0_a * (termA1+ termA2 + termA3 + termA4 + termA5 + termA6 + termA7 + termA8 + termA9)^2;

end

%Corregit afegint tant g5 com g6
function termM0_b = P1toSdwM0_b(I_Cos_R3, I_Sin_R2, I_Sin_R3, I_Sin_R4, I_Si_R0, I_Si_R, I_Si_R2,...
    I_Si_R3, I_Si_R4, m, M, DeltaE, lamb, eta, frac)

    s=sqrt(DeltaE^2 - M^2);

    constM0_b = 15 * sqrt(1 - (4 * m^2)/(M^2)) / (20480 * pi^2 * s^11 );
   
    termB1 = s^4 * (-8*frac*m^2*s^2 - 2*DeltaE^4*eta - 8*m^2*M^2*(eta+lamb) + ...
        2*M^4*(eta + 2*lamb) + 4*m^2*DeltaE^2*(3*eta + 2*lamb) - M^2*DeltaE^2*(eta + 4*lamb) ) * I_Si_R0;
    termB2 = 48 * pi * s * ( 2*M^4 - (12*m^2 + M^2)*DeltaE^2 + 2*DeltaE^4 ) * eta * I_Cos_R3;
    termB3 = -32 * pi^2 * s^2 * ( -4*m^2 + M^2 + s^2 ) * DeltaE * eta * I_Sin_R3;
    termB4 = - 8 * pi^3 * s^2 * DeltaE * ( -4*m^2 + DeltaE^2 ) * eta * I_Si_R3;
    termB5 = 8 * pi * s^4 * DeltaE * (- 4*m^2 + DeltaE^2 ) * eta * I_Si_R;
    termB6 = 12 * pi * (-8 + pi^2) * (2*s^4 + 3 * (-4*m^2 + M^2)*DeltaE^2 ) * eta * I_Sin_R4;
    termB7 = 3 * pi^4 * ( 2*M^4 - ( 12*m^2 + M^2) * DeltaE^2 + 2*DeltaE^4 ) * eta * I_Si_R4;
    termB8 = 4 * pi* s^2 * (8*frac*m^2*s^2 + 2*M^2*s^2*eta + M^2*DeltaE^2*eta + 2*s^2*DeltaE^2*eta + ...
        4*M^2*s^2*lamb - 4*m^2*(DeltaE^2*eta + 2*s^2*(eta + lamb) ) ) * I_Sin_R2;
    termB9 = 4* pi^2 * s^2 * ( 2*frac*m^2*s^2 + ( -2*M^4 + 2*m^2*(M^2 + 3*DeltaE^2) + M^2*DeltaE^2 - DeltaE^4)*eta - ...
        s^2*(2*m^2 - M^2)*lamb ) * I_Si_R2;

    termM0_b = constM0_b * (termB1+ termB2 + termB3 + termB4 + termB5 + termB6 + termB7 + termB8 + termB9)^2;

end

%Corregit afegint tant g5 com g6
function termM0_c = P1toSdwM0_c(I_Cos_R3, I_Sin_R2, I_Sin_R3, I_Sin_R4, I_Si_R0, I_Si_R, I_Si_R2,...
    I_Si_R3, I_Si_R4, m, M, DeltaE, lamb, eta, frac)

    s=sqrt(DeltaE^2 - M^2);

    constM0_c = - 10 * sqrt(1 - (4 * m^2)/(M^2))^3 * (2*M^2 + DeltaE^2) * eta ...
        / (20480 * pi^2 * s^11 );
  
    termA1 = 48 * pi * s * (2*DeltaE^2 + M^2) * I_Cos_R3;
    termA2 = 2*pi*(-8+pi^2)*(12*DeltaE^2 + 6*M^2) * I_Sin_R4;
    termA3 = -32 * pi^2 * s^2 * DeltaE * I_Sin_R3;
    termA4 = 2 * pi * s^2 * (3*s^2 + M^2 + DeltaE^2) * I_Sin_R2;
    termA5 = 3*pi^4 * (M^2 + 2*DeltaE^2) * I_Si_R4;
    termA6 = -8 * pi^3 * DeltaE * s^2 * I_Si_R3;
    termA7 = -4 * pi^2 * s^2 * (M^2 + DeltaE^2) * I_Si_R2;
    termA8 = 8 * DeltaE * s^4 * pi * I_Si_R;
    termA9 = - s^4 * (-M^2 + 2*DeltaE^2) * I_Si_R0;
    termB1 = s^4 * (-8*frac*m^2*s^2 - 2*DeltaE^4*eta - 8*m^2*M^2*(eta+lamb) + ...
        2*M^4*(eta + 2*lamb) + 4*m^2*DeltaE^2*(3*eta + 2*lamb) - M^2*DeltaE^2*(eta + 4*lamb) )*I_Si_R0;
    termB2 = 48 * pi * s * ( 2*M^4 - (12*m^2+M^2)*DeltaE^2 + 2*DeltaE^4 ) * eta * I_Cos_R3;
    termB3 = -32 * pi^2 * s^2 * ( -4*m^2 + M^2 + s^2 ) * DeltaE * eta * I_Sin_R3;
    termB4 = -8 * pi^3 * s^2 * DeltaE * ( -4*m^2 + DeltaE^2 ) * eta * I_Si_R3;
    termB5 = 8 * pi * s^4 * DeltaE * (- 4*m^2 + DeltaE^2 ) * eta * I_Si_R;
    termB6 = 12 * pi * (-8 + pi^2) * (2*s^4 + 3 * (-4*m^2 + M^2)*DeltaE^2 ) * eta * I_Sin_R4;
    termB7 = 3 * pi^4 * ( 2*M^4 - ( 12*m^2 + M^2) * DeltaE^2 + 2*DeltaE^4 ) * eta * I_Si_R4;
    termB8 = 4 * pi* s^2 * (8*frac*m^2*s^2 + 2*M^2*s^2*eta + M^2*DeltaE^2*eta + 2*s^2*DeltaE^2*eta + ...
        4*M^2*s^2*lamb - 4*m^2*(DeltaE^2*eta + 2*s^2*(eta + lamb) ) ) * I_Sin_R2;
    termB9 = 4* pi^2 * s^2 * ( 2*frac*m^2*s^2 + ( -2*M^4 + 2*m^2*(M^2 + 3*DeltaE^2) + M^2*DeltaE^2 - DeltaE^4)*eta - ...
        s^2*(2*m^2 - M^2)*lamb ) * I_Si_R2;

    termM0_c = constM0_c * (termA1+ termA2 + termA3 + termA4 + termA5 + termA6 + termA7 + termA8 + termA9) * ...
        (termB1 + termB2 + termB3 + termB4 + termB5 + termB6 + termB7 + termB8 + termB9);

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
