%Document from DecayWidth_P1toS but including uncertanties
%The uncertanties calculation ist done in a MonteCarlo procedure and uses
%the mathematica data from PhD_QuarkoniumUncertanties such as the mean, the
%standard deviation and the covariant matrices computed there.

% ============================================================
% Calculation of the transition
% ============================================================


%massa
load("dades.mat","m_c","m_b")
setm_q(m_c)

%Change this depending on the transition
ComputationIE = TransitionsP1toS('HQp1tS_full');      

[I_Cos_cell, I_Sin_cell, I_Si_cell, DeltaE, M] = ComputationIE(2, 2, 0.28, 0.467, 50);



file_name = 'GammaVSmass_HQcharm_1p1-2s_uncert_full.txt';


% ============================================================
% Data importation and preliminar definitions
% ============================================================

% Value of pion mass 140MeV
mpi = 0.14;
% Value of the global constant that divides everithing in GeV. Fpi=92MeV
FracFpi = 1/(pi^3 * 0.092^4 * 0.187);
%0.187GeV^2 is the string tension

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


% ============================================================
% Definition of the uncertanties in the parameters
% ============================================================

% M = Number of Monte Carlo iterations (e.g., 1000)
% N = 4 (Your 4 sets of parameters)
% num_val = length of x vector

gamma_mean = zeros(num_valid, 4);
gamma_std  = zeros(num_valid, 4);
int_mean   = zeros(4, 1);
int_std    = zeros(4, 1);

% Definition of the statistics from mathematica
% Means of the parameters from mathematica [lambda, eta, frac] 
means_matrix = [ 0.0346988,  0.0928852,  -0.44125;  % Red set
                -0.162547,   0.28906,    -0.439649;   % Green set
                 0.162547,  -0.28906,     0.439649;   % Green set
                -0.0346988, -0.0928852,   0.44125 ]; % Red set

% Covariance Matrices in a Cell Array
cov_matrices = cell(1, 4);

% Example for Set 1 (red)
cov_matrices{1} = [ 0.000443361, -0.000575928, 0.000229409;
                    -0.000575928, 0.00077305, -0.000360526;
                    0.000229409, -0.000360526, 0.000341399 ]; 
% Example for Set 2 (green)
cov_matrices{2} = [  0.000239351, -0.000169006, -0.000101436;
                    -0.000169006,  0.000165782, -0.0000336203;
                    -0.000101436, -0.0000336203, 0.00034744 ]; 
% Example for Set 3 (green)
cov_matrices{3} = cov_matrices{2}; 
% Example for Set 4 (red)
cov_matrices{4} = cov_matrices{1} ; 

% ============================================================
% MonteCarlo calculation
% ============================================================

for n = 1:4
    % 1. Get stats for this specific set
    mu = means_matrix(n, :);       % [lambda, eta, frac]for set n
    Sigma = cov_matrices{n};       % Your pre-computed 3x3 covariance matrix
    
    % 2. Generate M correlated samples
    nMC = 100000;
    samples = mvnrnd(mu, Sigma, nMC); 
    
    % 3. Pre-allocate temporary storage for this set's curves
    % Rows are iterations, Columns are x-points
    temp_gamma = zeros(nMC, num_valid);
    temp_integrals = zeros(nMC, 1);
    
    for m = 1:nMC
        % Extract the specific a, b, c for this iteration
        lamb = samples(m, 1); %c_pi
        eta =  samples(m, 2); %c_pipi
        frac = samples(m, 3); %c_E
        
        % Compute the gamma for all i points
        for i = 1:num_valid
            %multiplicat per eta^2
            dwMin1 = P1toSdwM1(I_Cos_R3(i), I_Sin_R2(i), I_Sin_R3(i), I_Sin_R4(i), I_Si_R(i), ...
            I_Si_R2(i), I_Si_R3(i), I_Si_R4(i), mpi, M_valid(i), DeltaE, eta);
    
            %multiplicat per eta^2
            dwMin0a = P1toSdwM0_a(I_Cos_R3(i), I_Sin_R2(i), I_Sin_R3(i), I_Sin_R4(i), I_Si_R0(i), I_Si_R(i), I_Si_R2(i),...
            I_Si_R3(i), I_Si_R4(i), mpi, M_valid(i), DeltaE, eta);
    
            %
            dwMin0b = P1toSdwM0_b(I_Cos_R3(i), I_Sin_R2(i), I_Sin_R3(i), I_Sin_R4(i), I_Si_R0(i), I_Si_R(i), I_Si_R2(i),...
            I_Si_R3(i), I_Si_R4(i), mpi, M_valid(i), DeltaE, lamb, eta, frac);

            dwMin0c = P1toSdwM0_c(I_Cos_R3(i), I_Sin_R2(i), I_Sin_R3(i), I_Sin_R4(i), I_Si_R0(i), I_Si_R(i), I_Si_R2(i),...
            I_Si_R3(i), I_Si_R4(i), mpi, M_valid(i), DeltaE, lamb, eta, frac);

            dwMin0 = (dwMin0a + dwMin0b + dwMin0c);

            % Total decay width is the average over initial states Min and addition
            % to final Min=0. Whe have Min=+1,-1,0 and Jin=1
            DW = ( dwMin1*2 + dwMin0 )/3;

            gamma_value = FracFpi * DW;

            temp_gamma(m, i) = gamma_value;

        end
        
        % Compute the numerical integral
        % Compute M^2
        M2 = M.^2;
        M2_valid = M2(valid_idx);
        % Perform numerical integration for each term
        temp_integrals(m) = trapz(M2_valid, temp_gamma(m, :) );
    end
    
    % 4. Aggregate statistics for Set N
    %Value of gamma and their uncertanties at each point
    gamma_mean(:, n) = mean(temp_gamma, 1)';
    gamma_std(:, n)  = std(temp_gamma, 0, 1)';
    
    %Value of the integrated gamma and their uncertanties 
    int_mean(n) = mean(temp_integrals);
    int_std(n)  = std(temp_integrals);
end


% Display results
for n = 1:4
    fprintf('Set%d %.15f %.15f\n', n, int_mean(n), int_std(n));
end



% ============================================================
% Data saving (Updated for Monte Carlo Results)
% ============================================================
folder_path = '/Users/sandra/Documents/Doctorat/Projectes PhD/Transicions a 2 pions/Lower order Lagrangian/Hibrids/Uncertanties';
full_path = fullfile(folder_path, file_name);
fileID = fopen(full_path, 'w');

if fileID == -1
    error('Could not open file for writing.');
end

% --- Section 1: Integrated Results (Total Decay Widths) ---
fprintf(fileID, '# INTEGRATED DECAY WIDTHS\n');
fprintf(fileID, '# Set_N | Mean_Int | Std_Int\n');
for n = 1:4
    fprintf(fileID, 'Set%d %.15f %.15f\n', n, int_mean(n), int_std(n));
end

fprintf(fileID, '\n'); % Blank line separator

% --- Section 2: Curve Data (Mean and Std Dev vs Mass) ---
% Using a "Long Format" makes Python plotting much easier
fprintf(fileID, '# CURVE DATA\n');
fprintf(fileID, 'Mass(GeV) Set_N Gamma_Mean(GeV) Gamma_Std(GeV)\n');

for n = 1:4
    for i = 1:num_valid
        % We write: Mass, Set Number, Mean, StdDev
        fprintf(fileID, '%.4f %d %.15f %.15f\n', ...
                M_valid(i), n, gamma_mean(i, n), gamma_std(i, n));
    end
end

fclose(fileID);
fprintf('Monte Carlo results saved to %s\n', full_path);


% ============================================================
% Functions to compute the decay width
% ============================================================



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

