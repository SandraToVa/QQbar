
%Document from DecayWidth_QQ but including uncertanties
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
ComputationIE = TransitionsAdded('QQStoS');  
%spin
s=1;

[I_if_square_cell, DeltaE, M] = ComputationIE(3, 2, 0.28, 0.453, 50, s);



file_name = 'GammaVSmass_QQcharm_3s-2s_uncert.txt';

% ============================================================
% Data importation and preliminar definitions
% ============================================================

%Obtain the I_if^2 from the cell
I_if_square = I_if_square_cell{1,1};
I_if_c_square = I_if_square_cell{1,2};
I_if_0c_square = I_if_square_cell{1,3};
I_if_1_square = I_if_square_cell{1,4};
I_if_2_square = I_if_square_cell{1,5};

% Value of pion mass 140MeV
m = 0.14;
% Value of the global constant that divides everithing in GeV. Fpi=92MeV
FracFpi = 1/(2^2 * pi^4 * 0.092^4);

% Filter values where DeltaE > M
valid_idx = (DeltaE > M);
M_valid = M(valid_idx);
I_if_valid = I_if_square(valid_idx);
I_if_c_valid = I_if_c_square(valid_idx);
I_if_0c_valid = I_if_0c_square(valid_idx);
I_if_1_valid = I_if_1_square(valid_idx);
I_if_2_valid = I_if_2_square(valid_idx);


num_valid = length(M_valid);
Cm2_values = zeros(1, num_valid);  % Preallocate for speed
CM2_values = zeros(1, num_valid);  
CECM_values = zeros(1, num_valid);  
CmCM_values = zeros(1, num_valid);  
CmCE_values = zeros(1, num_valid);  
CE2_values = zeros(1, num_valid);  


for i = 1:num_valid
    Cm2_values(i) = Cm2Term(I_if_valid(i), m, M_valid(i), DeltaE);
    CM2_values(i) = CM2Term(I_if_valid(i), m, M_valid(i), DeltaE);
    CECM_values(i) = CECMTerm(I_if_valid(i), I_if_0c_valid(i), m, M_valid(i), DeltaE);
    CmCM_values(i) = CmCMTerm(I_if_valid(i), m, M_valid(i), DeltaE);
    CmCE_values(i) = CmCETerm(I_if_valid(i), I_if_0c_valid(i), m, M_valid(i), DeltaE);
    CE2_values(i) = CE2Term(I_if_valid(i), I_if_0c_valid(i), I_if_c_valid(i), I_if_1_valid(i), I_if_2_valid(i), m, M_valid(i), DeltaE);
end

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
% Means of the parameters from mathematica [c_pi, c_pipi, c_E] 
means_matrix = [ -0.568834,  0.0405707,   0.0232213;  % Red set
                 -0.566163, -0.00900834,  0.0722651;   % Green set
                  0.566163,  0.00900834, -0.0722651;   % Green set
                  0.568834, -0.0405707,  -0.0232213 ]; % Red set

% Covariance Matrices in a Cell Array
cov_matrices = cell(1, 4);

% Example for Set 1 (red)
cov_matrices{1} = [ 0.000668187,  0.0000415763, -0.000139412;
                    0.0000415763, 0.0000151738, -0.0000236754;
                   -0.000139412, -0.0000236754,  0.0000483156 ]; 
% Example for Set 2 (green)
cov_matrices{2} = [ 0.000684674,   -0.0000934895, -7.59907*10^-6;
                   -0.0000934895,   0.0000279475, -0.0000107644;
                   -7.59907*10^-6, -0.0000107644,  0.0000103614 ]; 
% Example for Set 3 (green)
cov_matrices{3} = cov_matrices{2}; 
% Example for Set 4 (red)
cov_matrices{4} = cov_matrices{1} ; 


% ============================================================
% MonteCarlo calculation
% ============================================================

for n = 1:4
    % 1. Get stats for this specific set
    mu = means_matrix(n, :);       % [c_pi, c_pipi, c_E] for set n
    Sigma = cov_matrices{n};       % Your pre-computed 3x3 covariance matrix
    
    % 2. Generate nMC correlated samples
    nMC = 100000;
    samples = mvnrnd(mu, Sigma, nMC); 
    
    % 3. Pre-allocate temporary storage for this set's curves
    % Rows are iterations, Columns are x-points
    temp_gamma = zeros(nMC, num_valid);
    temp_integrals = zeros(nMC, 1);
    
    for m = 1:nMC
        % Extract the specific a, b, c for this iteration
        Cm = samples(m, 1); %c_pi
        CM = samples(m, 2); %c_pipi
        CE = samples(m, 3); %c_E
        
        % Compute the gamma for all i points
        for i = 1:num_valid
            CE2 = CE2_values(i);
            CmCE = CmCE_values(i);
            CmCM = CmCM_values(i);
            CECM = CECM_values(i);
            CM2 = CM2_values(i);
            Cm2 = Cm2_values(i);

            gamma_value = FracFpi * (CE^2 * CE2 + CE*Cm * CmCE + Cm*CM * CmCM + CM^2 * CM2 + Cm^2 * Cm2 + CE*CM* CECM);

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
folder_path = '/Users/sandra/Documents/Doctorat/Projectes PhD/Transicions a 2 pions/Lower order Lagrangian/Quarkonium/Uncertanties';
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


function term = Cm2Term(I_if_square, m, M, DeltaE)
    term = sqrt(DeltaE^2 - M^2) * sqrt(1 - (4 * m^2) / (M^2)) ...
        * (2* I_if_square * m^4);
end

function term = CM2Term(I_if_square, m, M, DeltaE)
    term = sqrt(DeltaE^2 - M^2) * sqrt(1 - (4 * m^2) / (M^2)) ...
        * (2* I_if_square * M^4);
end

function term = CECMTerm(I_if_square, I_if_0c_square, m, M, DeltaE)
    term1 = I_if_0c_square * (4*m^2*M^2 + 2*M^4 - 2*(2*m^2 + M^2)*DeltaE^2);
    term2 = I_if_square * (4*m^2*M^2 - M^4 + 2*(2*m^2 + M^2)*DeltaE^2);
    term = sqrt(DeltaE^2 - M^2) * sqrt(1 - (4 * m^2) / (M^2)) ...
        * (4 /3) * (term1 + term2);
end

function term = CmCMTerm(I_if_square, m, M, DeltaE)
    term = sqrt(DeltaE^2 - M^2) * sqrt(1 - (4 * m^2) / (M^2)) ...
        * (4 * I_if_square * m^2 * M^2);
end

function term = CmCETerm(I_if_square, I_if_0c_square, m, M, DeltaE)
    term1 = I_if_0c_square * (2*M^4 +4*m^2*M^2 - 2*(2*m^2 + M^2)*DeltaE^2 );
    term2 = I_if_square * ( -M^4 +4*m^2*M^2 + 2*(2*m^2 + M^2)*DeltaE^2  );
    term = sqrt(DeltaE^2 - M^2) * sqrt(1 - (4 * m^2) / (M^2)) ...
        * ((4 * m^2)/(3*M^2))  * (term1 + term2);
end

function term = CE2Term(I_if_square, I_if_0c_square, I_if_c_square, I_if_1_valid, I_if_2_square, m, M, DeltaE)
    term0= 2* ( (4/(15*M^2)) * I_if_1_valid * (-4*m^2 + M^4)^2 * DeltaE^2 );
    term1= (2/(15*M^4)) * I_if_c_square *( 32*m^4*M^4 + 24*m^2*M^6 + 7*M^8 - 4*M^2*DeltaE^2*(8*m^4 ...
        + 16*m^2*M^2 + 3*M^4) + 8*(6*m^4 + 2*m^2*M^2 + M^4)*DeltaE^4 );
    term2= (4/(15*M^4)) * I_if_square * ( M^4*(-4*m^2 + M^2)^2 - 4*M^2*(M^2 - 4*m^2)*(m^2+M^2)*DeltaE^2 ...
        + 4*(6*m^4 + 2*m^2*M^2 + M^4)*DeltaE^4 );
    term3= 2 * (1/(15)) * I_if_2_square * (-4*m^2 + M^2)^2  ;
    term4= (4/(15*M^4)) * I_if_0c_square * ( 32*m^4*M^4 + 4*m^2*M^6 - 3*M^8 + 10*M^4*(2*m^2 + M^2)*DeltaE^2 ...
        - 8*(6*m^4 + 2*m^2*M^2 + M^4)*DeltaE^4 );
    term = sqrt(DeltaE^2 - M^2) * sqrt(1 - (4 * m^2) / (M^2)) * ...
        (term0 + term1 + term2 + term3 + term4);
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



