%Program to plot the decay width with respect to the dipion mass for each
%transition and also for each set of coeficient values per a QQ->QQ
%transitions

%The coeficients have been computed separately with
%PhDSolvingDecayConstants.nb

%This program is special in the sence that insted that working with the
%functions that are given in the system of equations of solving Cm, CM and
%Ce. We work with the equations that are given for spin 1 in the Spin 1 for
%quarkonium document in the overleaf.

% Therefore the diferent parts in the program for transition added involving
% p->p are for the I_if without being squared! and diferent forn the other
% parts 

%Also is specific in the case of s=1 to transitions 1++ -> 1++

% We work in GeV!!

%Valor r0
setr0(3.964)
%massa
load("dades.mat","m_c","m_b")
setm_q(m_c)

%Change this depending on the transition

ComputationIE = TransitionsAdded('QQPtoP');  

%spin
s=1;

[I_if_cell, DeltaE, M] = ComputationIE(2, 1, 0.28, 0.474, 50, s);

%%

file_name = 'GammaVSmass_QQcharm_2p-1p_inf.txt';
base_filename = 'QQcharm_2p-1p';

%Obtain the I_if from the cell
I_if0 = I_if_cell{1,1};
I_if_c0 = I_if_cell{1,2};
I_if1 = I_if_cell{1,3};
I_if_c1 = I_if_cell{1,4};
I_if_1 = I_if_cell{1,5};
I_if_2 = I_if_cell{1,6};


% Value of pion mass 140MeV
m = 0.14;
% Value of the global constant that divides everithing in GeV. Fpi=92MeV
FracFpi = 1/(2^2 * pi^4 * 0.092^4);


% Filter values where DeltaE > M
valid_idx = (DeltaE > M);
M_valid = M(valid_idx);
I_if0_valid = I_if0(valid_idx);
I_if_c0_valid = I_if_c0(valid_idx);
I_if1_valid = I_if1(valid_idx);
I_if_c1_valid = I_if_c1(valid_idx);
I_if_1_valid = I_if_1(valid_idx);
I_if_2_valid = I_if_2(valid_idx);

%This belongs to the adition of transitions amplitude before squaring them
sa_if_01 = I_if0_valid + I_if1_valid; %T(1-1,1-1) + T(10,10) and %T(11,11) + T(10,10)
sa_if_1N = 2 .* I_if1_valid; %T(1-1,1-1) + T(11,11)
sa_if_c_01 = I_if_c0_valid + I_if_c1_valid; %T(1-1,1-1) + T(10,10) and %T(11,11) + T(10,10)
sa_if_c_1N = 2 .* I_if_c1_valid; %T(1-1,1-1) + T(11,11)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now we do the spin average for each square all with valid length
%this is particular for 1++ states!!!
% IMPORTANT!!!!!
I_if_square = ( sa_if_01.^2 + 2.* sa_if_1N.^2) ./ 12;
I_if_c_square = ( sa_if_c_01.^2 + 2.* sa_if_c_1N.^2 ) ./ 12;
I_if_0c_square = ( sa_if_01 .* sa_if_c_01 + 2.* sa_if_1N .* sa_if_c_1N ) ./ 12; 
I_if_1_square = (I_if_1_valid.^2 ) .* (4/12); %4 comes from all possible transitions involved
I_if_2_square = (I_if_2_valid.^2 ) .* (2/12); %2 comes from all possible transitions involved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_valid = length(M_valid);
Cm2_values = zeros(1, num_valid);  % Preallocate for speed
CM2_values = zeros(1, num_valid);  
CECM_values = zeros(1, num_valid);  
CmCM_values = zeros(1, num_valid);  
CmCE_values = zeros(1, num_valid);  
CE2_values = zeros(1, num_valid);  


for i = 1:num_valid
    Cm2_values(i) = Cm2Term(I_if_square(i), m, M_valid(i), DeltaE);
    CM2_values(i) = CM2Term(I_if_square(i), m, M_valid(i), DeltaE);
    CECM_values(i) = CECMTerm(I_if_square(i), I_if_0c_square(i), m, M_valid(i), DeltaE);
    CmCM_values(i) = CmCMTerm(I_if_square(i), m, M_valid(i), DeltaE);
    CmCE_values(i) = CmCETerm(I_if_square(i), I_if_0c_square(i), m, M_valid(i), DeltaE);
    CE2_values(i) = CE2Term(I_if_square(i), I_if_0c_square(i), I_if_c_square(i), I_if_1_square(i), I_if_2_square(i), m, M_valid(i), DeltaE);
end

%Param sets obtingut al resoldre el sistema d'equacions
param_sets = [
   -0.569395,  0.0413287,   0.0225432;
   -0.566654, -0.00973277,  0.0730556;
    0.566654,  0.00973277, -0.0730556;
    0.569395, -0.0413287,  -0.0225432
];


%Define a decaywidth matrix where each row is a gamma computed with the
%diferent 4 parameters. If weverything goes well 4 of the results should be
%the same 2 by 2
Gamma = zeros(4, num_valid);  % Rows = parameter sets, columns = data points

for i = 1:num_valid
    CE2 = CE2_values(i);
    CmCE = CmCE_values(i);
    CmCM = CmCM_values(i);
    CECM = CECM_values(i);
    CM2 = CM2_values(i);
    Cm2 = Cm2_values(i);
    
    for N = 1:4
        Cm = param_sets(N, 1);
        CM = param_sets(N, 2);
        CE = param_sets(N, 3);
        
        Gamma(N, i) = FracFpi * (CE^2 * CE2 + CE*Cm * CmCE + Cm*CM * CmCM + CM^2 * CM2 + Cm^2 * Cm2 + CE*CM* CECM);
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
Cm2_integral = trapz(M2_valid, Cm2_values);
CM2_integral = trapz(M2_valid, CM2_values);
CECM_integral = trapz(M2_valid, CECM_values);
CmCM_integral = trapz(M2_valid, CmCM_values);
CmCE_integral = trapz(M2_valid, CmCE_values);
CE2_integral = trapz(M2_valid, CE2_values);

GammaIntegral = zeros(1,4);
for N=1:4
    Cm = param_sets(N,1);
    CM = param_sets(N,2);
    CE = param_sets(N,3);

    GammaIntegral(N) = FracFpi * (CE^2 * CE2_integral + CE*Cm * CmCE_integral + Cm*CM * CmCM_integral + CM^2 * CM2_integral + Cm^2 * Cm2_integral + CE*CM* CECM_integral);
end



% Define file path and name
% The following code is to create .txt files where the variables are stored
%file_name = 'GammaVSmass_QQbottom_5s-3s.txt';
folder_path = '/Users/sandra/Documents/Doctorat/Projectes PhD/Transicions a 2 pions/Lower order Lagrangian/Quarkonium/DecayWidth';

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


%%
% ============================================================
% Write I_if_*_valid functions vs M_valid into .txt files
% ============================================================

% Example variables:
% M_valid, I_if_valid, I_if_c_valid, I_if_0c_valid,
% I_if_1_valid, I_if_2_valid

% ---- Set output folder and base filename ----
folder_path_new = '/Users/sandra/Documents/Doctorat/Projectes PhD/Transicions a 2 pions/Lower order Lagrangian/Quarkonium/TransitionData';
%base_filename = 'QQbottom_5s-3s';

% ---- Create a structure to simplify looping ----
dataStruct = {
    'I_if_valid',   '[R]^2',                    '.txt';
    'I_if_c_valid', '[R^0]^2',                  '_c.txt';
    'I_if_0c_valid','[R * R^0]',                '_0c.txt';
    'I_if_1_valid', '[R^1]^2 and [R^{-1}]^2',   '_1.txt';
    'I_if_2_valid', '[R^2]^2 and [R^{-2}]^2',   '_2.txt'
};

% ---- Loop over each dataset ----
for k = 1:size(dataStruct,1)
    varName  = dataStruct{k,1};
    header   = dataStruct{k,2};
    suffix   = dataStruct{k,3};

    % Get the variable dynamically
    if ~exist(varName, 'var')
        warning('Variable %s not found in workspace, skipping...', varName);
        continue;
    end
    Y = eval(varName);

    % Build full file path
    file_name = [base_filename, suffix];
    full_path = fullfile(folder_path_new, file_name);

    % Open file for writing
    fileID = fopen(full_path, 'w');
    if fileID == -1
        warning('Could not open %s for writing. Skipping...', full_path);
        continue;
    end

    % ---- Write header ----
    fprintf(fileID, 'M\t%s\n', header);

    % ---- Write data ----
    for i = 1:length(M_valid)
        fprintf(fileID, '%.8f\t%.15e\n', M_valid(i), Y(i));
    end

    % ---- Close file ----
    fclose(fileID);
    fprintf('File written: %s\n', full_path);
end


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

% VALOR DEL r0
function x3=r0
global v3
x3=v3;
end 

function setr0(val3)
global v3
v3 = val3;
end

