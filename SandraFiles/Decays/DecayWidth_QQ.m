
%Program to plot the decay width with respect to the dipion mass for each
%transition and also for each set of coeficient values per a QQ->QQ
%transitions

%The coeficients have been computed separately with
%PhDSolvingDecayConstants.nb

%The body of the program is the same of the others. I compute for each form
%factor the expected value of the angual integral (in ForFactor_ItoF)
%sandwiched with the radial part of the wave function normalized (using the
%ExpValFunctions in a way for each transition in TransitionsAdded). Example
%of use in ComputationExpVal where the FormFactors vs M are saved in a txt

%Now we do the same but instead of only computing the Form factors we
%multiply them for all the terms and constants to get the decay width and
%then add them up for each dipion mass and store them in a txt


% We work in GeV!!

%Valor r0
setr0(3.964)
%massa
load("dades.mat","m_c","m_b")
setm_q(m_b)

%Change this depending on the transition

ComputationIE = TransitionsAdded('QQStoS');      

[I_if_square_cell, DeltaE, M] = ComputationIE(5, 3, 0.28, 0.55, 50);

%Obtain the I_if^2 from the cell
I_if_square = I_if_square_cell{1,1};
I_if_c_square = I_if_square_cell{1,2};
I_if_0c_square = I_if_square_cell{1,3};
I_if_s_square = I_if_square_cell{1,4};


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
I_if_s_valid = I_if_s_square(valid_idx);


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
    CE2_values(i) = CE2Term(I_if_valid(i), I_if_0c_valid(i), I_if_c_valid(i), I_if_s_valid(i), m, M_valid(i), DeltaE);
end

% Now with each part of the gamma separated by parameters we define the
% parameters computed with mathematica:
% Each row is a set: [Cm, CM, CE]  -> corresponds to [Cp, Cpp, Ce]
param_sets = [
    -0.523739,  0.00249372,   0.0534711;
    -0.516711,  -0.0608308,   0.115396;
     0.516711,  0.0608308,    -0.115396;
     0.523739,  -0.00249372,  -0.0534711
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
file_name = 'GammaVSmass_QQbottom_5s-3s.txt';
folder_path = '/Users/sandra/Documents/Doctorat/Projectes PhD/Transicions a 2 pions/Lower order Lagrangian/DecayWidth';

full_path = fullfile(folder_path, file_name);

% Open the file for writing
fileID = fopen(full_path, 'w');

% Check if the file was opened successfully
if fileID == -1
    error('Could not open file for writing.');
end

% Write header
fprintf(fileID, 'Decay Width results for the different constant sets:\n');
fprintf(fileID, 'Mass (GeV)     Gamma1 (GeV)       Gamma2 (GeV)       Gamma3 (GeV)       Gamma4 (GeV)\n');
fprintf(fileID, '----------------------------------------------------------------------------------------\n');

for i = 1:num_valid
    fprintf(fileID,    '%.4f %.15f %.15f %.15f %.15f\n', M_valid(i), Gamma(1,i), Gamma(2,i), Gamma(3,i), Gamma(4,i));
end

% Close the file
fclose(fileID);

fprintf('Results saved to %s\n', full_path);




function term = Cm2Term(I_if_square, m, M, DeltaE)
    term = (1 / (15 * M^4)) * sqrt(DeltaE^2 - M^2) * sqrt(1 - (4 * m^2) / (M^2)) ...
        * (15 * I_if_square * m^4 * M^4);
end

function term = CM2Term(I_if_square, m, M, DeltaE)
    term = (1 / (15 * M^4)) * sqrt(DeltaE^2 - M^2) * sqrt(1 - (4 * m^2) / (M^2)) ...
        * (15 * I_if_square * M^8);
end

function term = CECMTerm(I_if_square, I_if_0c_square, m, M, DeltaE)
    term = (1 / (15 * M^4)) * sqrt(DeltaE^2 - M^2) * sqrt(1 - (4 * m^2) / (M^2)) ...
        * (10 * M^4 * (I_if_square * (8 * m^2 + M^2) * DeltaE^2 ...
        + I_if_0c_square * (4 * M^4 - M^2 * DeltaE^2 + 4 * m^2 * (2 * M^2 - 5 * DeltaE^2))));
end

function term = CmCMTerm(I_if_square, m, M, DeltaE)
    term = (1 / (15 * M^4)) * sqrt(DeltaE^2 - M^2) * sqrt(1 - (4 * m^2) / (M^2)) ...
        * (30 * I_if_square * m^2 * M^6);
end

function term = CmCETerm(I_if_square, I_if_0c_square, m, M, DeltaE)
    term = (1 / (15 * M^4)) * sqrt(DeltaE^2 - M^2) * sqrt(1 - (4 * m^2) / (M^2)) ...
        * (10 * m^2 * M^2 * (I_if_square * (8 * m^2 + M^2) * DeltaE^2 ...
        + I_if_0c_square * (4 * M^4 - M^2 * DeltaE^2 + 4 * m^2 * (2 * M^2 - 5 * DeltaE^2))));
end

function term = CE2Term(I_if_square, I_if_0c_square, I_if_c_square, I_if_s_square, m, M, DeltaE)
    term = (1 / (15 * M^4)) * sqrt(DeltaE^2 - M^2) * sqrt(1 - (4 * m^2) / (M^2)) * ...
        (2*2 * I_if_s_square * M^4 * (-4 * m^2 + M^2)^2 + I_if_square * (192 * m^4 * DeltaE^4 ...
        - 16 * m^2 * M^2 * DeltaE^4 + 7 * M^4 * DeltaE^4) ...
        + 2 * I_if_0c_square * DeltaE^2 * (4 * M^6 - 15 * M^4 * DeltaE^2 + ...
        32 * m^4 * (2 * M^2 - 15 * DeltaE^2) + m^2 * (88 * M^4 + 60 * M^2 * DeltaE^2)) ...
        + I_if_c_square * (7 * M^4 * (4 * M^4 + 5 * DeltaE^4) + 16 * m^4 * (8 * M^4 - 20 * M^2 * DeltaE^2 + 75 * DeltaE^4) ...
        + 8 * m^2 * (12 * M^6 - 50 * M^4 * DeltaE^2 - 25 * M^2 * DeltaE^4)));
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

