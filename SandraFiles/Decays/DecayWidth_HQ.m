%Program to plot the decay width with respect to the dipion mass for each
%transition and also for each set of coeficient values per a HQ->QQ
%transitions

%The coeficients have been computed separately with
%PhDSolvingDecayConstants.nb only one coificient is used here

%The body of the program is the same of the others. I compute for each form
%factor the expected value of the angual integral (in Hibrid_ItoF)
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
setm_q(m_c)

%Change this depending on the transition

ComputationIE = TransitionsAdded('HQsdtP');      

[I_if_square_cell, DeltaE, M] = ComputationIE(3, 2, 0.28, 0.387, 50);

%%

file_name = 'GammaVSmass_HQbottom_2sd1-2p.txt';

% Value of pion mass 140MeV
m = 0.14;
% Value of the global constant that divides everithing in GeV. Fpi=92MeV
FracFpi = 1/(pi * 0.092^4 * 0.187) * 4^2;
%0.187GeV^2 is the string tension
%4^2 per compensar que Ce=eta/4 i tenim eta^2

%Obtain the I_if^2 from the cell
I_if_01_square = I_if_square_cell{1, 1};
I_if_11_square = I_if_square_cell{1, 2};
I_if_0N_square = I_if_square_cell{1, 3};
I_if_NN_square = I_if_square_cell{1, 4};
I_if_10_square = I_if_square_cell{1, 5};
I_if_N0_square = I_if_square_cell{1, 6};
I_if_0N_10 = I_if_square_cell{1, 7};
I_if_01_N0 = I_if_square_cell{1, 8};
I_if_11_NN = I_if_square_cell{1, 9};

% Filter values where DeltaE > M
valid_idx = (DeltaE > M);
M_valid = M(valid_idx);

I_if_01_square = I_if_01_square(valid_idx);
I_if_11_square = I_if_11_square(valid_idx);
I_if_0N_square = I_if_0N_square(valid_idx);
I_if_NN_square = I_if_NN_square(valid_idx);
I_if_10_square = I_if_10_square(valid_idx);
I_if_N0_square = I_if_N0_square(valid_idx);
I_if_0N_10 = I_if_0N_10(valid_idx);
I_if_01_N0 = I_if_01_N0(valid_idx);
I_if_11_NN = I_if_11_NN(valid_idx);


num_valid = length(M_valid);
DW_values = zeros(1, num_valid);


for i = 1:num_valid
    term1 = HQdecaywidthTerm1(I_if_01_square(i), I_if_0N_square(i), I_if_0N_10(i), ...
    I_if_01_N0(i), I_if_10_square(i), I_if_N0_square(i), m, M_valid(i), DeltaE) ;

    term2 = HQdecaywidthTerm2(I_if_11_square(i), I_if_NN_square(i), I_if_11_NN(i), m, M_valid(i), DeltaE);

    DW_values(i) = term1 + term2;
end

% Now with each part of the gamma separated by parameters we define the
% parameters computed with mathematica:
% Each row is a set: [CE]  -> corresponds to [Ce]
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
    DW = DW_values(i);
    
    for N = 1:4
        CE = param_sets(N,3);
        
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




function term = HQdecaywidthTerm1(I_if_01_square, I_if_0N_square, I_if_0N_10, ...
    I_if_01_N0, I_if_10_square, I_if_N0_square, m, M, DeltaE)

    term = (1 / (240 * M^4)) * sqrt(DeltaE^2 - M^2)^3 * sqrt(1 - (4 * m^2) / (M^2)) ...
        * (M^2 * (M^2 - 4*m^2)^2 ) * ( 2*I_if_01_square + 2*I_if_0N_square - 2*sqrt(2)*I_if_0N_10 ...
        + 2*sqrt(2)*I_if_01_N0 + I_if_10_square + I_if_N0_square ) ;

end

function term = HQdecaywidthTerm2(I_if_11_square, I_if_NN_square, I_if_11_NN, m, M, DeltaE)

    term = (1 / (240 * M^4)) * sqrt(DeltaE^2 - M^2)^3 * sqrt(1 - (4 * m^2) / (M^2)) ...
        * (8 *(I_if_11_square + I_if_NN_square + 2*I_if_11_NN) * (6*m^4 + 2*m^2*M^2 + M^4) * DeltaE^2 ) ;

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

