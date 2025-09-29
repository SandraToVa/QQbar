
% Programa per fer les integrals de la dipionmass i obtindre les constants
% per les quals estan multiplicades C_pi, C_pipi, C_e i fer los ajustos per
% a QQ->QQ transitions

%Valor r0
setr0(3.964)
%massa
load("dades.mat","m_c","m_b")
setm_q(m_c)

%Change this depending on the transition

ComputationIE = TransitionsAdded('QQStoS');      

[I_if_square_cell, DeltaE, M] = ComputationIE(3, 2, 0.28, 0.454, 50);

%Obtain the I_if^2 from the cell
I_if_square = I_if_square_cell{1,1};
I_if_c_square = I_if_square_cell{1,2};
I_if_0c_square = I_if_square_cell{1,3};
I_if_s_square = I_if_square_cell{1,4};
%I_if_square = zeros(1,51);
%I_if_c_square = zeros(1,51);
%I_if_0c_square = zeros(1,51);
%I_if_s_square = zeros(1,51);


% Value of pion mass 140MeV
m = 0.14;

% Compute M^2
M2 = M.^2;

% Filter values where DeltaE > M
valid_idx = (DeltaE > M);
M_valid = M(valid_idx);
M2_valid = M2(valid_idx);
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

% Perform numerical integration for each term
Cm2_integral = trapz(M2_valid, Cm2_values);
CM2_integral = trapz(M2_valid, CM2_values);
CECM_integral = trapz(M2_valid, CECM_values);
CmCM_integral = trapz(M2_valid, CmCM_values);
CmCE_integral = trapz(M2_valid, CmCE_values);
CE2_integral = trapz(M2_valid, CE2_values);

% Display results
fprintf('Cm^2 Integral: %f\n', Cm2_integral);
fprintf('CM^2 Integral: %f\n', CM2_integral);
fprintf('CE*CM Integral: %f\n', CECM_integral);
fprintf('Cm*CM Integral: %f\n', CmCM_integral);
fprintf('Cm*CE Integral: %f\n', CmCE_integral);
fprintf('CE^2 Integral: %f\n', CE2_integral);


% Define file path and name
% The following code is to create .txt files where the variables are stored
file_name = 'QQcharm_3s-2s_integral.txt';
folder_path = '/Users/sandra/Documents/Doctorat/Projectes PhD/Transicions a 2 pions/Lower order Lagrangian/IntegratedConstants';

full_path = fullfile(folder_path, file_name);

% Open the file for writing
fileID = fopen(full_path, 'w');

% Check if the file was opened successfully
if fileID == -1
    error('Could not open file for writing.');
end

% Write results to the file
fprintf(fileID, 'Integral Results for only I_{i->f}:\n');
fprintf(fileID, '-----------------\n');
fprintf(fileID, 'Cm^2 Integral: %.50f\n', Cm2_integral);
fprintf(fileID, 'CM^2 Integral: %.50f\n', CM2_integral);
fprintf(fileID, 'CE*CM Integral: %.50f\n', CECM_integral);
fprintf(fileID, 'Cm*CM Integral: %.50f\n', CmCM_integral);
fprintf(fileID, 'Cm*CE Integral: %.50f\n', CmCE_integral);
fprintf(fileID, 'CE^2 Integral: %.50f\n', CE2_integral);

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
        (2 * 2 * I_if_s_square * M^4 * (-4 * m^2 + M^2)^2 + I_if_square * (192 * m^4 * DeltaE^4 ...
        - 16 * m^2 * M^2 * DeltaE^4 + 7 * M^4 * DeltaE^4) ...
        + 2 * I_if_0c_square * DeltaE^2 * (4 * M^6 - 15 * M^4 * DeltaE^2 + ...
        32 * m^4 * (2 * M^2 - 15 * DeltaE^2) + m^2 * (88 * M^4 + 60 * M^2 * DeltaE^2)) ...
        + I_if_c_square * (7 * M^4 * (4 * M^4 + 5 * DeltaE^4) + 16 * m^4 * (8 * M^4 - 20 * M^2 * DeltaE^2 + 75 * DeltaE^4) ...
        + 8 * m^2 * (12 * M^6 - 50 * M^4 * DeltaE^2 - 25 * M^2 * DeltaE^4)));
    term1=(1 / (15 * M^4)) * sqrt(DeltaE^2 - M^2) * sqrt(1 - (4 * m^2) / (M^2)) * 2 * 2 * I_if_s_square * M^4 * (-4 * m^2 + M^2)^2;
    term2=(1 / (15 * M^4)) * sqrt(DeltaE^2 - M^2) * sqrt(1 - (4 * m^2) / (M^2)) *I_if_square * (192 * m^4 * DeltaE^4 ...
        - 16 * m^2 * M^2 * DeltaE^4 + 7 * M^4 * DeltaE^4);
    term3=(1 / (15 * M^4)) * sqrt(DeltaE^2 - M^2) * sqrt(1 - (4 * m^2) / (M^2)) * 2 * I_if_0c_square * DeltaE^2 * (4 * M^6 - 15 * M^4 * DeltaE^2 + ...
        32 * m^4 * (2 * M^2 - 15 * DeltaE^2) + m^2 * (88 * M^4 + 60 * M^2 * DeltaE^2));
    term4=(1 / (15 * M^4)) * sqrt(DeltaE^2 - M^2) * sqrt(1 - (4 * m^2) / (M^2)) *I_if_c_square * (7 * M^4 * (4 * M^4 + 5 * DeltaE^4) + 16 * m^4 * (8 * M^4 - 20 * M^2 * DeltaE^2 + 75 * DeltaE^4) ...
        + 8 * m^2 * (12 * M^6 - 50 * M^4 * DeltaE^2 - 25 * M^2 * DeltaE^4));
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

