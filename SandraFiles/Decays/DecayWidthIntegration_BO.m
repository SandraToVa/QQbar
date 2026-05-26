%Program that reproduces the DecayWidthIntegration_QQ program in the BO
%limit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transition to compute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define parameters for the integration
DeltaE = zeros(1, 4); %3 transitions in GeV
% Our values
DeltaE(1,1) = 0.556; %4s-2s
DeltaE(1,2) = 0.876; %5s-2s
DeltaE(1,3) = 0.529; %5s-3s
DeltaE(1,4) = 0.347; %5s-3s

r_vev = zeros(1,4);
% Values form our calculations
r_vev(1,1) = -0.257764050621146; %4s-2s
r_vev(1,2) = -0.138643690152830; %5s-2s
r_vev(1,3) = -0.324767315124361; %5s-3s
r_vev(1,4) = -0.793741870480275; %3s-2s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dipion array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%value of the pion mass 
m=0.14;

%value for the dipion vector
Min = 2*m; % Double mass of the pion
k = 50; % Number of points for integration

% Preallocate a matrix to hold the 3 different M arrays
% 3 rows (for each DeltaE) and k+1 columns (for the points)
M = zeros(length(DeltaE), k + 1); 

% Loop through each DeltaE value to populate the matrix
for i = 1:length(DeltaE)
    Mfin = DeltaE(i);
    
    % linspace is highly recommended here over the colon operator (Min:nM:Mfin)
    % It automatically guarantees exactly k+1 points without floating-point errors.
    M(i, :) = linspace(Min, Mfin, k + 1);
end


file_name = 'QQbottom_4s-2s_integral_BOlimit.txt';
num=1;

M_valid = M(num, :);

I_if_valid = ( r_vev(num)./2 ).^2 ;
I_if_c_valid = ( r_vev(num)./6 ).^2 ;
I_if_0c_valid = ( r_vev(num)).^2 / 12 ;
I_if_1_valid = 0;
I_if_2_valid = 0;


Cm2_values = zeros(1, k+1);  % Preallocate for speed
CM2_values = zeros(1, k+1);  
CECM_values = zeros(1, k+1);  
CmCM_values = zeros(1, k+1);  
CmCE_values = zeros(1, k+1);  
CE2_values = zeros(1, k+1);  


for i = 1:k+1
    Cm2_values(i) = Cm2Term(I_if_valid, m, M_valid(i), DeltaE(num));
    CM2_values(i) = CM2Term(I_if_valid, m, M_valid(i), DeltaE(num));
    CECM_values(i) = CECMTerm(I_if_valid, I_if_0c_valid, m, M_valid(i), DeltaE(num));
    CmCM_values(i) = CmCMTerm(I_if_valid, m, M_valid(i), DeltaE(num));
    CmCE_values(i) = CmCETerm(I_if_valid, I_if_0c_valid, m, M_valid(i), DeltaE(num));
    CE2_values(i) = CE2Term(I_if_valid, I_if_0c_valid, I_if_c_valid, I_if_1_valid, I_if_2_valid, m, M_valid(i), DeltaE(num));
end

M2_valid = M_valid.^2;

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
%file_name = 'QQbottom_5s-3s_integral_prova.txt';
folder_path = '/Users/sandra/Documents/Doctorat/Projectes PhD/Transicions a 2 pions/Lower order Lagrangian/Quarkonium/IntegratedConstants';

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



%Remember: in the mathematica code i will multiply the decay width by
%(4 * fpi^4 pi^4) so this has to be taken out from here

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

