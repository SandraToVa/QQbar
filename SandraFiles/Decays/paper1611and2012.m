%Numerical integral of the decay width from papers 1611.00913v2, 2012.05034v2
%Mathematica document: PhD_ComparisonLambdaEta.nb

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transition to compute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define parameters for the integration
DeltaE = zeros(1, 3); %3 transitions in GeV
% Values form pdg
DeltaE(1,1) = 0.556; %4s-2s
DeltaE(1,2) = 0.8618; %5s-2s
DeltaE(1,3) = 0.5301; %5s-3s

r_vev = zeros(1,3);
% Values form our calculations
r_vev(1,1) = -0.257764050621146; %4s-2s
r_vev(1,2) = -0.138643690152830; %5s-2s
r_vev(1,3) = -0.324767315124361; %5s-3s


% Constants for each transition
% c1 and c2 for the 3 transitions in GeV^-1
c1 = zeros(1, 3);
c2 = zeros(1, 3);
%4s-2s
c1(1,1) = 1.2 / 10;
c2(1,1) = -1.0 / 10;
%5s-2s
c1(1,2) = 21.1 / 10000;
c2(1,2) = -12.6 / 10000;
%5s-3s
c1(1,3) = -17.8 / 10000;
c2(1,3) = 16.1 / 10000;


% Constants for each transition our
% c1 and c2 for the 3 transitions in GeV^-1
cm = zeros(1, 3);
cM = zeros(1, 3);
ce = zeros(1, 3);
%4s-2s
cm(1,1) = 0.275051;
cM(1,1) = -0.653245;
ce(1,1) = 0.51572;
%5s-2s
cm(1,2) = 0.0217332;
cM(1,2) = -0.0229477;
ce(1,2) = 0.0120811;
%5s-3s
cm(1,3) = -0.00185559;
cM(1,3) = 0.00751786;
ce(1,3) = -0.00659007;



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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Invariant mass arrey
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate dGamma to perfectly match the size of M (3 rows, k+1 columns)
dGamma_their = zeros(size(M)); 
dGamma_our = zeros(size(M));

% Get the number of transitions (rows) and points (columns)
num_transitions = size(M, 1);
num_points = size(M, 2); % This is k + 1

% Outer loop: Iterate through each of the 3 transitions
for i = 1:num_transitions
    
    % Inner loop: Iterate through all 51 mass points for the current transition
    for j = 1:num_points
        
        % Calculate dGamma using the specific DeltaE and M value
        dGamma_their(i, j) = dw(m, M(i, j), DeltaE(i), c1(i), c2(i));

        %All our functions
        R2 = ( r_vev(i)./2 ).^2;
        R02 = ( r_vev(i)./6 ).^2;
        RR0 = ( r_vev(i)).^2 / 12;

        Cm2 = Cm2Term(R2, m, M(i, j), DeltaE(i));
        CM2 = CM2Term(R2, m, M(i, j), DeltaE(i));
        CECM = CECMTerm(R2, RR0, m, M(i, j), DeltaE(i));
        CmCM = CmCMTerm(R2, m, M(i, j), DeltaE(i));
        CmCE = CmCETerm(R2, RR0, m, M(i, j), DeltaE(i));
        CE2 = CE2Term(R2, RR0, R02, 0, 0, m, M(i, j), DeltaE(i));

        FracFpi = 1/(2^2 * pi^4 * 0.092^4);

        dGamma_our(i, j) = FracFpi * (ce(i)^2 * CE2 + ce(i)*cm(i) * CmCE + cm(i)*cM(i) * CmCM + cM(i)^2 * CM2 + cm(i)^2 * Cm2 + ce(i)*cM(i)* CECM);

        
    end
end

% Now dGamma is a 3x51 matrix holding the results for all transitions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decay width integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute M^2
M2 = M.^2;

Gamma_their=zeros(1,3);
Gamma_our=zeros(1,3);

for i = 1:num_transitions

    Gamma_their(i) = trapz(M2(i,:), dGamma_their(i,:));
    Gamma_our(i) = trapz(M2(i,:), dGamma_our(i,:));

end

% Display results
fprintf('Their results \n');
fprintf('Gamma 4s-2s (GeV): %.15f\n', Gamma_their(1));
fprintf('Gamma 5s-2 (GeV): %.15f\n', Gamma_their(2));
fprintf('Gamma 5s-3 (GeV): %.15f\n', Gamma_their(3));
%
fprintf('Our results \n');
% Display results
fprintf('Gamma 4s-2s (GeV): %.15f\n', Gamma_our(1));
fprintf('Gamma 5s-2 (GeV): %.15f\n', Gamma_our(2));
fprintf('Gamma 5s-3 (GeV): %.15f\n', Gamma_our(3));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Their

function dw_c1c2 = dw(m, M, DeltaE, c1, c2)
    
    fpi=0.092;
    const = 1 / (480 * fpi^4 * pi^3 * M^4 );
    sqrts = sqrt(DeltaE^2 - M^2) * sqrt( 1 - 4*m^2/M^2 );

    term1M4 = 3 * c2^2 * ( M^2 - 4*m^2 )^2;
    term2M4 = 60 * c1^2 * ( M^2 - 2*m^2 )^2;
    term3M4 = 20 * c1 * c2 * ( 8*m^4 - 6*m^2*M^2 + M^4 );

    term1M2 = -8 * ( 5*c1 + 3*c2 ) * m^4;
    term2M2 = 2 * c2 * m^2 * M^2;
    term3M2 = M^4 * ( 10*c1 + c2 );

    termD4 = 6*m^4 + 2*m^2*M^2 + M^4;

    dw_c1c2 = const * sqrts * ( M^4 * ( term1M4 + term2M4 + term3M4 ) + ...
        4 * c2 * M^2 * ( term1M2 + term2M2 + term3M2  ) * DeltaE^2 + ...
        8 * c2^2 * DeltaE^4 * termD4 );

end


%Ours
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
