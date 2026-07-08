%Functions to compute spectrum and w.f. of quarkonium with Jcal=1
% Quarkonium S=1, Jcal=1 they can be several states:
% - S-States J=1, L=0
% - P-states J=1, L=1
% - D-states J=1, L=2
% Autor: Sandra Tomàs
% Creador subrutines i font: Ruben Oncala

function [E,W,x]=QQbarS1J1(m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh and potential matrix                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % the endpoints of the integration interval:
system.a=0.001;   
system.b=22; 
% parameters of the boundary conditions:

%---------Coose n!----------
% without mixing just quarkonium spectrum codes
% with mixing of hybrids spin 0 and quarkonium spin 1: n=6
n=6;
system.A1= eye(n);
system.A2= zeros(n,n);
system.B1= eye(n);
system.B2= zeros(n,n);
% function handle to the function returning the potential matrix (capture mass m)
system.V=@(x) potentialMatrix(x,m);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of states we compute               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kmax=6; 
%EigvData.eigenvalues returens 3 times the same state. We only need it 1

%number of actual states
E=zeros(1,kmax+1);
tol0=5e-6; 
[EigvData,meshData]=computeEigenvalues(system,0,kmax,tol0);
E(:) = EigvData.eigenvalues(:) / m;

for i=1:length(E)
    [x,Y,~]=computeEigenfunction(system,meshData,E(i)*m,1);
    W(:,:,i)=Y;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System parameters                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TOTAL ANGULAR MOMENTUM (Jcal=J for quarkonium)
function [j]=parameters1
j=1;
end 


% Spectrum true value
% The potentias as they are they compute E(GeV) to obtain the spectrum we
% need to add 2mQ + Eg
function [Eg]=parameters2
Eg=false;
end

% Mixing parameter
function [mix]=parameters3
mix=true;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Potential Matrix                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function r=potentialMatrix(x,m) 
% returns the potential matrix evaluated in x

[j]=parameters1;
[mix]=parameters3;

if mix == true
    r = zeros(6,6,4);

    for i=1:4
        % Get diagonal and off-diagonal components
        v11 = v0Diagonal(x(i), j, m, 'v11');
        v22 = v0Diagonal(x(i), j, m, 'v22');
        v33 = v0Diagonal(x(i), j, m, 'v33');
        v44 = v0Diagonal(x(i), j, m, 'v44');
        v55 = v0Diagonal(x(i), j, m, 'v55');
        v66 = v0Diagonal(x(i), j, m, 'v66');
    
        [v12, v35, v36, v45, v46] = v0OffDiagonal(x(i), j, m);
        v56 = v0Coupling(x(i), j, m);
    
    
        % Assembly with safety checks
        r(1,1,i) = v11;
        r(1,2,i) = v12; 
        r(2,1,i) = v12;
    
        r(2,2,i) = v22;

        r(3,3,i) = v33;
        r(3,5,i) = v35;
        r(5,3,i) = v35;
        r(3,6,i) = v36;
        r(6,3,i) = v36;
    
        r(4,4,i) = v44;
        r(4,5,i) = v45;
        r(5,4,i) = v45;
        r(4,6,i) = v46;
        r(6,4,i) = v46;

        r(5,5,i) = v55;
        r(5,6,i) = v56;
        r(6,5,i) = v56;

        r(6,6,i) = v66;
    
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Potential functions                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = VSigG(r,m)
    %We define if we want to compute the true spectrum or just the energy E
    [Eg]=parameters2;

    % Perturvation part
    V0 = V0Pert_Calc(r);
    
    % Fit potential VSigG
    num_term = 0.724838892832331 - 2.5561383325357028 * r.^3 + ...
               8.017804077992961 * r.^5 + 41.54269365639395 * ...
               (-0.004879614124924147 - 0.16318790345716533 ./ r.^3 - pi ./ (12 * r) + 0.21 * r) .* r.^9 ...
               + V0;  
               
    den_term = 1 + 41.54269365639395 * r.^9;
    
    y = 0.004879614124924147 + (num_term ./ den_term);

    if Eg == true
        if m == 1.496
            y = y + 2*m - 0.44;
        else
            y = y + 2*m + 0;
        end
    end
end
function y = VPiU(r,m)
    %We define if we want to compute the true spectrum or just the energy E
    [Eg]=parameters2;

    % Perturvation part
    V0 = V0Pert_Calc(r);
    
    % Fit potential VPiU
    num_term = 1.1142023723639383 + 1.0032809918713739 * r.^3 + ...
               0.14416332343387012 * r.^5 + 0.018733976201780835 * r.^9 .* ...
               (-0.007999908701569785 + sqrt(1.2095131716320702 + 0.0441 * r.^2)) ...
               - (V0 / 8); 
               
    den_term = 1 + 0.555086331387556 * r.^3 + 0.33981203544081023 * r.^4 + ...
               0.018733976201780835 * r.^9;
               
    y = 0.007999908701569785 + (num_term ./ den_term);

    if Eg == true
        if m == 1.496
            y = y + 2*m - 0.44;
        else
            y = y + 2*m + 0;
        end
    end
end
function y = VSigU(r,m)
    %We define if we want to compute the true spectrum or just the energy E
    [Eg]=parameters2;

    % Perturvation part
    V0 = V0Pert_Calc(r);
    
    % Fit potential VPiU
    num_term = 1.1131475711643908 + 1.2576502448714892 * r.^3 + ...
               1.5668236020705248 * r.^4 + 0.34642902531803854 * r.^5 + ...
               0.0000348777989761813 * r.^9 .* ...
               (-0.007999908701569785 + sqrt(3.848451000647496 + 0.0441 * r.^2)) ...
               - (V0 / 8);  
               
    den_term = 1 + 0.6407603775634713 * r.^3 + 1.7884530121585636 * r.^4 + ...
               0.0000348777989761813 * r.^9;
               
    y = 0.007999908701569785 + (num_term ./ den_term);

    if Eg == true
        if m == 1.496
            y = y + 2*m - 0.44;
        else
            y = y + 2*m + 0;
        end
    end
end
function y=Vq(r,m)
    y=VPiU(r,m)-VSigU(r,m);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mixing functions                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = VPiMix(r, m)
    % Constant
    if m == 1.496
        cF = 1.12;
    else
        cF = 0;
    end
    
    num = 0.253 + 0.00841297265408362 * r.^2 + 4.980786383914752e-6 * r.^7;
    den = 1 + 0.3513451629416728 * r.^2 + 8.19809128976659e-6 * r.^9;
    
    y = (cF / m ) * (num ./ den);
end
function y = VSigMix(r, m)
    % Constant
    if m == 1.496
        cF = 1.12;
    else
        cF = 0;
    end
    
    num = 0.253 + 0.14378213128924067 * r.^2 + 0.00018940820727945146 * r.^6;
    den = 1 + 1.19191862857913 * r.^2 - 0.29186710310745323 * r.^3 + ...
          0.13264734605500272 * r.^4 + 0.000020150616940783642 * r.^9;
    
    y = (cF / m ) * (num ./ den);
end
function y = VMixq(r, m)
    y = VSigMix(r, m) - VPiMix(r, m);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliar functions                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function V0Pert_val = V0Pert_Calc(r)
    % V0Pert_Calc - This computes the perturvative part shared by all
    % potentials
    
    % --- Global Constants ---
    EulerGamma = 0.5772156649015328;
    zeta3 = 1.202056903159594; 
    zeta5 = 1.036927755143370; 
    CA = 3;
    CF = 4/3;
    
    beta0_pot = 9;
    beta1_pot = 64;
    beta2_pot = 3863/6;
    
    a1 = 7;
    a2 = 100/9 - 2*(55/3 - 16*zeta3) + 9*(4343/162 + (16*pi^2 - pi^4)/4 + (22*zeta3)/3) - (9*(1798/81 + (56*zeta3)/3))/2;
    a3 = 5199.767568902478;
    
    Nm = 0.5626;
    b = 32/81;
    c0 = 1;
    c1p = -3397/20736;
    
    c2p = (-322289305/4 + (6561 * (1349963/54 + 3564 * zeta3 + ...
           9 * (50065/162 + (6472 * zeta3) / 81) - ...
           3 * (1078361/162 + (6508 * zeta3) / 27))) / 128) / 329204736;
           
    c3 = (-39200051773/4251528 + ...
          (3397 * (-479552/729 + (1349963/54 + 3564 * zeta3 + ...
               9 * (50065/162 + (6472 * zeta3) / 81) - 3 * (1078361/162 + ...
                 (6508 * zeta3) / 27)) / 2304)) / 3 + ...
          648 * (-83397479/26244 + (1349963/54 + 3564 * zeta3 + ...
              9 * (50065/162 + (6472 * zeta3) / 81) - ...
              3 * (1078361/162 + (6508 * zeta3) / 27)) / 162 + ...
            (-8157455/16 + (9801 * pi^4)/20 - 81 * (1205/2916 - (152 * zeta3) / 81) - ...
              (621885 * zeta3) / 2 - 9 * (25960913/1944 - (5263 * pi^4) / 405 + ...
                (698531 * zeta3) / 81 - (381760 * zeta5) / 81) + 288090 * zeta5 - ...
              27 * (-630559/5832 + (809 * pi^4) / 1215 - (48722 * zeta3) / 243 + ...
                (460 * zeta5) / 9) - 3 * (-336460813/1944 + (6787 * pi^4) / 108 - ...
                (4811164 * zeta3) / 81 + (1358995 * zeta5) / 27)) / 9216)) / 78274560;

    S = @(n) c0 * (gamma(n + 1 + b) / gamma(1 + b)) + ...
             c1p * (gamma(n + 1 + b - 1) / gamma(1 + b - 1)) + ...
             c2p * (gamma(n + 1 + b - 2) / gamma(1 + b - 2)) + ...
             c3 * (gamma(n + 1 + b - 2) / gamma(1 + b - 2));

    % --- Evaluation ---
    nuval = max(1 ./ r, 1);
    as_nuval = alphaSAnalytic(nuval, zeta3);
    nuusval = max(3 * as_nuval ./ r, 1);
    
    log_term = log(nuval .* exp(EulerGamma) .* r);
    
    anu1 = a1 + 2 * beta0_pot * log_term;
    anu2 = a2 + (pi^2 / 3) * beta0_pot^2 + (4 * a1 * beta0_pot + 2 * beta1_pot) * log_term + 4 * beta0_pot^2 * log_term.^2;
    anu3 = a3 + a1 * beta0_pot^2 * pi^2 + (5 * pi^2 / 6) * beta0_pot * beta1_pot + 13 * zeta3 * beta0_pot^3 + ...
           (2 * pi^2 * beta0_pot^3 + 6 * a2 * beta0_pot + 4 * a1 * beta1_pot + 2 * beta2_pot + (16/3) * CA^3 * pi^2) * log_term + ...
           (12 * a1 * beta0_pot^2 + 10 * beta0_pot * beta1_pot) * log_term.^2 + 8 * beta0_pot^3 * log_term.^3;
           
    delta_a3us = (16/3) * CA^3 * pi^2 * log(nuusval ./ nuval);
    
    d0 = (beta0_pot / 2) * log(nuval);
    d1 = (beta1_pot / 8) * log(nuval);
    
    delta_m2 = Nm * (beta0_pot / (2*pi)) * (S(1) * (2*d0 / pi) + (beta0_pot / (2*pi)) * S(2));
    delta_m3 = Nm * (beta0_pot / (2*pi)) * (S(1) * ((3*d0.^2 + 2*d1) / pi^2) + ...
               (beta0_pot / (2*pi)) * S(2) * (3*d0 / pi) + (beta0_pot / (2*pi))^2 * S(3));
               
    delta_mRSp = as_nuval.^2 * 1.207713407019128 + as_nuval.^3 .* delta_m2 + as_nuval.^4 .* delta_m3;
    
    bracket = 1 + (as_nuval / (4*pi)) .* anu1 + (as_nuval / (4*pi)).^2 .* anu2 + (as_nuval / (4*pi)).^3 .* (anu3 + delta_a3us);
    Vs0_val = -(CF * as_nuval ./ r) .* bracket + 2 * delta_mRSp;
    
    V0Pert_val = real(Vs0_val);
end

% --- Local Functions ---
function vals = alphaSAnalytic(mu, zeta3)
    mc = 1.27; mb = 4.18; mt = 175;
    lambda3 = 0.32920945283353376; lambda4 = 0.28904534271282883;
    lambda5 = 0.20834647918878832; lambda6 = 0.08768634615464756;
    
    vals = zeros(size(mu));
    mask1 = mu < mc; if any(mask1), vals(mask1) = alphaS4loop(mu(mask1), lambda3, 3, zeta3); end
    mask2 = (mu >= mc) & (mu < mb); if any(mask2), vals(mask2) = alphaS4loop(mu(mask2), lambda4, 4, zeta3); end
    mask3 = (mu >= mb) & (mu < mt); if any(mask3), vals(mask3) = alphaS4loop(mu(mask3), lambda5, 5, zeta3); end
    mask4 = mu >= mt; if any(mask4), vals(mask4) = alphaS4loop(mu(mask4), lambda6, 6, zeta3); end
end

function res = alphaS4loop(mu, lambda, nf, zeta3)
    b0 = 11 - (2/3)*nf; b1 = 102 - (38/3)*nf;
    b2 = 2857/2 - (5033/18)*nf + (325/54)*nf^2;
    b3 = (149753/6 + 3564*zeta3) - (1078361/162 + (6508/27)*zeta3)*nf + ...
         (50065/162 + (6472/81)*zeta3)*nf^2 + (1093/729)*nf^3;
    L = log(mu.^2 / lambda^2);
    res = ((4*pi) ./ (b0 * L)) .* (1 - (b1 * log(L)) ./ (b0^2 * L) + ...
          (b1^2 * (log(L).^2 - log(L) - 1) + b0 * b2) ./ (b0^4 * L.^2) - ...
          (b1^3 * (log(L).^3 - 2.5*log(L).^2 - 2*log(L) + 0.5) + 3*b0*b1*b2*log(L) - 0.5*b0^2*b3) ./ (b0^6 * L.^3));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Potential matrix functions                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Quarkonium potentials %%%%%%%%%%%%
function y=v(x,j,m)
    y=(j*(j+1))/(x.^2)+m*VSigG(x,m);
end

% Mixing potentials %%%%%%%%%%%%
% The ones from the program of matrices of Ruben but with this m factor
% extra. Important point
% This potentials are different thet for the case of mixing spin 0
% quarkonium with spin 1 hybrids
function val = v0Diagonal(x, j, m, type)
    % Calculates specific diagonal entries for MixingH0 with element-wise ops
    switch type
        case 'v11', val = (j*(j+1))./(x.^2) + m.*VSigG(x,m);
        case 'v22', val = (j*(j+1))./(x.^2) + m.*VPiU(x,m);
        case 'v33', val = ((j+1)*(j+2))./(x.^2) + m.*VSigG(x,m);
        case 'v44', val = ((j-1)*j)./(x.^2) + m.*VSigG(x,m);
        case 'v55', val = ((j+1)*(j+2))./(x.^2) + m.*VSigU(x,m) + m.*Vq(x,m).*j./(2*j+1);
        case 'v66', val = ((j-1)*j)./(x.^2) + m.*VSigU(x,m) + m.*Vq(x,m).*(j+1)./(2*j+1);
    end
end

function [v12, v35, v36, v45, v46] = v0OffDiagonal(x, j, m)
    % Calculations for off-diagonal mixing terms
    v12 = 2*m.*VSigMix(x,m); 
    
    v35 = 2*m.*(VPiMix(x,m) + ((j+1)/(2*j+1)).*VMixq(x,m));
    v36 = -2*m.*VMixq(x,m).*sqrt(j*(j+1))/(2*j+1);
    
    v45 = v36; 
    v46 = 2*m.*(VPiMix(x,m) + (j/(2*j+1)).*VMixq(x,m));
end

function val = v0Coupling(x, j, m)
    % Calculates the v56 interaction coupling 
    val = m.*Vq(x,m).*sqrt(j*(j+1))/(2*j+1);
end

