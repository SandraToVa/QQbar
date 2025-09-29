%Special sandwitches for the trensition of 1p1 to 2s using the full
%expression with g4+g5 terms

%This file is equibalent to the hibrid_ItoF and Transitions added
%alltogeder only for thet special transiction

% PART OF THE CALL FOR OTHER PROGRAMMS

%Create a function to call the last function to compute the I_cells in DecawWidth_P1toS:
function res = TransitionsP1toS(x)
    disp(['Input value: ', x]);

    if strcmp(x, 'HQp1tS_full')
        res = @P1toStransFULL;

    end
    disp(['Output function handle: ', func2str(res)]);
end


% PART OF THE SANDWITCH
%Second we directly use the functions and make the sandwitch with the pre
%made HQ sandwitch method

%Special case of transition 1p1 to 2s where we take into acount the g5 term
%contribution

% 1p1 -> 2s (SPECIAL CASE)
function [I_Cos_cell, I_Sin_cell, I_Si_cell, DeltaE, M] = P1toStransFULL(Nin, Nfin, Min, Mfin, length)

k = length;

% Dipion mass row
nM = (Mfin - Min) / k;
M = Min:nM:Mfin;

% Initialize variables
I_Cos_cell = cell(1, 1);
I_Sin_cell = cell(1, 3);
I_Si_cell = cell(1, 4);

%Sandwitches that we need
I_Cos_R3 = zeros(1, k + 1);
I_Sin_R2 = zeros(1, k + 1);
I_Sin_R3 = zeros(1, k + 1);
I_Sin_R4 = zeros(1, k + 1);
I_Si_R1 = zeros(1, k + 1);
I_Si_R2 = zeros(1, k + 1);
I_Si_R3 = zeros(1, k + 1);
I_Si_R4 = zeros(1, k + 1);


E = zeros(1, k + 1);

% Precompute ExpValFunc

ExpValFunc = ExpValFunctions('HQ');


% Compute first E value and check if consistent
[I_Cos_R3(1), E(1)] = ExpValFunc(Nin, Nfin, M(1), @Hp1toQstrans_cos3, 1, 0, false);
[I_Sin_R2(1), ~] = ExpValFunc(Nin, Nfin, M(1), @Hp1toQstrans_sin2, 1, 0, false);
[I_Sin_R3(1), ~] = ExpValFunc(Nin, Nfin, M(1), @Hp1toQstrans_sin3, 1, 0, false);
[I_Sin_R4(1), ~] = ExpValFunc(Nin, Nfin, M(1), @Hp1toQstrans_sin4, 1, 0, false);
[I_Si_R1(1), ~] = ExpValFunc(Nin, Nfin, M(1), @Hp1toQstrans_si1, 1, 0, false);
[I_Si_R2(1), ~] = ExpValFunc(Nin, Nfin, M(1), @Hp1toQstrans_si2, 1, 0, false);
[I_Si_R3(1), ~] = ExpValFunc(Nin, Nfin, M(1), @Hp1toQstrans_si3, 1, 0, false);
[I_Si_R4(1), ~] = ExpValFunc(Nin, Nfin, M(1), @Hp1toQstrans_si4, 1, 0, false);

% Fill the rest of the arrays
for n = 2:k+1
    [I_Cos_R3(n), E(n)] = ExpValFunc(Nin, Nfin, M(n), @Hp1toQstrans_cos3, 1, 0, false);
    [I_Sin_R2(n), ~] = ExpValFunc(Nin, Nfin, M(n), @Hp1toQstrans_sin2, 1, 0, false);
    [I_Sin_R3(n), ~] = ExpValFunc(Nin, Nfin, M(n), @Hp1toQstrans_sin3, 1, 0, false);
    [I_Sin_R4(n), ~] = ExpValFunc(Nin, Nfin, M(n), @Hp1toQstrans_sin4, 1, 0, false);
    [I_Si_R1(n), ~] = ExpValFunc(Nin, Nfin, M(n), @Hp1toQstrans_si1, 1, 0, false);
    [I_Si_R2(n), ~] = ExpValFunc(Nin, Nfin, M(n), @Hp1toQstrans_si2, 1, 0, false);
    [I_Si_R3(n), ~] = ExpValFunc(Nin, Nfin, M(n), @Hp1toQstrans_si3, 1, 0, false);
    [I_Si_R4(n), ~] = ExpValFunc(Nin, Nfin, M(n), @Hp1toQstrans_si4, 1, 0, false);

    % Check if E is constant
    if E(n) ~= E(1)
        error('E values are not consistent across iterations.');
    end
end


% Store in cell array
I_Cos_cell{1, 1} = I_Cos_R3;
I_Sin_cell{1, 1} = I_Sin_R2;
I_Sin_cell{1, 2} = I_Sin_R3;
I_Sin_cell{1, 3} = I_Sin_R4;
I_Si_cell{1, 1} = I_Si_R1;
I_Si_cell{1, 2} = I_Si_R2;
I_Si_cell{1, 3} = I_Si_R3;
I_Si_cell{1, 4} = I_Si_R4;

% Return E as a single scalar if consistent
DeltaE = E(1);

end

% PART OF THE FUNCTIONS
%We write the functions that will be sandwitched

%Trans p1 -> s : cos/r^3
function Hp1s_cos3=Hp1toQstrans_cos3(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r_vec=x;

vector = cos(0.5 .* r_vec .* Sq) ./ (r_vec.^3) ;

matrix = diag(vector);

Hp1s_cos3 = matrix;

end

%Trans p1 -> s : sin/r^2
function Hp1s_sin2=Hp1toQstrans_sin2(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r_vec=x;

vector = sin(0.5 .* r_vec .* Sq) ./ (r_vec.^2) ;

matrix = diag(vector);

Hp1s_sin2 = matrix;

end

%Trans p1 -> s : sin/r^3
function Hp1s_sin3=Hp1toQstrans_sin3(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r_vec=x;

vector = sin(0.5 .* r_vec .* Sq) ./ (r_vec.^3) ;

matrix = diag(vector);

Hp1s_sin3 = matrix;

end

%Trans p1 -> s : sin/r^4
function Hp1s_sin4=Hp1toQstrans_sin4(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r_vec=x;

vector = sin(0.5 .* r_vec .* Sq) ./ (r_vec.^4) ;

matrix = diag(vector);

Hp1s_sin4 = matrix;

end

%Trans p1 -> s : (si - si)/r
function Hp1s_si1=Hp1toQstrans_si1(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r_vec=x;

vector = ( sinint(0.5 .* (pi - r_vec .* Sq)) - sinint(0.5 .* (pi + r_vec .* Sq)) ) ./ (r_vec) ;

matrix = diag(vector);

Hp1s_si1 = matrix;

end

%Trans p1 -> s : (si - si)/r^2
function Hp1s_si2=Hp1toQstrans_si2(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r_vec=x;

vector = ( sinint(0.5 .* (pi - r_vec .* Sq)) - sinint(0.5 .* (pi + r_vec .* Sq)) ) ./ (r_vec.^2) ;

matrix = diag(vector);

Hp1s_si2 = matrix;

end

%Trans p1 -> s : (si - si)/r^3
function Hp1s_si3=Hp1toQstrans_si3(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r_vec=x;

vector = ( sinint(0.5 .* (pi - r_vec .* Sq)) - sinint(0.5 .* (pi + r_vec .* Sq)) ) ./ (r_vec.^3) ;

matrix = diag(vector);

Hp1s_si3 = matrix;

end

%Trans p1 -> s : (si - si)/r^4
function Hp1s_si4=Hp1toQstrans_si4(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r_vec=x;

vector = ( sinint(0.5 .* (pi - r_vec .* Sq)) - sinint(0.5 .* (pi + r_vec .* Sq)) ) ./ (r_vec.^4) ;

matrix = diag(vector);

Hp1s_si4 = matrix;

end
