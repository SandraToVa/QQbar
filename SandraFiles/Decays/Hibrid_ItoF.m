% FITXER ACTUALITZAT les funcions estan calculades per als dos operadors
% hibrids \partial_0\partial_i i \partial_z\partial_i


%These functions are the angular dependent part (theta and fi) integrated
%with the corresponendt spherical harmonics of the states we want to
%compute the decay with of. They apear in the mathematica file
%"PhD_I_i-f_hibrid.nb"
%We have still to do the numeric integration which we do in the other
%programs

%Create a function to call all the functions in here:

function res = Hibrid_ItoF(x)
    disp(['Input value: ', x]);
    %%%%%%%%%% HQ %%%%%%%%%%%%
    % Transition p1 -> s
    if strcmp(x, 'HQp1m1tos0I01')
        res = @Hp1M1toQs0trans_I01;
    elseif strcmp(x, 'HQp1m1tos0IN0')
        res = @Hp1M1toQs0trans_IN0;
    elseif strcmp(x, 'HQp1m0tos0I11')
        res = @Hp1M0toQs0trans_I11;

    % Transition p1 -> d
    %m final =2
    elseif strcmp(x, 'HQp1m1tod2I10')
        res = @Hp1M1toQd2trans_I10;
    %m final =1
    elseif strcmp(x, 'HQp1m1tod1I11')
        res = @Hp1M1toQd1trans_I11;
    elseif strcmp(x, 'HQp1m0tod1I0N')
        res = @Hp1M0toQd1trans_I0N;
    %m final =0
    elseif strcmp(x, 'HQp1m1tod0I01')
        res = @Hp1M1toQd0trans_I01;
    elseif strcmp(x, 'HQp1m1tod0IN0')
        res = @Hp1M1toQd0trans_IN0;
    elseif strcmp(x, 'HQp1m0tod0I11')
        res = @Hp1M0toQd0trans_I11;

    % Transition (s/d)1 -> p
    %m final =1
    elseif strcmp(x, 'HQsd1m1tp1I11')
        res = @Hsd1M1toQp1trans_I11;
    elseif strcmp(x, 'HQsd1m1tp1INN')
        res = @Hsd1M1toQp1trans_INN;
    elseif strcmp(x, 'HQsd1m0tp1I10')
        res = @Hsd1M0toQp1trans_I10;
    elseif strcmp(x, 'HQsd1m0tp1I0N')
        res = @Hsd1M0toQp1trans_I0N;
    %m final =0
    elseif strcmp(x, 'HQsd1m1tp0I01')
        res = @Hsd1M1toQp0trans_I01;
    elseif strcmp(x, 'HQsd1m1tp0IN0')
        res = @Hsd1M1toQp0trans_IN0;
    elseif strcmp(x, 'HQsd1m0tp0I11')
        res = @Hsd1M0toQp0trans_I11;
    end
    disp(['Output function handle: ', func2str(res)]);
end

%Recordo que en totes les expressions s'inclou un factor 1/r que fa que els
%exponencials de r sigin diferents als del document de mathematica
%"PhD_I_i-t_hibrid.nb"

%%%%%%%%% p0  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For p0->p transitions

%No hi ha transicions des de p0

%%%%%%%%% p1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For p1->s transitions

%Trans p1 (J=1,L=1,M=1) -> s (L=0,M=0)per a I01
function Hp1M1s_I01=Hp1M1toQs0trans_I01(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r_vec=x;
%Aproximation
%r=diag(x);
%Hp1M1s_I01=1/(sqrt(6)*pi^2) - sqrt(3/2)*(-8+pi^2).*r.^2.*Sq^2./(40*pi^4);
%Hole function
A_vector = 1 ./ ((r_vec.^3) .* (Sq.^3));
B = - sqrt(3/2) / 2;
C_vector = 4 .* sin(0.5 .* r_vec .* Sq) + ...
               pi .* sinint(0.5 .* (pi - r_vec .* Sq)) - ...
               pi .* sinint(0.5 .* (pi + r_vec .* Sq));

A = diag(A_vector);
C = diag(C_vector);
Hp1M1s_I01= A * B * C;
end

%Trans p1 (J=1,L=1,M=1) -> s (L=0,M=0) per a IN0
function Hp1M1s_IN0=Hp1M1toQs0trans_IN0(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r_vec=x;
%Aproximation
%r=diag(x);
%Hp1M1s_IN0=1/(sqrt(3)*pi^2) + (-8+pi^2).*r.^2.*Sq^2./(40*sqrt(3)*pi^4);
%Hole function
A_vector = 1 ./ ((r_vec.^3) .* (Sq.^3));
B = sqrt(3) / (4 * pi);
C_vector = 4 * pi .* sin(0.5 .* r_vec .* Sq) + ...
               (pi - r_vec .* Sq) .* (pi + r_vec .* Sq) .* ...
               (sinint(0.5 .* (pi - r_vec .* Sq)) - sinint(0.5 .* (pi + r_vec .* Sq)));

A = diag(A_vector);
C = diag(C_vector);
Hp1M1s_IN0= A * B * C;
end

% Trans p1 (J=1,L=1,M=0) -> s0 (L=0,M=0) for I11
function Hp1M0s0_I11 = Hp1M0toQs0trans_I11(x, Ei, Ef, M)
% x is the system = r in the string
% E = Ef - Ei of the transition
% M = dipion invariant mass
E = Ef - Ei;
Sq = sqrt(E^2 - M^2);
r_vec = x;

% Approximation:
% Hp1M0s0_I11 = 1/(sqrt(3)*pi^2) + (-8 + pi^2)*r.^2.*Sq^2 / (40*sqrt(3)*pi^4);

% Full function
A_vector = 1 ./ ((r_vec.^3) .* (Sq.^3));
B = sqrt(3) / (4 * pi);
C_vector = 4 * pi .* sin(0.5 .* r_vec .* Sq) + ...
    (pi - r_vec .* Sq) .* (pi + r_vec .* Sq) .* ...
    (sinint(0.5 .* (pi - r_vec .* Sq)) - sinint(0.5 .* (pi + r_vec .* Sq)));

A = diag(A_vector);
C = diag(C_vector);
Hp1M0s0_I11 = A * B * C;
end

% Trans p1 (J=1,L=1,M=0) -> s0 (L=0,M=0) for INN
%SAME TRANSITION AS P1 M=1 TO S0 I11 SAME SIGN


% Trans p1 (J=1,L=1,M=N) -> s0 (L=0,M=0) for I0N
%SAME TRANSITION AS P1 M=1 TO S0 I01 DIFERENT SIGN


% Trans p1 (J=1,L=1,M=N) -> s0 (L=0,M=0) for I10
%SAME TRANSITION AS P1 M=1 TO S0 IN0 SAME SIGN



% For p1->d transitions

% Trans p1 (J=1,L=1,M=1) -> d2 (L=2,M=2) for I10
function Hp1M1d2_I10 = Hp1M1toQd2trans_I10(x, Ei, Ef, M)
E = Ef - Ei;
Sq = sqrt(E^2 - M^2);
r_vec = x;

% Approximation:
% Hp1M1d2_I10 = 1/(sqrt(5/2)*pi^2) - sqrt(10)*(pi^2 - 8)*r.^2.*Sq.^2 / (28*pi^4);

A_vector = 1 ./ ((r_vec.^5) .* (Sq.^5));
B = - (3 * sqrt(5/2)) / (8 * pi);
C_vector = ...
    16 * pi * r_vec .* Sq .* cos(0.5 * r_vec .* Sq) + ...
    4 * pi * (-8 + pi^2 - r_vec.^2 .* Sq.^2) .* sin(0.5 * r_vec .* Sq) + ...
    (pi^2 - r_vec.^2 .* Sq.^2).^2 .* ...
    (sinint(0.5 * (pi - r_vec .* Sq)) - sinint(0.5 * (pi + r_vec .* Sq)));

A = diag(A_vector);
C = diag(C_vector);
Hp1M1d2_I10 = A * B * C;
end

% Trans p1 (J=1,L=1,M=1) -> d1 (L=2,M=1) for I11
function Hp1M1d1_I11 = Hp1M1toQd1trans_I11(x, Ei, Ef, M)
E = Ef - Ei;
Sq = sqrt(E^2 - M^2);
r_vec = x;

% Approximation:
% Hp1M1d1_I11 = -1/(sqrt(5)*pi^2) + 3*(pi^2 - 8)*r.^2.*Sq.^2 / (56*sqrt(5)*pi^4);

A_vector = 1 ./ ((r_vec.^5) .* (Sq.^5));
B = - (3 * sqrt(5)) / 4;
C_vector = ...
    16 * r_vec .* Sq .* cos(0.5 * r_vec .* Sq) + ...
    4 * (-8 + pi^2) .* sin(0.5 * r_vec .* Sq) + ...
    pi * (pi - r_vec .* Sq) .* (pi + r_vec .* Sq) .* ...
    (sinint(0.5 * (pi - r_vec .* Sq)) - sinint(0.5 * (pi + r_vec .* Sq)));

A = diag(A_vector);
C = diag(C_vector);
Hp1M1d1_I11 = A * B * C;
end

% Trans p1 (J=1,L=1,M=0) -> d1 (L=2,M=1) for I0N
function Hp1M0d1_I0N = Hp1M0toQd1trans_I0N(x, Ei, Ef, M)
E = Ef - Ei;
Sq = sqrt(E^2 - M^2);
r_vec = x;

% Approximation:
% Hp1M0d1_I0N = -1/(sqrt(10)*pi^2) + 3*(pi^2 - 8)*r.^2.*Sq.^2 / (56*sqrt(10)*pi^4);

A_vector = 1 ./ ((r_vec.^5) .* (Sq.^5));
B = - (3 * sqrt(5/2)) / 4;
C_vector = ...
    16 * r_vec .* Sq .* cos(0.5 * r_vec .* Sq) + ...
    4 * (-8 + pi^2) .* sin(0.5 * r_vec .* Sq) + ...
    pi * (pi - r_vec .* Sq) .* (pi + r_vec .* Sq) .* ...
    (sinint(0.5 * (pi - r_vec .* Sq)) - sinint(0.5 * (pi + r_vec .* Sq)));

A = diag(A_vector);
C = diag(C_vector);
Hp1M0d1_I0N = A * B * C;
end

% Trans p1 (J=1,L=1,M=1) -> d0 (L=2,M=1) for I01
function Hp1M1d0_I01 = Hp1M1toQd0trans_I01(x, Ei, Ef, M)
E = Ef - Ei;
Sq = sqrt(E^2 - M^2);
r_vec = x;

% Approximation:
% Hp1M1d0_I01 = 1/(sqrt(15/2)*pi^2) - sqrt(3/10)*(pi^2 - 8)*r.^2.*Sq.^2 / (14*pi^4);

A_vector = 1 ./ ((r_vec.^5) .* (Sq.^5));
B = sqrt(15/2) / 4;
C_vector = ...
    -48 * r_vec .* Sq .* cos(0.5 * r_vec .* Sq) - ...
    4 * (-24 + 3 * pi^2 + 2 * r_vec.^2 .* Sq.^2) .* sin(0.5 * r_vec .* Sq) - ...
    pi * (3 * pi^2 - r_vec.^2 .* Sq.^2) .* ...
    (sinint(0.5 * (pi - r_vec .* Sq)) - sinint(0.5 * (pi + r_vec .* Sq)));

A = diag(A_vector);
C = diag(C_vector);
Hp1M1d0_I01 = A * B * C;
end

% Trans p1 (J=1,L=1,M=1) -> d0 (L=2,M=1) for IN0
function Hp1M1d0_IN0 = Hp1M1toQd0trans_IN0(x, Ei, Ef, M)
E = Ef - Ei;
Sq = sqrt(E^2 - M^2);
r_vec = x;

% Approximation:
% Hp1M1d0_IN0 = -1/(sqrt(15)*pi^2) - (pi^2 - 8)*r.^2.*Sq.^2 / (56*sqrt(15)*pi^4);

A_vector = 1 ./ ((r_vec.^5) .* (Sq.^5));
B = sqrt(15) / (8 * pi);
C_vector = ...
    48 * pi * r_vec .* Sq .* cos(0.5 * r_vec .* Sq) + ...
    4 * pi * (-24 + 3 * pi^2 - r_vec.^2 .* Sq.^2) .* sin(0.5 * r_vec .* Sq) + ...
    (3 * pi^4 - 4 * pi^2 * r_vec.^2 .* Sq.^2 + r_vec.^4 .* Sq.^4) .* ...
    (sinint(0.5 * (pi - r_vec .* Sq)) - sinint(0.5 * (pi + r_vec .* Sq)));

A = diag(A_vector);
C = diag(C_vector);
Hp1M1d0_IN0 = A * B * C;
end

% Trans p1 (J=1,L=1,M=0) -> d0 (L=2,M=0) for I11
function Hp1M0d0_I11 = Hp1M0toQd0trans_I11(x, Ei, Ef, M)
E = Ef - Ei;
Sq = sqrt(E^2 - M^2);
r_vec = x;

% Approximation:
% Hp1M0d0_I11 = -1/(sqrt(15)*pi^2) - (pi^2 - 8)*r.^2.*Sq.^2 / (56*sqrt(15)*pi^4);

A_vector = 1 ./ ((r_vec.^5) .* (Sq.^5));
B = sqrt(15) / (8 * pi);
C_vector = ...
    48 * pi * r_vec .* Sq .* cos(0.5 * r_vec .* Sq) + ...
    4 * pi * (-24 + 3 * pi^2 - r_vec.^2 .* Sq.^2) .* sin(0.5 * r_vec .* Sq) + ...
    (3 * pi^4 - 4 * pi^2 * r_vec.^2 .* Sq.^2 + r_vec.^4 .* Sq.^4) .* ...
    (sinint(0.5 * (pi - r_vec .* Sq)) - sinint(0.5 * (pi + r_vec .* Sq)));

A = diag(A_vector);
C = diag(C_vector);
Hp1M0d0_I11 = A * B * C;
end

% Trans p1 (J=1,L=1,M=0) -> d0 (L=2,M=0) for INN
%SAME TRANSITION AS P1 M=0 TO D0 I11 SAME SIGN

% Trans p1 (J=1,L=1,M=-1) -> d0 (L=2,M=0) for I0N
%SAME TRANSITION AS P1 M=1 TO D0 I01 DIFERENT SIGN

% Trans p1 (J=1,L=1,M=-1) -> d0 (L=2,M=0) for I10
%SAME TRANSITION AS P1 M=1 TO D0 IN0 SAME SIGN

% Trans p1 (J=1,L=1,M=0) -> d-1 (L=2,M=-1) for I01
%SAME TRANSITION AS P1 M=0 TO D1 I0N DIFERENT SIGN

% Trans p1 (J=1,L=1,M=-1) -> d-1 (L=2,M=-1) for INN
%SAME TRANSITION AS P1 M=1 TO D1 I11 SAME SIGN

% Trans p1 (J=1,L=1,M=N) -> d-2 (L=2,M=-2) for IN0
%SAME TRANSITION AS P1 M=1 TO D2 I10 SAME SIGN


%%%%%%%%% (s/d)1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For (s/d)1->p transitions

%Trans (s/d)1 (J=1,L=0,2,M=1) -> p1 (L = 1, M = 1)
%Sandwitch: 1st row L=0 (s) and 2nd row L=2 (d) part of the hyrbid
function Hsd1M1p1_I11 = Hsd1M1toQp1trans_I11(x,Ei,Ef,M)
E = Ef - Ei;
Sq = sqrt(E^2 - M^2);
r_vec = x;
len=length(x);
Hsd1M1p1_I11 = zeros(len,len,2);

%L=0 (part s del hybrid)

% Approximation:
% Hsd1M1p1_I11(:,:,1) = sqrt(2/3)/(pi^2) - (pi^2 - 8)*r.^2.*Sq.^2 / (20*sqrt(6)*pi^4);

A_vector = 1 ./ ((r_vec.^3) .* (Sq.^3));
B = sqrt(3/2) / (2 * pi);
C_vector = ...
    4 * pi .* sin(0.5 * r_vec .* Sq) + ...
    (pi - r_vec .* Sq) .* (pi + r_vec .* Sq) .* ...
    (sinint(0.5 * (pi - r_vec .* Sq)) - sinint(0.5 * (pi + r_vec .* Sq)));

A = diag(A_vector);
C = diag(C_vector);
Hsd1M1p1_I11(:,:,1) = A * B * C;

%L=2 (part d del hybrid)

% Approximation:
% Hsd1M1p1_I11(:,:,2) = -1/(5*sqrt(3)*pi^2) - (pi^2 - 8)*r.^2.*Sq.^2 / (280*sqrt(3)*pi^4);

A_vector = 1 ./ ((r_vec.^5) .* (Sq.^5));
B = sqrt(3) / (8 * pi);
C_vector = ...
    48 * pi * r_vec .* Sq .* cos(0.5 * r_vec .* Sq) + ...
    4 * pi * (-24 + 3 * pi^2 - r_vec.^2 .* Sq.^2) .* sin(0.5 * r_vec .* Sq) + ...
    (3 * pi^4 - 4 * pi^2 * r_vec.^2 .* Sq.^2 + r_vec.^4 .* Sq.^4) .* ...
    (sinint(0.5 * (pi - r_vec .* Sq)) - sinint(0.5 * (pi + r_vec .* Sq)));

A = diag(A_vector);
C = diag(C_vector);
Hsd1M1p1_I11(:,:,2) = A * B * C;

end

%Trans (s/d)1 (J=1,L=0,2,M=1) -> p1 (L = 1, M = 1)
%Sandwitch: 1st row L=0 (s) and 2nd row L=2 (d) part of the hyrbid
function Hsd1M1p1_INN = Hsd1M1toQp1trans_INN(x,Ei,Ef,M)
E = Ef - Ei;
Sq = sqrt(E^2 - M^2);
r_vec = x;
len=length(x);
Hsd1M1p1_INN = zeros(len,len,2);

%L=0 (part s del hybrid)

% Approximation:
% Hsd1M1p1_INN(:,:,1) = 0;

%L=2 (part d del hybrid)

% Approximation:
% Hsd1M1p1_I11(:,:,2) = 2*sqrt(3)/(5*pi^2) - sqrt(3)*(pi^2 - 8)*r.^2.*Sq.^2 / (140*pi^4);

A_vector = 1 ./ ((r_vec.^5) .* (Sq.^5));
B = -3* sqrt(3) / (8 * pi);
C_vector = ...
    16 * pi * r_vec .* Sq .* cos(0.5 * r_vec .* Sq) + ...
    4 * pi * (-8 + pi^2 - r_vec.^2 .* Sq.^2) .* sin(0.5 * r_vec .* Sq) + ...
    (pi^2 - r_vec.^2 .* Sq.^2).^2 .* ...
    (sinint(0.5 * (pi - r_vec .* Sq)) - sinint(0.5 * (pi + r_vec .* Sq)));

A = diag(A_vector);
C = diag(C_vector);
Hsd1M1p1_INN(:,:,2) = A * B * C;

end


% Transition: (s/d)1 (J=1,L=0,2,M=0) -> p1 (L=1, M=1)
% Sandwich: 1st row L=0 (s) and 2nd row L=2 (d) part of the hybrid
function Hsd1M0p1_I10 = Hsd1M0toQp1trans_I10(x,Ei,Ef,M)
E = Ef - Ei;
Sq = sqrt(E^2 - M^2);
r_vec = x;
len = length(x);
Hsd1M0p1_I10 = zeros(len,len,2);

% ----- L=0 (s-wave part of hybrid) -----
% Approximation:
% Hsd1M0p1_I10(:,:,1) = sqrt(2/3)/(pi^2) - (pi^2 - 8)*r.^2.*Sq.^2 / (20*sqrt(6)*pi^4);

A_vector = 1 ./ ((r_vec.^3) .* (Sq.^3));
B = sqrt(3/2) / (2 * pi);
C_vector = ...
    4 * pi .* sin(0.5 * r_vec .* Sq) + ...
    (pi - r_vec .* Sq) .* (pi + r_vec .* Sq) .* ...
    (sinint(0.5 * (pi - r_vec .* Sq)) - sinint(0.5 * (pi + r_vec .* Sq)));

A = diag(A_vector);
C = diag(C_vector);
Hsd1M0p1_I10(:,:,1) = A * B * C;

% ----- L=2 (d-wave part of hybrid) -----
% Approximation:
% Hsd1M0p1_I10(:,:,2) = -sqrt(3/2)/(5*pi^2) + 3*(pi^2 - 8)*r.^2.*Sq.^2 / (280*sqrt(6)*pi^4);

A_vector = 1 ./ ((r_vec.^5) .* (Sq.^5));
B = sqrt(3) / (4 * pi);
C_vector = ...
    -48 * pi * r_vec .* Sq .* cos(0.5 * r_vec .* Sq) + ...
    4 * pi * (24 - 3 * pi^2 + r_vec.^2 .* Sq.^2) .* sin(0.5 * r_vec .* Sq) - ...
    (3 * pi^4 - 4 * pi^2 * r_vec.^2 .* Sq.^2 + r_vec.^4 .* Sq.^4) .* ...
    (sinint(0.5 * (pi - r_vec .* Sq)) - sinint(0.5 * (pi + r_vec .* Sq)));

A = diag(A_vector);
C = diag(C_vector);
Hsd1M0p1_I10(:,:,2) = A * B * C;

end

% Transition: (s/d)1 (J=1,L=0,2,M=0) -> p1 (L=1, M=1)
% Sandwich: 1st row L=0 (s) and 2nd row L=2 (d) part of the hybrid
function Hsd1M0p1_I0N = Hsd1M0toQp1trans_I0N(x,Ei,Ef,M)
E = Ef - Ei;
Sq = sqrt(E^2 - M^2);
r_vec = x;
len = length(x);
Hsd1M0p1_I0N = zeros(len,len,2);

% ----- L=0 (s-wave part of hybrid) -----
% Approximation:
% Hsd1M0p1_I0N(:,:,1) = 0;

% ----- L=2 (d-wave part of hybrid) -----
% Approximation:
% Hsd1M0p1_I0N(:,:,2) = 2/(5*sqrt(3)*pi^2) + (pi^2 - 8)*r.^2.*Sq.^2 / (140*sqrt(3)*pi^4);

A_vector = 1 ./ ((r_vec.^5) .* (Sq.^5));
B = -3 * sqrt(3/2) / (4);
C_vector = ...
    16 * r_vec .* Sq .* cos(0.5 * r_vec .* Sq) + ...
    4 * (-8 + pi^2) .* sin(0.5 * r_vec .* Sq) + ...
    pi * (pi - r_vec .* Sq) .* (pi + r_vec .* Sq) .* ...
    (sinint(0.5 * (pi - r_vec .* Sq)) - sinint(0.5 * (pi + r_vec .* Sq)));

A = diag(A_vector);
C = diag(C_vector);
Hsd1M0p1_I0N(:,:,2) = A * B * C;

end


% Transition: (s/d)1 (J=1,L=0,2,M=1) -> p0 (L=1, M=0)
% Sandwich: 1st row L=0 (s) and 2nd row L=2 (d) part of the hybrid
function Hsd1M1p0_I01 = Hsd1M1toQp0trans_I01(x,Ei,Ef,M)
E = Ef - Ei;
Sq = sqrt(E^2 - M^2);
r_vec = x;
len = length(x);
Hsd1M1p0_I01 = zeros(len,len,2);

% ----- L=0 (s-wave part of hybrid) -----
% Approximation:
% Hsd1M1p0_I01(:,:,1) = -1/sqrt(3)/pi^2 + (pi^2 - 8)*r.^2.*Sq.^2 / (40*sqrt(3)*pi^4);

A_vector = 1 ./ ((r_vec.^3) .* (Sq.^3));
B = sqrt(3) / 2;
C_vector = ...
    4 * sin(0.5 * r_vec .* Sq) + ...
    pi * (sinint(0.5 * (pi - r_vec .* Sq)) - sinint(0.5 * (pi + r_vec .* Sq)));

A = diag(A_vector);
C = diag(C_vector);
Hsd1M1p0_I01(:,:,1) = A * B * C;

% ----- L=2 (d-wave part of hybrid) -----
% Approximation:
% Hsd1M1p0_I01(:,:,2) = -sqrt(2/3)/5/pi^2 + (pi^2 - 8)*r.^2.*Sq.^2 / (70*sqrt(6)*pi^4);

A_vector = 1 ./ ((r_vec.^5) .* (Sq.^5));
B = sqrt(3/2) / 4;
C_vector = ...
    48 * r_vec .* Sq .* cos(0.5 * r_vec .* Sq) + ...
    4 * (-24 + 3 * pi^2 + 2 * r_vec.^2 .* Sq.^2) .* sin(0.5 * r_vec .* Sq) + ...
    pi * (3 * pi^2 - r_vec.^2 .* Sq.^2) .* ...
    (sinint(0.5 * (pi - r_vec .* Sq)) - sinint(0.5 * (pi + r_vec .* Sq)));

A = diag(A_vector);
C = diag(C_vector);
Hsd1M1p0_I01(:,:,2) = A * B * C;

end

% Transition: (s/d)1 (J=1,L=0,2,M=1) -> p0 (L=1, M=0)
% Sandwich: 1st row L=0 (s) and 2nd row L=2 (d) part of the hybrid
function Hsd1M1p0_IN0 = Hsd1M1toQp0trans_IN0(x,Ei,Ef,M)
E = Ef - Ei;
Sq = sqrt(E^2 - M^2);
r_vec = x;
len = length(x);
Hsd1M1p0_IN0 = zeros(len,len,2);

% ----- L=0 (s-wave part of hybrid) -----
% Approximation:
% Hsd1M1p0_IN0(:,:,1) = 0;

% ----- L=2 (d-wave part of hybrid) -----
% Approximation:
% Hsd1M1p0_IN0(:,:,2) = -sqrt(3)/5/pi^2 + 3*(pi^2 - 8)*r.^2.*Sq.^2 / (280*sqrt(3)*pi^4);

A_vector = 1 ./ ((r_vec.^5) .* (Sq.^5));
B = -3 * sqrt(3) / 4;
C_vector = ...
    16 * r_vec .* Sq .* cos(0.5 * r_vec .* Sq) + ...
    4 * (-8 + pi^2) .* sin(0.5 * r_vec .* Sq) + ...
    pi * (pi - r_vec .* Sq) .* (pi + r_vec .* Sq) .* ...
    (sinint(0.5 * (pi - r_vec .* Sq)) - sinint(0.5 * (pi + r_vec .* Sq)));

A = diag(A_vector);
C = diag(C_vector);
Hsd1M1p0_IN0(:,:,2) = A * B * C;

end


% Transition: (s/d)1 (J=1,L=0,2,M=0) -> p0 (L=1, M=0)
% Sandwich: 1st row L=0 (s) and 2nd row L=2 (d) part of the hybrid
function Hsd1M0p0_I11 = Hsd1M0toQp0trans_I11(x,Ei,Ef,M)
E = Ef - Ei;
Sq = sqrt(E^2 - M^2);
r_vec = x;
len = length(x);
Hsd1M0p0_I11 = zeros(len,len,2);

% ----- L=0 (s-wave part of hybrid) -----
% Approximation: 0
% Hsd1M0p0_I11(:,:,1) = 0;

% ----- L=2 (d-wave part of hybrid) -----
% Approximation:
% Hsd1M0p0_I11(:,:,2) = -sqrt(3)/5/pi^2 + 3*(pi^2 - 8)*r.^2.*Sq.^2 / (280*sqrt(3)*pi^4);

A_vector = 1 ./ ((r_vec.^5) .* (Sq.^5));
B = -3 * sqrt(3) / 4;
C_vector = ...
    16 * r_vec .* Sq .* cos(0.5 * r_vec .* Sq) + ...
    4 * (-8 + pi^2) .* sin(0.5 * r_vec .* Sq) + ...
    pi * (pi - r_vec .* Sq) .* (pi + r_vec .* Sq) .* ...
    (sinint(0.5 * (pi - r_vec .* Sq)) - sinint(0.5 * (pi + r_vec .* Sq)));

A = diag(A_vector);
C = diag(C_vector);
Hsd1M0p0_I11(:,:,2) = A * B * C;

end

% Transition: (s/d)1 (J=1,L=0,2,M=0) -> p0 (L=1, M=0) for INN
%SAME TRANSITION AS (s/d)1 (J=1,L=0,2,M=0) -> p0 (L=1, M=0) for I11 (=SIGN)

% Transition: (s/d)1 (J=1,L=0,2,M=-1) -> p0 (L=1, M=0) for I0N
%SAME TRANSITION AS (s/d)1 (J=1,L=0,2,M=1) -> p0 (L=1, M=0) for I01 (=SIGN)

% Transition: (s/d)1 (J=1,L=0,2,M=-1) -> p0 (L=1, M=0) for I10
%SAME TRANSITION AS (s/d)1 (J=1,L=0,2,M=1) -> p0 (L=1, M=0) for IN0 (-SIGN)

% Transition: (s/d)1 (J=1,L=0,2,M=0) -> p0 (L=1, M=-1) for IN0
%SAME TRANSITION AS (s/d)1 (J=1,L=0,2,M=0) -> p0 (L=1, M=1) for I10 (-SIGN)

% Transition: (s/d)1 (J=1,L=0,2,M=0) -> p0 (L=1, M=-1) for I01
%SAME TRANSITION AS (s/d)1 (J=1,L=0,2,M=0) -> p0 (L=1, M=1) for I0N (=SIGN)

% Transition: (s/d)1 (J=1,L=0,2,M=-1) -> p0 (L=1, M=-1) for INN
%SAME TRANSITION AS (s/d)1 (J=1,L=0,2,M=1) -> p0 (L=1, M=1) for I11 (-SIGN)

% Transition: (s/d)1 (J=1,L=0,2,M=-1) -> p0 (L=1, M=-1) for I11
%SAME TRANSITION AS (s/d)1 (J=1,L=0,2,M=1) -> p0 (L=1, M=1) for INN (-SIGN)



