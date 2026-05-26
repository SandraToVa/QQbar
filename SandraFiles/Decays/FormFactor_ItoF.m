
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The update with all the possible transitions is missin!!!
%s->s and s->d and p->p are updated with transitions including F1 and F-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FITXER ACTUALITZAT després de l'error per falta d'un operador

%All the form factors integrated
%These functions are the angular dependent part (theta and fi) integrated
%with the corresponendt spherical harmonics of the states we want to
%compute the decay with of.
%We have still to do the numeric integration which we do in the other
%programs

%Create a function to call all the functions in here:
function res = FormFactor_ItoF(x)
    disp(['Input value: ', x]);
    %%%%%%%%%% QQ %%%%%%%%%%%%
    if strcmp(x, 'QQS0toS0_F/')
        res = @StoStransF;
    elseif strcmp(x, 'QQS0toS0_Fc')
        res = @StoStransFC;
        %
    elseif strcmp(x, 'QQP0toP0_F/')
        res = @P0toP0transF;
    elseif strcmp(x, 'QQP0toP0_Fc')
        res = @P0toP0transFC;
    elseif strcmp(x, 'QQP1toP1_F/')
        res = @P1toP1transF;
    elseif strcmp(x, 'QQP1toP1_Fc')
        res = @P1toP1transFC;
    elseif strcmp(x, 'QQPNtoP0_F1')
        res = @P1toP0transF1;
    elseif strcmp(x, 'QQP1toP0_FN')
        res = @P1toP0transFN;
    elseif strcmp(x, 'QQP0toP1_F1')
        res = @P0toP1transF1;
    elseif strcmp(x, 'QQP0toPN_FN')
        res = @P0toP1transFN;
    elseif strcmp(x, 'QQP1toP1_Fs')
        res = @P1toP1transFS;
    elseif strcmp(x, 'QQP1toP1_Fx')
        res = @P1toP1transFSX;
        %
    elseif strcmp(x, 'QQD0toS0_F/')
        res = @D0toStransF;
    elseif strcmp(x, 'QQD0toS0_Fc')
        res = @D0toStransFC;    
    elseif strcmp(x, 'QQD2toS0_Fs')
        res = @D2toStransFS;
    elseif strcmp(x, 'QQD2toS0_Fx')
        res = @D2toStransFSX;
    elseif strcmp(x, 'QQD2toS0_F1')
        res = @D2toStransF1;
    elseif strcmp(x, 'QQD2toS0_FN')
        res = @D2toStransFN;
        %
    elseif strcmp(x, 'QQD0toD0_F/')
        res = @D0toD0transF;
    elseif strcmp(x, 'QQD0toD0_Fc')
        res = @D0toD0transFC;    
    elseif strcmp(x, 'QQD1toD1_F/')
        res = @D1toD1transF;
    elseif strcmp(x, 'QQD1toD1_Fc')
        res = @D1toD1transFC;
    elseif strcmp(x, 'QQD2toD2_F/')
        res = @D2toD2transF;
    elseif strcmp(x, 'QQD2toD2_Fc')
        res = @D2toD2transFC;
    end
    disp(['Output function handle: ', func2str(res)]);
end


%First all the functions in Quarkonium -> Quarkonium +2Pions transitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%With the functions below we want to create a matrix "operator"
%This matrix is diagonal, however I sum and substract elements to each
%component of the matrix resuylting in a non diagonal matrix.
%In the function ComputeExpVal I only use the matrix elements of these
%matrix and hence it does not matter that the matrix from the functions are
%not diagonal
%But IMPORTANT to remember!!!

%This results are the ones form the PhD_itof_Quarkonium mathmetaica
%notebook
%To obtain the I_ItoF for each form factor we still need to integrate
%numerically with the wave functions


%For s->s transitions
%In here we have contributions from F and Fc
%Totes les constants estan fora dels form factors

%F (l=0,m=0 -> l=0,m=0)
function Fss=StoStransF(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%Fss=r/2-(r^3*Sq^2)/144;
%Hole function
Fss=sinint((Sq/2).*r)/Sq;
end

%Fc (l=0,m=0 -> l=0,m=0)
function FCss=StoStransFC(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%FCss=r/6 - (r^3 * Sq^2)/240;
%Hole function
A= diag(diag((1 ./ (r.^2 .* Sq.^3))));
B = -2 * r .* Sq .* cos((Sq / 2) .* r) + 4 .* sin((Sq / 2) .* r);
FCss=A*B;
end


%For p->p transitions
%In here we have contributions from F, Fc and Fs i Fsx
%Depenent de la M utilitzada

%F (l=1,m=0 -> l=1,m=0)
function Fpp00=P0toP0transF(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%Fpp00=r/2 - (r^3 * Sq^2)/80;
%Hole function
A= diag(diag((1 ./ (r.^2 .* Sq.^3))));
B = -6 * (r .* Sq .* cos((Sq / 2) .* r) - 2 .* sin((Sq / 2) .* r));
Fpp00=A*B;
end

%Fc (l=1,m=0 -> l=1,m=0)
function FCpp00=P0toP0transFC(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%FCpp00=3*r/10 - (r^3 * Sq^2)/112;
%Hole function
A= diag(diag((1 ./ (r.^4 .* Sq.^5))));
B = -6 * r .* Sq .* (-24 + r.^2 * Sq^2) .* cos((Sq / 2) .* r) + 36 .* (-8 + r.^2 * Sq^2) .* sin((Sq / 2) .* r);
FCpp00=A*B;
end

%F (l=1,m=+-1 -> l=1,m=+-1)
function Fpp11=P1toP1transF(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%Fpp11=r/2 - (r^3 * Sq^2)/240;
%Hole function
A= diag(diag((1 ./ (2 * r.^2 .* Sq.^3))));
B = 3 * (2 * r .* Sq .* cos((Sq / 2) .* r) - 4 .* sin((Sq / 2) .* r) + r.^2 .* Sq^2 .* sinint((Sq / 2) .* r));
Fpp11=A*B;
end

%Fc (l=1,m=+-1 -> l=1,m=+-1)
function FCpp11=P1toP1transFC(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%FCpp11=r/10 - (r^3 * Sq^2)/560;
%Hole function
A= diag(diag((1 ./ (r.^4 .* Sq.^5))));
B = -12 * (6 * r .* Sq .* cos((Sq / 2) .* r) + (-12 + r.^2 .* Sq^2) .* sin((Sq / 2) .* r));
FCpp11=A*B;
end

%F1 (l=1,m=-1 -> l=1,m=0)
function F1ppN0=P1toP0transF1(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%FCpp11=r/10 - (r^3 * Sq^2)/560;
%Hole function
A= diag(diag((1 ./ (r.^4 .* Sq.^5))));
B = -12 * sqrt(2) * (6 * r .* Sq .* cos((Sq / 2) .* r) + (-12 + r.^2 .* Sq^2) .* sin((Sq / 2) .* r));
F1ppN0=A*B;
end

%FN (l=1,m=+1 -> l=1,m=0)
function FNpp10=P1toP0transFN(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%FCpp11=r/10 - (r^3 * Sq^2)/560;
%Hole function
A= diag(diag((1 ./ (r.^4 .* Sq.^5))));
B = 12 * sqrt(2) * (6 * r .* Sq .* cos((Sq / 2) .* r) + (-12 + r.^2 .* Sq^2) .* sin((Sq / 2) .* r));
FNpp10=A*B;
end

%F1 (l=1,m=0 -> l=1,m=1)
function F1pp01=P0toP1transF1(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%FCpp11=r/10 - (r^3 * Sq^2)/560;
%Hole function
A= diag(diag((1 ./ (r.^4 .* Sq.^5))));
B = 12 * sqrt(2) * (6 * r .* Sq .* cos((Sq / 2) .* r) + (-12 + r.^2 .* Sq^2) .* sin((Sq / 2) .* r));
F1pp01=A*B;
end

%FN (l=1,m=+1 -> l=1,m=0)
function FNpp0N=P0toP1transFN(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%FCpp11=r/10 - (r^3 * Sq^2)/560;
%Hole function
A= diag(diag((1 ./ (r.^4 .* Sq.^5))));
B = - 12 * sqrt(2) * (6 * r .* Sq .* cos((Sq / 2) .* r) + (-12 + r.^2 .* Sq^2) .* sin((Sq / 2) .* r));
FNpp0N=A*B;
end

%Fs (l=1,m=-1 -> l=1,m=+1)
function FSpp11=P1toP1transFS(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%FSpp11= -2*r/5 + (r^3 * Sq^2)/420;
%Hole function
A= diag(diag((1 ./ (2 * r.^4 .* Sq.^5))));
B = -3 * (2 * r .* Sq .* (24 + r.^2 .* Sq^2) .* cos((Sq / 2) .* r) + 4 .* (-24 + r.^2 .* Sq^2) .* sin((Sq / 2) .* r) + r.^4 .* Sq^4 .* sinint((Sq / 2) .* r));
FSpp11=A*B;
end

%Fsx (l=1,m=+1 -> l=1,m=-1)
function FSXpp11=P1toP1transFSX(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%FSXpp11= -2*r/5 + (r^3 * Sq^2)/420;
%Hole function
A= diag(diag((1 ./ (2 * r.^4 .* Sq.^5))));
B = -3 * (2 * r .* Sq .* (24 + r.^2 .* Sq^2) .* cos((Sq / 2) .* r) + 4 .* (-24 + r.^2 .* Sq^2) .* sin((Sq / 2) .* r) + r.^4 .* Sq^4 .* sinint((Sq / 2) .* r));
FSXpp11=A*B;
end


%For d->s transitions
%In here we have contributions from F, Fc and Fs i Fsx
%Depenent de la M utilitzada


%F (l=2,m=0 -> l=0,m=0)
function Fds00=D0toStransF(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%Fds00=-(r^3*Sq^2)/(72*sqrt(5));
%Hole function
A=diag(diag((1 ./ (r.^2 .* Sq.^3))));
B = -(sqrt(5) / 2) * (6 * r .* Sq .* cos((Sq / 2) .* r) - 12 .* sin((Sq / 2) .* r) + r.^2 .* Sq^2 .* sinint((Sq / 2) .* r));
Fds00=A*B;
end

%Fc (l=2,m=0 -> l=0,m=0)
function FCds00=D0toStransFC(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%FCds00= r/(3*sqrt(5)) - (r^3*Sq^2)/(84*sqrt(5));
%Hole function
A=diag(diag((1 ./ (r.^4 .* Sq.^5))));
B = -2 * sqrt(5) * ( r .* Sq .* (-36 + r.^2 * Sq^2) .* cos((Sq / 2) .* r) ...
      - 8 * (-9 + r.^2 * Sq^2) .* sin((Sq / 2) .* r) );
FCds00=A*B;
end

%F1 (l=2,m=-1 -> l=0,m=0)
function F1dsN0=D2toStransF1(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%FSds20= sqrt(2/15)*r -(r^3*Sq^2)/(84*sqrt(30));
%Hole function
A=diag(diag((1 ./ (r.^4 .* Sq.^5))));
B = - 4* sqrt(30) *  (6 .* r .* Sq .* cos((Sq / 2) .* r) + (-12 + r.^2 .* Sq^2) .* sin((Sq / 2) .* r) );
F1dsN0=A*B;
end

%FSX (l=2,m=1 -> l=0,m=0)
function FNds10=D2toStransFN(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%FSXds20= sqrt(2/15)*r -(r^3*Sq^2)/(84*sqrt(30));
%Hole function
A=diag(diag((1 ./ (r.^4 .* Sq.^5))));
B = 4* sqrt(30) *  (6 .* r .* Sq .* cos((Sq / 2) .* r) + (-12 + r.^2 .* Sq^2) .* sin((Sq / 2) .* r) );
FNds10=A*B;
end


%FS (l=2,m=-2 -> l=0,m=0)
function FSds20=D2toStransFS(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%FSds20= sqrt(2/15)*r -(r^3*Sq^2)/(84*sqrt(30));
%Hole function
A=diag(diag((1 ./ (2 .* r.^4 .* Sq.^5))));
B = sqrt(15 / 2) *  (2 .* r .* Sq .* (24 + r.^2 .* Sq^2) .* cos((Sq / 2) .* r) + 4 .* (-24 + r.^2 .* Sq^2) .* sin((Sq / 2) .* r) + r.^4 .* Sq^4 .* sinint((Sq / 2) .* r));
FSds20=A*B;
end

%FSX (l=2,m=2 -> l=0,m=0)
function FSXds20=D2toStransFSX(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%FSXds20= sqrt(2/15)*r -(r^3*Sq^2)/(84*sqrt(30));
%Hole function
A=diag(diag((1 ./ (2 * r.^4 .* Sq.^5))));
B = sqrt(15 / 2) * (2 * r .* Sq .* (24 + r.^2 .* Sq^2) .* cos((Sq / 2) .* r) + 4 .* (-24 + r.^2 .* Sq^2) .* sin((Sq / 2) .* r) + r.^4 .* Sq^4 .* sinint((Sq / 2) .* r));
FSXds20=A*B;
end


%For d->d transitions
%In here we have contributions from F and Fc

%F (l=2,m=0 -> l=2,m=0)
function Fdd00=D0toD0transF(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%Fdd00= r/2 - (11 * r^3 * Sq^2)/1008;
%Hole function
A=diag(diag((1 ./ (4 * r.^4 .* Sq.^5))));
B = 5 * ( -6 * r .* Sq .* (-72 + r.^2 .* Sq^2) .* cos((Sq / 2) .* r) + 12 .* (-72 + 7 * r.^2 .* Sq^2) .* sin((Sq / 2) .* r) + r.^4 .* Sq.^4 .* sinint((Sq / 2) .* r));
Fdd00=A*B;
end

%Fc (l=2,m=0 -> l=2,m=0)
function FCdd00=D0toD0transFC(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%FCdd00= (11/42)*r -(r^3*Sq^2)/112;
%Hole function
A=diag(diag((1 ./ (r.^6 .* Sq.^7))));
B = -10 * (r .* Sq .* (4320 - 144 * r.^2 .* Sq.^2 + r.^4 .* Sq.^4) .* cos((Sq / 2) .* r) - 2 .* (4320 - 504 * r.^2 .* Sq.^2 + 7 * r.^4 .* Sq.^4) .* sin((Sq / 2) .* r));
FCdd00=A*B;
end

%F (l=2,m=+-1 -> l=2,m=+-1)
function Fdd11=D1toD1transF(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%Fdd11= r/2 - (r^3 * Sq^2)/112;
%Hole function
A=diag(diag((1 ./ (r.^4 .* Sq.^5))));
B = -60 * (6 * r .* Sq .* cos((Sq / 2) .* r) + (-12 + r.^2 .* Sq.^2) .* sin((Sq / 2) .* r));
Fdd11=A*B;
end

%Fc (l=2,m=+-1 -> l=2,m=+-1)
function FCdd11=D1toD1transFC(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%FCdd11= (3/14)*r -(5*r^3*Sq^2)/1008;
%Hole function
A=diag(diag((1 ./ (r.^6 .* Sq.^7))));
B = -60 * (2 * r .* Sq .* (-240 + 7 * r.^2 .* Sq.^2) .* cos((Sq / 2) .* r) + (960 - 108 * r.^2 .* Sq.^2 + r.^4 .* Sq.^4) .* sin((Sq / 2) .* r));
FCdd11=A*B;
end

%F (l=2,m=+-2 -> l=2,m=+-2)
function Fdd22=D2toD2transF(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%Fdd22= r/2 - (r^3 * Sq^2)/336;
%Hole function
A=diag(diag((1 ./ (8 * r.^4 .* Sq.^5))));
B = 15 * (2 * r .* Sq .* (24 + r.^2 .* Sq.^2) .* cos((Sq / 2) .* r) + 4 .* (-24 + r.^2 .* Sq.^2) .* sin((Sq / 2) .* r) + r.^4 .* Sq.^4 .* sinint((Sq / 2) .* r));
Fdd22=A*B;
end

%Fc (l=2,m=+-2 -> l=2,m=+-2)
function FCdd22=D2toD2transFC(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%FCdd22= r/14 - (r^3*Sq^2)/1008;
%Hole function
A=diag(diag((1 ./ (r.^6 .* Sq.^7))));
B = 120 * (r .* Sq .* (-60 + r.^2 .* Sq.^2) .* cos((Sq / 2) .* r) - 12 .* (-10 + r.^2 .* Sq.^2) .* sin((Sq / 2) .* r));
FCdd22=A*B;
end