%Code to compute the expected value of r, 1/r, 1/r^2 etc


% One of its uses is to check if for a transition from híbrid to quarkonium the value of
 %<Psi_QQ| 1/r |Psi_HQ> with r = x in the code is >> than DeltaE of the QQ
 %and HQ states.

 %We use hibrid states but not splited to the hiperfine so we only need r0,
 %mq and spin

setr0(3.964)
load("dades.mat","m_c","m_b")
setm_q(m_c)

%I will try 2 transitions 1(s/d)->1p and 1p1->2s for bottomonium
% Initialize cell array to store variable-sized matrices
wf_hq = cell(1, 2);         %for each cell there will be a matrix(wf made by diferent states at diferent x positions)
wf_qq = cell(1, 2);  
xvector = cell(1, 2);

E = zeros(2,2);
%E = zeros(2,3);
%Rvalues = zeros(1,2);
%R2values = zeros(1,2);
%DEvalues = zeros(1,2);
%DERvalues = zeros(1,2);


%We compute the hibrid functions
setspin(1)

%Hibrid states without the hiperfine:
[aux,auxwf,x]=QuarkoniumS0J1(m_q,spin); %1(s/d)1 i 1p1
% 1--
E(1,1)=aux(1);
E(1,2)=aux(2);
wf_hq{1} = auxwf(1:2, :, 1); %Primer element 1:2 si son (s/d) o (p/f)
wf_hq{2} = auxwf(3, :, 2); %Primer element 1:2 si son (s/d) o (p/f)
xvector{1}=x;


%We compute the quarkonium functions
setspin(0)

%Hibrid states without the hiperfine:
[aux,auxwf,x]=QuarkoniumS0J1(m_q,spin);
% (0,1,2)++ p
E(2,1)=aux(1);
wf_qq{1} = auxwf(1, :, 1);
xvector{2}=x;

[aux,auxwf,x]=QuarkoniumS0J0(m_q,spin);
% 1-- s
E(2,2)=aux(2);
wf_qq{2} = auxwf(1, :, 2);
xvector{2}=x;




%After computing the quarkonium with same x we proceed

[sizeXhq,lengthXhq] = size(xvector{1});
[sizeXqq,lengthXqq] = size(xvector{2});

if sizeXhq ~= sizeXqq
    disp('Number of rows in x not equal')
    disp(sizeXhq);
    disp(sizeXqq);
end
if lengthXhq ~= lengthXqq
    disp('Number of columns in x not equal')
    disp(lengthXhq);
    disp(lengthXqq);
end


dim=2;
%The matrix is the same for both transitions
Rmatrix = 1 ./ x ;
R2matrix = 1 ./ x.^2 ;


%Values of the diference of energy
%DEvalues(1) = E(2,1)-E(1,1); %1p - 1(s/d)1
DEvalues(2) = E(2,2)-E(1,2); %2s - 1p1


%The value of the 1/r sandwitch
%Rvalues(1) = computeSandwitch(Rmatrix,wf_hq{1},wf_qq{1},x);
Rvalues(2) = computeSandwitch(Rmatrix,wf_hq{2},wf_qq{2},x);
%The value of the 1/r^2 sandwitch
%R2values(1) = computeSandwitch(R2matrix,wf_hq{1},wf_qq{1},x);
R2values(2) = computeSandwitch(R2matrix,wf_hq{2},wf_qq{2},x);
%Values of the diference of energy/r
%DERvalues(1) = DEvalues(1)*Rvalues(1);
DERvalues(2) = DEvalues(2)*Rvalues(2);



%The result should be the same as using with Op the Rmatirx 3 times
%Op=zeros(lengthXhq,lengthXhq,3);
%Op(:,:,1)=Rmatrix;
%Op(:,:,2)=Rmatrix;
%R2sdto1p=ExpectedValue(3,1,Op,@QuarkoniumS0J1,@QuarkoniumS0J1,m_q,1,0);
%R2sdto2p=ExpectedValue(3,2,Op,@QuarkoniumS0J1,@QuarkoniumS0J1,m_q,1,0);



 function Rvalue=computeSandwitch(O,wf_hq,wf_qq,x)
%Input:     Rmatrix = matrix that will be sanwihed between the
%                       wavefunctions
%           wf_hq = 2 row wavefunction 1st row is s and 2n is d for hibrids
%           wf_qq = 3 row wavefunction only 1 row is useful for quarkonium
%
%
%Output:    Rvalue = result

    [elements,~] = size(wf_hq);

    Rvector=zeros(1,elements);
    Rvalue=0;

    for N=1:elements

        Yf=wf_qq; %row vector fianl
        Yi=wf_hq(N,:); %row vector inicial


        result_matrix =  trapz(x, conj(Yf) .* (O .* Yi));

        normf = sqrt(trapz(x,abs(Yf).^2));
        normi = sqrt(trapz(x,abs(Yi).^2));

        Rvector(N) = result_matrix / (normf * normi);


        Rvalue = Rvalue + Rvector(N);

    end 


 end




%-----------------------------------------------------------------------

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




% VALOR DEL Vhf
function x1=k1
global v1
x1=v1;
end 

function setk1(val1)
global v1
v1 = val1;
end

% VALOR DEL Vhf2
function x2=k2
global v2
x2=v2;
end 

function setk2(val2)
global v2
v2 = val2;
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

% VALOR DEL -g\Lambda' (el signe negatiu ja esta en la Vsp)
function x4=L1
global v4
x4=v4;
end 

function setL1(val4)
global v4
v4 = val4;
end

% VALOR DEL -g\Lambda''' (per a la Vpp)
function x5=L3
global v5
x5=v5;
end 

function setL3(val5)
global v5
v5 = val5;
end

%Usant la interpolació o només les llargues distàncies
%l=1; llargues distàncies
%l=0; interpolació
function x6=l
global v6
x6=v6;
end

function setl(val6)
global v6
v6=val6;
end