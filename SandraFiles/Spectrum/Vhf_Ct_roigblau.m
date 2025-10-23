%Valor r0
setr0(3.964)
%Constants ajustades
setk1(-0.0700)
setk2(0.0117)
%Valors Lambda g
setL1(-0.059)
setL3(-0.230)
% l= interpolació (0), llagures distncies (1), bad long distances (2)
setl(0)
%massa
load("dades.mat","m_c","m_b")
setm_q(m_c)
setspin(1) %s=1 for hybrids

%Totes les E en GeV
%Vector de valors de E calculades
a=zeros(1,14);

res=zeros(14,1);
%Vector de valors de E teòriques
t = [4.0296 3.8976 3.9286 4.0746 4.1106 4.1756 4.2386 4.4396 ...
     4.1116 4.1786 4.5136 4.1436 4.2306 4.2516];
%Vector error de energies teòriques
e = [0.0176 0.0186 0.0236 0.0216 0.0276 0.0186 0.0266 0.0466 ...
     0.0236 0.0276 0.0536 0.0256 0.0326 0.0346];


%Vector en los valors de la energia que necesito
[aux,~,~]=QuarkoniumS0J1(m_q,spin);
a(1)=aux(1);
a(12)=aux(2);
[aux,~,~]=Spin1Jcal0_1(m_q,spin);
a(2)=aux(1);
[aux,~,~]=Spin1Jcal1_2(m_q,spin);
a(3)=aux(1);
[aux,~,~]=Spin1Jcal2_1(m_q,spin);
a(4)=aux(1);
[aux,~,~]=Spin1Jcal0_2(m_q,spin);
a(5)=aux(1);
[aux,~,~]=Spin1Jcal2_2(m_q,spin);
a(6)=aux(1);
%[aux,W,x]=Spin1Jcal2_2(m_q,spin);
a(7)=aux(2);
[aux,~,~]=QuarkoniumS0J0(m_q,spin);
a(8)=aux(1);
[aux,~,~]=Spin1Jcal1_1(m_q,spin);
a(9)=aux(1);
%[aux,W,x]=Spin1Jcal1_1(m_q,spin);
a(10)=aux(2);
%[aux,W,x]=Spin1Jcal1_1(m_q,spin);
if m_q==1.4702
    a(11)=aux(3); %If charmonium
end
if m_q==4.8802
    a(11)=aux(5); %If bottomium
end
%[aux,W,x]=QuarkoniumS0J1(m_q,spin);
%a(12)=aux(2);
[aux,~,~]=QuarkoniumS0J2(m_q,spin);
a(13)=aux(1);
[aux,W,x]=Spin1Jcal3_1(m_q,spin);
a(14)=aux(1);

if m_q==1.4702
    chi= sum( ( (a - t) ./ e ).^2 );
end

%Final value ordered:
%(s/d)1
res(1,1)=a(1,1);
res(2,1)=a(1,2);
res(3,1)=a(1,3);
res(4,1)=a(1,4);
%(p1)
res(5,1)=a(1,12);
res(6,1)=a(1,5);
res(7,1)=a(1,9);
res(8,1)=a(1,6);
%(p/f)2
res(9,1)=a(1,13);
res(10,1)=a(1,10);
res(11,1)=a(1,7);
res(12,1)=a(1,14);
%p0
res(13,1)=a(1,8);
res(14,1)=a(1,11);





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
