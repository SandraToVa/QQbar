%Valor r0
setr0(3.964)
%Constants ajustades
setk1(-0.019)
setk2(0.0028)
%Valors Lambda g
setL1(-0.059)
setL3(-0.230)
load("dades.mat","m_c","m_b")
setm_q(m_b)
% l= interpolació (0), llagures distncies (1), bad long distances (2)
setl(0)
setspin(1) %s=1 for hybrids


%Programa que calcula los valors del espectre per als nivells del c
% 2(s/d)2 and 1d2
c=zeros(8,1);
b=zeros(8,1);

%Vector en los valors de la energia que necesito
[aux,~,~]=QuarkoniumS0J1(m_q,spin);
c(1)=aux(3);
if m_q==4.8802
    b(1)=aux(4);
end
[aux,~,~]=Spin1Jcal0_1(m_q,spin);
c(2)=aux(2);
[aux,~,~]=Spin1Jcal1_2(m_q,spin);
c(3)=aux(3);
c(6)=aux(2);
[aux,~,~]=Spin1Jcal2_1(m_q,spin);
c(4)=aux(3);
c(7)=aux(2);
[aux,~,~]=QuarkoniumS0J2(m_q,spin);
c(5)=aux(2);
if m_q==4.8802
    b(5)=aux(3);
end
%[aux,~,~]=Spin1Jcal1_2(m_q,spin);
%c(6)=aux(2);
%[aux,~,~]=Spin1Jcal2_1(m_q,spin);
%c(7)=aux(2);
[aux,~,~]=Spin1Jcal3_2(m_q,spin);
c(8)=aux(1);

c;



%Per al bottom son los estats c + los estats b

if m_q==4.8802
    %[aux,~,~]=QuarkoniumS0J1(m_q,spin);
    %b(1)=aux(4);
    [aux,~,~]=Spin1Jcal0_2(m_q,spin);
    b(2)=aux(2);
    [aux,~,~]=Spin1Jcal1_1(m_q,spin);
    b(3)=aux(3);
    b(6)=aux(4);
    [aux,~,~]=Spin1Jcal2_2(m_q,spin);
    b(4)=aux(3);
    b(7)=aux(5);
    %[aux,~,~]=QuarkoniumS0J2(m_q,spin);
    %b(5)=aux(3);
    %[aux,~,~]=Spin1Jcal1_1(m_q,spin);
    %b(6)=aux(4);
    %[aux,~,~]=Spin1Jcal2_2(m_q,spin);
    %b(7)=aux(5);
    [aux,~,~]=Spin1Jcal3_1(m_q,spin);
    b(8)=aux(3);

    b;
    
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

% VALOR DEL SPIN
function x7=spin
global v7
x7=v7;
end 

function setspin(val7)
global v7
v7 = val7;
end
