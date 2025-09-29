function [E,W,x]=Spin1Jcal1_2(m,s)
% Matriu 7x7 (capsa 2: 3x3) Acoblament spin per a s=1 j=1. Autor: Sandra Tomàs (Font: Ruben Oncala)
% Estats n(s/d)_J nd_J etc
% INCORPORACIÓ DELS POTENCIALS QUE DEPENEN DE L'SPIN EN ELS HYBRIDS
% No incorpora mixing. Mirar QuarkoinumAndHybrid.m si es vol
% 
%disp('**************************************************')
%disp('Hybrid-Quarkonium meson states with spin coupling s=1(j=1):')
%disp('**************************************************')
%format long;


  % the endpoints of the integration interval:
system.a=0.001;   
system.b=22;



% parameters of the boundary conditions:
system.A1= eye(3);
system.A2= zeros(3,3);
system.B1= eye(3);
system.B2= zeros(3,3);

% function handle to the function returning the potential matrix
system.V=@potentialMatrix;


%disp(' ');
%disp('Calculating eigenvalues with indices between 0 and 10:');
%transicions 0-5
t=cputime;
tol=3e-9; %x=314
[EigvData,meshData]=computeEigenvalues(system,0,4,tol); 

% disp(['Number of intervals in the mesh: ' num2str(length(meshData.h))]); % number of intervals in the mesh
%[lam0]=parameters2
% disp(['Time: ' num2str (cputime - t) ' sec']);
% disp(sprintf('k \t  E_k \t\t\t\t ErrEst \t status'))

E=EigvData.eigenvalues/m;

for i=1:length(EigvData.eigenvalues)
    %disp(sprintf('%-3.0f\t %-16.12f\t %-2.0e\t\t % d', EigvData.indices(i),EigvData.eigenvalues(i)/m,EigvData.errorEstimations(i),EigvData.status(i)));
    [x,Y,YP]=computeEigenfunction(system,meshData,EigvData.eigenvalues(i),1);
    W(:,:,i)=Y;
end


end

% aqui escollirem els parametres L i m i S del sistema ---------------------------------------------

function [m]=parameters
 m=m_q;
end

% TOTAL ANGULAR MOMENTUM (Jcal=J EN LO MIXING DESACTIVAT)
function [j]=parameters1
j=1;
end 

% Escollim cas quarkonium o cas hybrid quarkonium
function [s]=parameters3
 s=spin;     
%s=0      %Quarkonium
%s=1;     %Hybrid
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


%function [lam0,lam1,lam3,k,pi]=parameters2
% PERCENTATGE DE HYBRID (lam0=0 desactivem mixing)
%lam0=0.6;
%lam1=0.1182;
%lam3=0.2299;
%k=0.187;
%pi=3.1416;
%end
%-----------------------------------------------------------------------------------------------



function r=potentialMatrix(x) % returns the potential matrix evaluated in x

[j]=parameters1;
[s]=parameters3;

  for i=1:4 

       %El (-2) és el que es factor comú en la eq (4) del paper
  %Les matrius (8) i (10) son prporcionals a aquest -2Vhf
  %Les matrius (9) i (11) son prporcionals a +2Vhf2 per tant en les
  %expressions de les matrius hi ha un - extra que compensa al (-2)
  
  if s==1 % CAS HYBRIT QUARKONIUM P^(+,0,-)
    if j==1
      r(1,1,i)=A(x(i),j)+ (-2)*(-DD(x(i),j)+2*E9(x(i),j));
      r(2,2,i)=CC(x(i),j)+ (-2)*(+BB(x(i),j)+2*E1(x(i),j));
      r(3,3,i)=F(x(i),j+1)+ (-2)*(-BB(x(i),j)+2*A5(x(i),j));

      r(1,2,i)=B(x(i),j)+ (-2)*(+E3(x(i),j)+E7(x(i),j));
      r(2,1,i)=r(1,2,i);
      r(2,3,i)=+(-2)*(BC(x(i),j)+B4(x(i),j)+D2(x(i),j));
      r(3,2,i)=r(2,3,i);
            
      r(1,3,i)=+(-2)*(B6(x(i),j));
      r(3,1,i)=r(1,3,i);        
              
         
    end

  end
  

  end
end

  
    %---------------------------------------------------------------------------------------------
function f1=Vg(x)
  [m]=parameters;
  if m==1.4702
  Eo=2.6984;
  else
  Eo=9.5325;    
  end
  f1=-0.489/x+0.187*x+Eo; 
  end

% EQUACIÓ 8
  function f2=Vs(x)
  [m]=parameters;
  if m==1.4702
  Eo=3.499539;
  else
  Eo=10.333696;    
  end 
  f2=0.0611/x+0.187*x+Eo;
  end
%EQUACIÓ 11
  function f3=Vp(x)
  [m]=parameters;
  if m==1.4702
  Eo=3.491169;
  else
  Eo=10.325269;    
  end 
  b1=0.0696430221609656;
  b2=-1.4593432845775876;
  a1=-0.06732994686962318;
  a2=0.014330609468130364;
  f3=0.187*x+(0.0611/x)*(1.0+b1*x+b2*x.^2)/(1.0+a1*x+a2*x.^2)+Eo;
  end

  function f4=Vq(x)
  f4=Vp(x)-Vs(x);
  end
  
  function M1=F(x,j)
  [m]=parameters;
  M1=(j*(j+1))/(x.^2)+m*Vp(x);
  end

  function M2=R(x,j)
  [m]=parameters;
  M2=(j*(j+1))/(x.^2)+m*Vg(x);
  end

  function M3=A(x,j)
  [m]=parameters;
  M3=(j*(j-1))/(x.^2)+m*Vq(x)*(j+1)/(2*j+1)+m*Vs(x);
  end

  function M4=B(x,j)
  [m]=parameters;
  M4=m*Vq(x)*sqrt (j*(j+1))/(2*j+1);
  end

  function M5=CC(x,j)
  [m]=parameters;
  M5=((j+1)*(j+2))/(x.^2)+m*Vq(x)*(j)/(2*j+1)+m*Vs(x);
  end

  function M6=G(x,j)
  [m]=parameters;
  M6=2/(x.^2)+m*Vs(x);
  end

  %----------------------------------------------------------------

  
 % EFECTES HIPERFINS:
  %Corregits!!
  function HFP=Vsa(~) 
   % Potencial V_pp ---- Este és el Vsa no el Vpp
  [m]=parameters;
  %El signe - és el que porta la pròpia formula de Vsa
  %El signe que surti del L3 serà el signe del +- del propi glambda'''
  if m==1.4702
  cf=1.12155;
  k=0.187;
  HFP=-2*(cf*L3*pi^2)/(m*k*r0^3);
  else
  cf=0.87897;
  k=0.187;
  HFP=-2*(cf*L3*pi^2)/(m*k*r0^3);
  end
  
  end
  function HFS=Vsb(~) 
   % Potencial Vsp ---- Este és el Vsb no el Vsp
  [m]=parameters;
  %EL signe - és el signe que porta L1 (glambda') que sabem que es negatiu
  %El signe que tingue L1 serà el del +- del potencial Vsb
  %El dos del error està inclos aquí
  if m==1.4702
  cf=1.12155;
  k=0.187;
  HFS=-2*(cf*L1*pi^2)/(m*sqrt(pi*k)*r0^2);
  else
  cf=0.87897;
  k=0.187;
  HFS=-2*(cf*L1*pi^2)/(m*sqrt(pi*k)*r0^2);
  end
  
  end
  
  % PRIMER POTENCIAL DE SPIN ----------------------------------------------
 
  function HF1=Vhf(x) 
   % Control de l'efecte del 1r potencial
   % El -2 que acompanya als potencials ja etsà inclos en la matriu de
   % sobre. Aquí només la pròpia expressió del Vhf
   if l==0 %interpolació
       HF1=(k1+((x./r0).^2*((Vsa/6)-(1/3)*(x./r0)*Vsb)))./(1+(x./r0).^5);
   elseif l==1 %llargues distàncies
       hevi=heaviside(sym(x) - r0);
       HF1=((Vsa/6)*(r0^3/x.^3)-(1/3)*Vsb*(r0^2/x.^2))*hevi;
   elseif l==2 %Bad long distances:
       HF1=(Vsa/6)-(1/3)*Vsb;
   end
  end  
  
  function H1=BB(x,j)
  H1=Vhf(x)/(j+1);
  end
  function H2=BC(x,j)
  H2=Vhf(x)*((sqrt(j*(j+2)))/(j+1));
end
function H3=DD(x,j)
  H3=Vhf(x)/j;
  end
  function H4=DE(x,j)
  H4=Vhf(x)*(sqrt((j^2)-1)/j);
  end
  function H5=FF(x,j)
  H5=Vhf(x)*((j-1)/j);
  end
  function H6=FG(x,j)
  H6=Vhf(x)*(((j+1)*(sqrt(2*j-1)))/(j*(sqrt(2*j+1))));
  end
  function H7=GG(x,j)
  H7=Vhf(x)/(j*(j+1));
end
function H8=GH(x,j)
  H8=Vhf(x)*((j*sqrt(2*j+3))/((j+1)*sqrt(2*j+1)));
end  
function H9=HH(x,j)
  H9=Vhf(x)*((j+2)/(j+1));
end

% SEGON POTENCIAL DE SPIN -------------------------------------------------
 function HF2=Vhf2(x) 
   % Control de l'efecte del 1r potencial
   if l==0 %interpolació
       HF2=((k2*(x.^2)) - ( (x./r0).^5 * ( (1/2)*(r0./x)*Vsa + (1/2)*Vsb ) )) ./ (1+(x./r0).^7);
   elseif l==1 %llargues distàncies
       hevi=heaviside(sym(x) - r0);
       HF2=-(1/2)*( Vsa*(r0^3/x.^3) + Vsb*(r0^2/x.^2) )*hevi;
   elseif l==2 %Bad long distances:
       HF2=-(1/2)*( Vsa + Vsb );
   end
  end 
  
  function AA1=A1(x,j)
  AA1=Vhf2(x)*(((3+j)/(9+6*j)));
 end
   function AA3=A3(x,j)
  AA3=-Vhf2(x)*((sqrt((1+j)*(2+j)))/(3+2*j));
   end
   function AA5=A5(x,j)
  AA5=Vhf2(x)/(3+3*j);
   end
   function AA7=A7(x,j)
  AA7=Vhf2(x)*(((2+j)^(3/2))/(sqrt(1+j)*(3+2*j)));
   end
   function AA9=A9(x,j)
  AA9=-Vhf2(x)*((j*(2+j))/(9+15*j+6*j^2));
   end
   function BB4=B4(x,j)
  BB4=Vhf2(x)*((sqrt(j)*((2+j)^(3/2)))/(3*(1+j)*(1+2*j)));
   end
   function BB6=B6(x,j)
  BB6=-Vhf2(x)*((j*sqrt(1+1/(1+j)))/(1+2*j));
   end
   function BB8=B8(x,j)
  BB8=-Vhf2(x)*(j/(3*(j+1)))*(((2*j+3)/(2*j+1))^(1/2));
   end
   function DD2=D2(x,j)
  DD2=-Vhf2(x)*((sqrt(j*(2+j)))/(3+3*j));
  end
   function DD4=D4(x,j)
  DD4=-Vhf2(x)*(j*(2+j))/(sqrt((1+j)*(2+j)*(1+2*j)*(3+2*j)));
  end
  function DD6=D6(x,j)
  DD6=Vhf2(x)*(j^2)/(3*(1+j)*sqrt(3+4*j*(2+j)));
   end
   function EE1=E1(x,j)
  EE1=Vhf2(x)*(2+j)/(3+9*j+6*(j^2));
   end
   function EE3=E3(x,j)
  EE3=-Vhf2(x)*(sqrt(j/(1+j)))/(1+2*j);
   end
   function EE5=E5(x,j)
  EE5=Vhf2(x)/(3*j*(1+j));
  end
   function EE7=E7(x,j)
  EE7=Vhf2(x)*(sqrt(1+1/j))/(1+2*j);
   end
  function EE9=E9(x,j)
  EE9=Vhf2(x)/(3*j)-Vhf2(x)/(1+2*j);
   end
   function FF4=F4(x,j)
  FF4=Vhf2(x)*((1+j)^2)/(3*j*(sqrt(4*(j^2)-1)));
   end
   function FF6=F6(x,j)
  FF6=-Vhf2(x)*(j^2-1)/sqrt(j*(j-1)*(2*j-1)*(2*j+1));
  end
   function FF8=F8(x,j)
  FF8=-Vhf2(x)*(sqrt(j^2-1))/(3*j);
   end
   function GG2=G2(x,j)
  GG2=-Vhf2(x)*((1+j)*sqrt(1-2/(1+2*j)))/(3*j);
   end
   function GG4=G4(x,j)
  GG4=-Vhf2(x)*(sqrt((j-1)/j)*(1+j))/(1+2*j);
   end
   function GG6=G6(x,j)
  GG6=Vhf2(x)*(sqrt((j-1)/(j+1))*(j^2-1))/(3*j*(1+2*j));
  end
   function II1=I1(x,j)
  II1=-Vhf2(x)*(1-j^2)/(3*j-6*j^2);
   end
  function II3=I3(x,j)
  II3=Vhf2(x)*((j-1)^(3/2))/(sqrt(j)*(2*j-1));
   end
   function II5=I5(x,j)
  II5=-Vhf2(x)/(3*j);
   end
   function II7=I7(x,j)
  II7=-Vhf2(x)*sqrt((j-1)*j)/(2*j-1);
  end
   function II9=I9(x,j)
  II9=-Vhf2(x)*(2-j)/(-3+6*j);
  end