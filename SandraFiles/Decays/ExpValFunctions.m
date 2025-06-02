%Different functions for QQ and HQ and HH because the operators i diferent
%For Quarkonium-Quarkonium the operators used is one because the wave
%functions are all 1 row vectors
%For Hybrid-Quarkonium the operators are 2 because the hybrids have 2 rows
%and for each row the operator is diferent so we need to import 2 operators


function res = ExpValFunctions(x)
    if strcmp(x, 'QQ')
        res = @QQComputeExpValM;
    elseif strcmp(x, 'HQ')
        res = @HQComputeExpValM;
    end
end

%Prepares the elements
function [EV, DeltaE]=QQComputeExpValM(i,f,M,trans,Ji,Jf)
%With a initial al final states we compute expectations values
%The operators of the expectation depend on M (dipion mass)
%We want to find the relationship between the result and the M used
%trans=transition function StoS, PtoP, etc

%initail
setspin(0);
sin=spin;
sfin=spin;

%Define the funtions to use
QuarkoniumS0Ji=str2func(['QuarkoniumS0J' num2str(Ji)]);
QuarkoniumS0Jf=str2func(['QuarkoniumS0J' num2str(Jf)]);

[Einicial,~,xi]=QuarkoniumS0Ji(m_q,sin);

Ei=Einicial(i);

%final
[Efinal,~,xf]=QuarkoniumS0Jf(m_q,sfin);

Ef=Efinal(f);

DeltaE=abs(Ef-Ei);
if M>DeltaE
    EV=0;
    disp(['M>DeltaE it cannot be', num2str(DeltaE)])

else
%Operator that we want to find the expectated value
%The x must be the same. It is a problem of tolerance if not
    if length(xi)~=length(xf)
     disp('Change tolerance:');
     disp([length(xi), length(xf), length(xi)-length(xf)]);
    else
     x=xi; %The same for xi and xf same x
     Op=trans(x,Ei,Ef,M);
     Op=diag(diag(Op));

     EV=ExpectedValue(i,f,Op,QuarkoniumS0Ji,QuarkoniumS0Jf,m_q,sin,sfin);

    end
end

end




%IMPORTANT: I HAVE TO REWRITE THIS

function EV=HQComputeExpValM(i,f,M,trans,Ji,Jf,double)
%With a initial al final states we compute expectations values
%The operators of the expectation depend on M (dipion mass)
%We want to find the relationship between the result and the M used
%trans=transition function from I_thetaFunctions respository
%double= if the initial state (the hybrid) is a state composed by two
%(s/d)1 or (p/f)2 double is true and the operator has 2 components and the
%w.f. is not only one like but two

%Define the funtions to use
QuarkoniumS0Ji=str2func(['QuarkoniumS0J' num2str(Ji)]);
QuarkoniumS0Jf=str2func(['QuarkoniumS0J' num2str(Jf)]);

%initail
setspin(1);
sin=spin;

[Einicial,~,xi]=QuarkoniumS0Ji(m_q,sin);

Ei=Einicial(i);

%final
setspin(0)
sfin=spin;

[Efinal,~,xf]=QuarkoniumS0Jf(m_q,sfin);

Ef=Efinal(f);

E1=abs(Ef-Ei);
if M>E1
    EV=0;
    disp(['M>E1 cannot be', num2str(E1)])

else
%Operator that we want to find the expectated value
%The x must be the same. It is a problem of tolerance if not
    if length(xi)~=length(xf)
     disp('Change tolerance:');
     disp(length(xi)-length(xf));
    else
    %We must direnciate for hybrids with 2component wf and hybrids with
    %1component wf
    %create a variable: "double". If true, it has double wf. If false it has
    %not
     x=xi; %The same for xi and xf same x

     if double
         %Double hybrids have 3 rows and the 1st and 2nd are treated
         %searately and the 3rd is zeros

         %Important: The 1st row is the one that will multiply the 1st row
         %of the wf. They are s from (s/d)1 and p from (p/f)2
         Op_d=trans(x,Ei,Ef,M);  %double
         Op_d1=Op_d(:,:,1);
         Op_d2=Op_d(:,:,2);
         [numR,numC] = size(Op_d1);

         Op=zeros(numR,numC,3);

         Op(:,:,1)=diag(diag(Op_d1));
         Op(:,:,2)=diag(diag(Op_d2));
         Op(:,:,3)=zeros(numR,numC);


     else
         %Non double hybrids have 3 rows but only the third/first has
         %information (or one or the other)
         Op_s=trans(x,Ei,Ef,M);   %single
         [numR,numC] = size(Op_s);
         Op_s=diag(diag(Op_s));

         Op=zeros(numR,numC,3);

         Op(:,:,1)=Op_s;
         Op(:,:,2)=zeros(numR,numC);
         Op(:,:,3)=Op_s;
     end

     %PROBLEM WITH THE FUNCTION, x for hybrid and x for quarkonium have diferent
     %dimension. Due to ConstrucMesh?????
     V=ExpectedValue(i,f,Op,QuarkoniumS0Ji,QuarkoniumS0Jf,m_q,sin,sfin);

     %SHA DE CANVIAR PER A FER LA SUMA SI DOUBLE I NO SI SINGLE!! o meibi no
     %fa falta
     EV=sum(V);

    end
end

end


function V=ExpectedValue(i,f,Op,initial,final,m,sin,sfin)
%#codegen

%Function to compute espected values using wavefunctions and any observable
%Output:        - Returns the expectation value of
%               final wavefunction * operator * initial wavefunction
%Input:         - i=position in the spectum of eneries of the hybrid or
%               quarkonium that we want the inital wavefunction to have
%               - f= same as i but for final wavefucntion
%               - Op= operator in matrix form (it has to be a diagonal
%               matrix) the rowas and columns the same as the
%               wavefunctions. This is the operator that we want to find
%               the expectation value of
%               - initial= initial wavefunction form the scripts of "Spin..."
%               for hybrid and "Quarkonium..." for quarkonium
%               - final= same as initial for final state
%               - m= mass of the quark. Charm or Bottom

%I have to modify the tolerance in "initial" and "final" for the mesh to
%coincide. And for the x and lenght of the 2nd dimension of W to coincide!
setspin(sin)
[Ei,Wi,~]=initial(m,spin);

setspin(sfin)
[Ef,Wf,x]=final(m,spin);

%The vale of our initial state
%disp('Initial state i')
%i;
%disp('With energy:')
Ei(i);
%disp('Wavefunction:')
Yi=Wi(:,:,i);

%The length of the wavefunction i 
[numRi,~] = size(Yi);


%The vale of our final state
%disp('Initial final f')
%f;
%disp('With energy:')
Ef(f);
%disp('Wavefunction:')
Yf=Wf(:,:,f);

%The length of the wavefunction f 

[numRf,~] = size(Yf);

if numRi ~= numRf
    disp('Number of rows in wave functions not equal')
end

%If we work with quarkonium the wave functions only have on relevant row:
%the first one
if sin==0 && sfin ==0
    Yi_f=Yi(1,:); %row vector inicial
    Yi_c=Yi_f'; %column vector inicial
    
    Yf_f=Yf(1,:); %row vector fianl

    norm = Yf_f * Yi_c;
    
    %sandwitch wave function en operador
    result = Yf_f * Op * Yi_c; 

    V = result/norm;

end %end of if
        
%If I work with hybrid-to-qurakonium, the hybrid wave functions have 2 rows
%for each row I have to apply different operators between the quarkonium
%and the bybrid wavefunctions

if sin==1 && sfin==0
    %If the initial states is hybrid, (should be) numRi=3 if the hybrid is
    %a double the two first rows are treated separately and if the hybrid
    %is not a double the only relevant row is the last one
    %If the final state is quarkonium, (should be) numRi=3 but we only
    %want the first one
        
    Vvector=zeros(numRi,1);
      
    for N=1:numRi
        %Muliply each row for any row of the other wavefunction normalized
        
        Yi_f=Yi(N,:); %row vector inicial
        Yi_c=Yi_f'; %column vector inicial

        Yf_f=Yf(1,:); %row vector final

        %The operator Op in this case is three elements, one for the upper row and
        %one for the bottom row (in the case of double) and the third for
        %the case of single
        OpRow=Op(:,:,N);

        %compute the norm for each row of the hibrid state
        norm = Yf_f * Yi_c;

        %sandwitch wave function with operator for each row of the hibrid
        result = Yf_f * OpRow * Yi_c; 

        Vvector(N,1) = result/norm; %store for each row

     end %end for

     %to get a unic value ignoring the NaN because of a function having the
     %1st and 2nd or the last component 3 depending on whether they are
     %(s/d) or p0 we do

     %V=nansum(Vector)

     V=Vvector;
    
end %end if


end %end function


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

