%Code to compute the expected value of r, 1/r, 1/r^2 etc


% One of its uses is to check if for a transition from híbrid to quarkonium the value of
 %<Psi_QQ| 1/r |Psi_HQ> with r = x in the code is >> than DeltaE of the QQ
 %and HQ states.

 %We use hibrid states but not splited to the hiperfine so we only need r0,
 %mq and spin

setr0(3.964)
load("dades.mat","m_c","m_b")
setm_q(m_c)

%I will try 2 transitions 1(s/d)->3s and 1(s/d)->2s for bottomonium
% Initialize cell array to store variable-sized matrices
wf_hq = cell(1, 1);         %for each cell there will be a matrix(wf made by diferent states at diferent x positions)
wf_qq = cell(1, 2); 
xvector = cell(1, 2);

E = zeros(2,2);
Rvalues = zeros(1,2);
R2values = zeros(1,2);
DEvalues = zeros(1,2);
DERvalues = zeros(1,2);


%We compute the hibrid functions
setspin(1)

%Hibrid states without the hiperfine:
[aux,auxwf,x]=QuarkoniumS0J1(m_q,spin);
% 1--
E(1,1)=aux(2);
E(1,2)=aux(2);
wf_hq{1} = auxwf(1, :, 1); %Primer element 1:2 si son (s/d) o (p/f)
xvector{1}=x;

%We compute the quarkonium functions
setspin(0)

%Hibrid states without the hiperfine:
[aux,auxwf,x]=QuarkoniumS0J1(m_q,spin);
% 1--
E(2,1)=aux(2);
E(2,2)=aux(3);
wf_qq{1} = auxwf(1, :, 2);
wf_qq{2} = auxwf(1, :, 3);
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
Rmatrix = compute1Rmatrix(xvector{1});
R2matrix = compute2Rmatrix(xvector{1});
%The value of the 1/r sandwitch
Rvalues(1) = computeSandwitch(Rmatrix,wf_hq{1},wf_qq{1});
Rvalues(2) = computeSandwitch(Rmatrix,wf_hq{1},wf_qq{2});
%The value of the 1/r^2 sandwitch
R2values(1) = computeSandwitch(R2matrix,wf_hq{1},wf_qq{1});
R2values(2) = computeSandwitch(R2matrix,wf_hq{1},wf_qq{2});
%Values of the diference of energy
DEvalues(1) = E(2,1)-E(1,1);
DEvalues(2) = E(2,2)-E(1,2);
%Values of the diference of energy/r
DERvalues(1) = DEvalues(1)*Rvalues(1);
DERvalues(2) = DEvalues(2)*Rvalues(2);



%The result should be the same as using with Op the Rmatirx 3 times
%Op=zeros(lengthXhq,lengthXhq,3);
%Op(:,:,1)=Rmatrix;
%Op(:,:,2)=Rmatrix;
%R2sdto1p=ExpectedValue(3,1,Op,@QuarkoniumS0J1,@QuarkoniumS0J1,m_q,1,0);
%R2sdto2p=ExpectedValue(3,2,Op,@QuarkoniumS0J1,@QuarkoniumS0J1,m_q,1,0);


% Computes the 1/r matrix for diferent híbrid states

 function Rmatrix=compute1Rmatrix(x)
%Input:     dim = dimension of the rellevant! wf =1 for p_n and d_n states
%                 and 2 for (s/d) and (p/f)
%           x = array of meshpoints in where the matrix are computed

%Output:    Rmatrix = matrix corresponding to 1/r that will be multiplied only for the
%                       rellevant wf

    Rvector = zeros(length(x),length(x));
    for x_el=1:length(x)
        Rvector(x_el,x_el) = 1 / x(x_el);
    end 
    Rmatrix = Rvector;
 end

function R2matrix=compute2Rmatrix(x)
%Input:     dim = dimension of the rellevant! wf =1 for p_n and d_n states
%                 and 2 for (s/d) and (p/f)
%           x = array of meshpoints in where the matrix are computed

%Output:    Rmatrix = matrix corresponding to 1/r that will be multiplied only for the
%                       rellevant wf

    Rvector = zeros(length(x),length(x));
    for x_el=1:length(x)
        Rvector(x_el,x_el) = 1 / x(x_el).^2;
    end 
    R2matrix = Rvector;
 end

 function Rvalue=computeSandwitch(Rmatrix,wf_hq,wf_qq)
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

        Yf_f=wf_qq; %row vector fianl
        Yi_f=wf_hq(N,:); %row vector inicial
        Yi_c = Yi_f';

        result_matrix =  Yf_f * Rmatrix * Yi_c;
        norm = Yf_f *  Yi_c;

        Rvector(N) = result_matrix / norm;

        Rvalue = Rvalue + Rvector(N);

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