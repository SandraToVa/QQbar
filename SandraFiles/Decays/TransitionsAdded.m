
%IMPORTANT: This code is updated (only the quarkonium part) with the 
% inclusion of the new operator


%As I don't know the best way of making this type of code. This is simply a
%document that acts as a gide of what I have to use inside the loop of
%ComputationExpVal.m

% Uses the functions from FormFactor_ItoF.m and gives the final function
% to compute the transition
%Example: in l=1->l'=1 there are multiple transitions from m->m' that apear
%in I_thetaFuntions. Here we add and obtain the final transitions l->l' in
%terms of funcions m->m'

%This functions will be called in ComputationExpVal in order to directly
%compute the I_ItoF (\, c, s) for each M

%Create a function to call all the functions in here:
function res = TransitionsAdded(x)
    disp(['Input value: ', x]);

    %%%%%%%%%% QQ %%%%%%%%%%%%
    %s->s
    if strcmp(x, 'QQStoS')
        res = @StoStrans;
    %p->p
    elseif strcmp(x, 'QQPtoP')
        res = @PtoPtrans;
    %d->s
    elseif strcmp(x, 'QQDtoS')
        res = @DtoStrans;
    %s->d
    elseif strcmp(x, 'QQStoD')
        res = @StoDtrans;
    %d->d
    elseif strcmp(x, 'QQDtoD')
        res = @DtoDtrans;
   
    %%%%%%%%%%%% HQ %%%%%%%%%%%%%%
    % p1 -> s
    elseif strcmp(x, 'HQp1tS') 
        res = @P1toStrans;
     % p1 -> d
    elseif strcmp(x, 'HQp1tD') 
        res = @P1toDtrans;


    % (s/d)1 -> p
    elseif strcmp(x, 'HQsdtP') 
        res = @SD1toPtrans;

    end
    disp(['Output function handle: ', func2str(res)]);
end

%First all the functions in Quarkonium -> Quarkonium +2Pions transitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% s->s
function [I_if_square_cell, DeltaE, M] = StoStrans(Nin, Nfin, Min, Mfin, length)

k = length;

% Dipion mass row
nM = (Mfin - Min) / k;
M = Min:nM:Mfin;

% Initialize variables
I_if_square_cell = cell(1, 4);
%I_if_square = zeros(1, k + 1);
%I_if_c_square = zeros(1, k + 1);
%I_if_0c_square = zeros(1, k + 1);
I_if_s_square = zeros(1, k + 1);
%m=0
I_if0 = zeros(1, k + 1);
I_if_c0 = zeros(1, k + 1);

E = zeros(1, k + 1);

% Precompute ExpValFunc
%només 1 contribució Ji=0,m=0->Jf=0,m'=0
%al spin average contribueixen els termes I^2, I_c^2 i I*I_c
ExpValFunc = ExpValFunctions('QQ');
transition0 = FormFactor_ItoF('QQS0toS0_F/');
transitionC0 = FormFactor_ItoF('QQS0toS0_Fc');

% Compute first E value and check if consistent
[I_if0(1), E(1)] = ExpValFunc(Nin, Nfin, M(1), transition0, 0, 0);
[I_if_c0(1), ~] = ExpValFunc(Nin, Nfin, M(1), transitionC0, 0, 0);

% Fill the rest of the arrays
for n = 2:k+1
    [I_if0(n), E(n)] = ExpValFunc(Nin, Nfin, M(n), transition0, 0, 0);
    [I_if_c0(n), ~] = ExpValFunc(Nin, Nfin, M(n), transitionC0, 0, 0);

    % Check if E is constant
    if E(n) ~= E(1)
        error('E values are not consistent across iterations.');
    end
end

% Compute spin average = com nomes hi ha 1 m possible per a la transicio
%directament per a cada
I_if_square = I_if0.^2;
I_if_c_square = I_if_c0.^2;
I_if_0c_square = I_if0 .* I_if_c0;

% Store in cell array
I_if_square_cell{1, 1} = I_if_square;
I_if_square_cell{1, 2} = I_if_c_square;
I_if_square_cell{1, 3} = I_if_0c_square;
I_if_square_cell{1, 4} = I_if_s_square;

% Return E as a single scalar if consistent
DeltaE = E(1);

end


% p->p
function [I_if_square_cell, DeltaE, M] = PtoPtrans(Nin, Nfin, Min, Mfin, length)

k = length;

% Dipion mass row
nM = (Mfin - Min) / k;
M = Min:nM:Mfin;

% Initialize variables
I_if_square_cell = cell(1, 4);
%I_if_square = zeros(1, k + 1);
%I_if_c_square = zeros(1, k + 1);
%I_if_0c_square = zeros(1, k + 1);
%I_if_s_square = zeros(1, k + 1);
%m=0
I_if0 = zeros(1, k + 1);
I_if_c0 = zeros(1, k + 1);
%m=1
I_if1=zeros(1,k+1);
I_if_c1=zeros(1,k+1);
I_if_s1=zeros(1,k+1);

E = zeros(1, k + 1);

% Precompute ExpValFunc
% 3 contributions: m=0->m'=0, m=1->m'=1 (estes dos son la mateixa)
% m=-1->m'=-1 for the I, Ic form factors. We will ned a spin average
% (sumatori sobre estats finals i average sobre inicals)

% 1 contribution: m=+1->m'=-1 for the Ix and the same m=-1->m'=+1 for the
% Is case
ExpValFunc = ExpValFunctions('QQ');
%m=0
transition0=FormFactor_ItoF('QQP0toP0_F/');
transitionC0=FormFactor_ItoF('QQP0toP0_Fc');
%m=1
transition1=FormFactor_ItoF('QQP1toP1_F/');
transitionC1=FormFactor_ItoF('QQP1toP1_Fc');
transitionS1=FormFactor_TtoF('QQP1toP1_Fs');


% Compute first E value and check if consistent
[I_if0(1), E(1)] = ExpValFunc(Nin, Nfin, M(1), transition0, 1, 1);
[I_if_c0(1), ~] = ExpValFunc(Nin, Nfin, M(1), transitionC0, 1, 1);
[I_if1(1), ~] = ExpValFunc(Nin, Nfin, M(1), transition1, 1, 1);
[I_if_c1(1), ~] = ExpValFunc(Nin, Nfin, M(1), transitionC1, 1, 1);
[I_if_s1(1), ~] = ExpValFunc(Nin, Nfin, M(1), transitionS1, 1, 1);


% Fill the rest of the arrays
for n = 2:k+1
    [I_if0(n), E(n)] = ExpValFunc(Nin, Nfin, M(n), transition0, 1, 1);
    [I_if_c0(n), ~] = ExpValFunc(Nin, Nfin, M(n), transitionC0, 1, 1);
    [I_if1(n), ~] = ExpValFunc(Nin, Nfin, M(n), transition1, 1, 1);
    [I_if_c1(n), ~] = ExpValFunc(Nin, Nfin, M(n), transitionC1, 1, 1);
    [I_if_s1(n), ~] = ExpValFunc(Nin, Nfin, M(n), transitionS1, 1, 1);

    % Check if E is constant
    if E(n) ~= E(1)
        error('E values are not consistent across iterations.');
    end
end

% Compute spin average = com nomes hi ha 1 m possible per a la transicio
%directament per a cada
I_if_square = (I_if0.^2 + 2*I_if1.^2) / 3; % *2 because trans from +-1 to +-1
I_if_c_square = (I_if_c0.^2 + 2*I_if_c1.^2) / 3;
I_if_0c_square = (I_if0 .* I_if_c0 + 2*I_if1 .* I_if_c1) / 3;
I_if_s_square = (2*I_if_s1.^2) / 3; % *2 because trans from +-1 to -+1

% Store in cell array
I_if_square_cell{1, 1} = I_if_square;
I_if_square_cell{1, 2} = I_if_c_square;
I_if_square_cell{1, 3} = I_if_0c_square;
I_if_square_cell{1, 4} = I_if_s_square;

% Return E as a single scalar if consistent
DeltaE = E(1);  

end


% d->s
function [I_if_square_cell, DeltaE, M] = DtoStrans(Nin, Nfin, Min, Mfin, length)

k = length;

% Dipion mass row
nM = (Mfin - Min) / k;
M = Min:nM:Mfin;

% Initialize variables
I_if_square_cell = cell(1, 4);
%I_if_square = zeros(1, k + 1);
%I_if_c_square = zeros(1, k + 1);
%I_if_0c_square = zeros(1, k + 1);
%I_if_s_square = zeros(1, k + 1);
%m=0
I_if0 = zeros(1, k + 1);
I_if_c0 = zeros(1, k + 1);
%m=2
I_if_s2 = zeros(1, k + 1);

E = zeros(1, k + 1);

% Precompute ExpValFunc
% 1 contribution: m=0->m'=0 for the I, Ic form factors
% 1 contribution: m=+2->m'=0 for the Ix and the same m=-2->m'=0 for the
% Is case
ExpValFunc = ExpValFunctions('QQ');
%m=0
transition0 = FormFactor_ItoF('QQD0toS0_F/');
transitionC0 = FormFactor_ItoF('QQD0toS0_Fc');
%m=2
transitionS2 = FormFactor_ItoF('QQD2toS0_Fs');

% Compute first E value and check if consistent
[I_if0(1), E(1)] = ExpValFunc(Nin, Nfin, M(1), transition0, 2, 0);
[I_if_c0(1), ~] = ExpValFunc(Nin, Nfin, M(1), transitionC0, 2, 0);
[I_if_s2(1), ~] = ExpValFunc(Nin, Nfin, M(1), transitionS2, 2, 0);

% Fill the rest of the arrays
for n = 2:k+1
    [I_if0(n), E(n)] = ExpValFunc(Nin, Nfin, M(n), transition0, 2, 0);
    [I_if_c0(n), ~] = ExpValFunc(Nin, Nfin, M(n), transitionC0, 2, 0);
    [I_if_s2(n), ~] = ExpValFunc(Nin, Nfin, M(n), transitionS2, 2, 0);


    % Check if E is constant
    if E(n) ~= E(1)
        error('E values are not consistent across iterations.');
    end
end


% Compute el spin average = com nomes hi ha 1 m possible per a la transicio
%directament per a cada
I_if_square = (I_if0.^2) / 5;
I_if_c_square = (I_if_c0.^2) / 5;
I_if_0c_square = (I_if0 .* I_if_c0) / 5;
I_if_s_square = (I_if_s2.^2) / 5;

% Store in cell array
I_if_square_cell{1, 1} = I_if_square;
I_if_square_cell{1, 2} = I_if_c_square;
I_if_square_cell{1, 3} = I_if_0c_square;
I_if_square_cell{1, 4} = I_if_s_square;

% Return E as a single scalar if consistent
DeltaE = E(1);  

end


% s->d
function [I_if_square_cell, DeltaE, M] = StoDtrans(Nin, Nfin, Min, Mfin, length)

k = length;

% Dipion mass row
nM = (Mfin - Min) / k;
M = Min:nM:Mfin;

% Initialize variables
I_if_square_cell = cell(1, 4);
%I_if_square = zeros(1, k + 1);
%I_if_c_square = zeros(1, k + 1);
%I_if_0c_square = zeros(1, k + 1);
%I_if_s_square = zeros(1, k + 1);
%m=0
I_if0 = zeros(1, k + 1);
I_if_c0 = zeros(1, k + 1);
%m=2
I_if_s2 = zeros(1, k + 1);

E = zeros(1, k + 1);

% Precompute ExpValFunc
% 1 contribution: m=0->m'=0 for the I, Ic form factors
% 1 contribution: m=+2->m'=0 for the Ix and the same m=-2->m'=0 for the
% Is case
ExpValFunc = ExpValFunctions('QQ');
%m=0
transition0 = FormFactor_ItoF('QQD0toS0_F/');
transitionC0 = FormFactor_ItoF('QQD0toS0_Fc');
%m=2
transitionS2 = FormFactor_ItoF('QQD2toS0_Fs');

% Compute first E value and check if consistent
[I_if0(1), E(1)] = ExpValFunc(Nin, Nfin, M(1), transition0, 0, 2);
[I_if_c0(1), ~] = ExpValFunc(Nin, Nfin, M(1), transitionC0, 0, 2);
[I_if_s2(1), ~] = ExpValFunc(Nin, Nfin, M(1), transitionS2, 0, 2);

% Fill the rest of the arrays
for n = 2:k+1
    [I_if0(n), E(n)] = ExpValFunc(Nin, Nfin, M(n), transition0, 0, 2);
    [I_if_c0(n), ~] = ExpValFunc(Nin, Nfin, M(n), transitionC0, 0, 2);
    [I_if_s2(n), ~] = ExpValFunc(Nin, Nfin, M(n), transitionS2, 0, 2);


    % Check if E is constant
    if E(n) ~= E(1)
        error('E values are not consistent across iterations.');
    end
end


% Compute el spin average = com nomes hi ha 1 m possible per a la transicio
%directament per a cada
%For spin =0
I_if_square = (I_if0.^2);
I_if_c_square = (I_if_c0.^2);
I_if_0c_square = (I_if0 .* I_if_c0);
I_if_s_square = (I_if_s2.^2);

%For spin = 1
%I_if_square = (2*(1/2)*I_if0.^2)/3;
%I_if_c_square = (I_if_c0.^2)/3;
%I_if_0c_square = (I_if0 .* I_if_c0)/3;
%I_if_s_square = ((2/3)*I_if_s2.^2 + (1/3)*I_if_s2.^2)/3;


% Store in cell array
I_if_square_cell{1, 1} = I_if_square;
I_if_square_cell{1, 2} = I_if_c_square;
I_if_square_cell{1, 3} = I_if_0c_square;
I_if_square_cell{1, 4} = I_if_s_square;

% Return E as a single scalar if consistent
DeltaE = E(1);  

end


% d->d
function [I_if_square_cell, DeltaE, M] = DtoDtrans(Nin, Nfin, Min, Mfin, length)

k = length;

% Dipion mass row
nM = (Mfin - Min) / k;
M = Min:nM:Mfin;

% Initialize variables
I_if_square_cell = cell(1, 4);
%I_if_square = zeros(1, k + 1);
%I_if_c_square = zeros(1, k + 1);
%I_if_0c_square = zeros(1, k + 1);
I_if_s_square = zeros(1, k + 1);
%m=0
I_if0=zeros(1,k+1);
I_if_c0=zeros(1,k+1);
%m=1
I_if1=zeros(1,k+1);
I_if_c1=zeros(1,k+1);
%I_if_s1=zeros(1,k+1);
%m=2
I_if2=zeros(1,k+1);
I_if_c2=zeros(1,k+1);
%I_if_s2=zeros(1,k+1);

E = zeros(1, k + 1);

% Precompute ExpValFunc
% 5 contributions: m=0->m'=0, m=1->m'=1 (estes dos son la mateixa)
% m=-1->m'=-1, m=+-2->m'=+-2 for the I, Ic form factors. We will ned a spin average
% (sumatori sobre estats finals i average sobre inicals)
ExpValFunc=ExpValFunctions('QQ');
%m=0
transition0=FormFactor_ItoF('QQD0toD0_F/');
transitionC0=FormFactor_ItoF('QQD0toD0_Fc');
%m=1
transition1=FormFactor_ItoF('QQD1toD1_F/');
transitionC1=FormFactor_ItoF('QQD1toD1_Fc');
%m=1
transition2=FormFactor_ItoF('QQD2toD2_F/');
transitionC2=FormFactor_ItoF('QQD2toD2_Fc');

% Compute first E value and check if consistent
[I_if0(1), E(1)] = ExpValFunc(Nin, Nfin, M(1), transition0, 2, 2);
[I_if_c0(1), ~] = ExpValFunc(Nin, Nfin, M(1), transitionC0, 2, 2);
[I_if1(1), ~] = ExpValFunc(Nin, Nfin, M(1), transition1, 2, 2);
[I_if_c1(1), ~] = ExpValFunc(Nin, Nfin, M(1), transitionC1, 2, 2);
[I_if2(1), ~] = ExpValFunc(Nin, Nfin, M(1), transition2, 2, 2);
[I_if_c2(1), ~] = ExpValFunc(Nin, Nfin, M(1), transitionC2, 2, 2);

% Fill the rest of the arrays
for n = 2:k+1
    [I_if0(n), E(n)] = ExpValFunc(Nin, Nfin, M(n), transition0, 2, 2);
    [I_if_c0(n), ~] = ExpValFunc(Nin, Nfin, M(n), transitionC0, 2, 2);
    [I_if1(n), ~] = ExpValFunc(Nin, Nfin, M(n), transition1, 2, 2);
    [I_if_c1(n), ~] = ExpValFunc(Nin, Nfin, M(n), transitionC1, 2, 2);
    [I_if2(n), ~] = ExpValFunc(Nin, Nfin, M(n), transition2, 2, 2);
    [I_if_c2(n), ~] = ExpValFunc(Nin, Nfin, M(n), transitionC2, 2, 2);

    % Check if E is constant
    if E(n) ~= E(1)
        error('E values are not consistent across iterations.');
    end
end

% Compute el spin average
I_if_square = (I_if0.^2 + 2*I_if1.^2 + 2*I_if2.^2) / 5; % *2 because trans from +-1 to +-1 and from +-2 to +-2
I_if_c_square = (I_if_c0.^2 + 2*I_if_c1.^2 + 2*I_if_c2.^2) / 5;
I_if_0c_square = (I_if0.*I_if_c0 + 2*I_if1.*I_if_c1 + 2*I_if2.*I_if_c2) / 5;

% Store in cell array
I_if_square_cell{1, 1} = I_if_square;
I_if_square_cell{1, 2} = I_if_c_square;
I_if_square_cell{1, 3} = I_if_0c_square;
I_if_square_cell{1, 4} = I_if_s_square;

% Return E as a single scalar if consistent
DeltaE = E(1);  

end



%Now the functions in Quarkonium -> Quarkonium +2Pions transitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% p1 -> s
function [I_if_square_cell, DeltaE, M] = P1toStrans(Nin, Nfin, Min, Mfin, length)

k = length;

% Dipion mass row
nM = (Mfin - Min) / k;
M = Min:nM:Mfin;

% Initialize variables
I_if_square_cell = cell(1, 9);

%m inical=1
I_01_m1 = zeros(1, k + 1);
I_N0_m1 = zeros(1, k + 1);
%m inical=0
I_11_m0 = zeros(1, k + 1);
%fins aquí son rellevants de calcular la resta son repetides

E = zeros(1, k + 1);

% Precompute ExpValFunc

ExpValFunc = ExpValFunctions('HQ');
transition01 = Hibrid_ItoF('HQp1m1tos0I01');
transitionN0 = Hibrid_ItoF('HQp1m1tos0IN0');
transition11 = Hibrid_ItoF('HQp1m0tos0I11');

% Compute first E value and check if consistent
[I_01_m1(1), E(1)] = ExpValFunc(Nin, Nfin, M(1), transition01, 1, 0, false);
[I_N0_m1(1), ~] = ExpValFunc(Nin, Nfin, M(1), transitionN0, 1, 0, false);
[I_11_m0(1), ~] = ExpValFunc(Nin, Nfin, M(1), transition11, 1, 0, false);

% Fill the rest of the arrays
for n = 2:k+1
    [I_01_m1(n), E(n)] = ExpValFunc(Nin, Nfin, M(n), transition01, 1, 0, false);
    [I_N0_m1(n), ~] = ExpValFunc(Nin, Nfin, M(n), transitionN0, 1, 0, false);
    [I_11_m0(n), ~] = ExpValFunc(Nin, Nfin, M(n), transition11, 1, 0, false);

    % Check if E is constant
    if E(n) ~= E(1)
        error('E values are not consistent across iterations.');
    end
end

% un cop calculades les transicions rellevants calculem les repetides
    %m inical =0
I_NN_m0 = I_11_m0;
    %m inical=-1
I_0N_mN = - I_01_m1;
I_10_mN = I_N0_m1;


% Compute spin average de cada terme = com nomes hi ha 1 m possible per a la transicio
%directament per a cada
I_if_01_square = (I_01_m1.^2)./3;
I_if_11_square = (I_11_m0.^2)./3;
I_if_0N_square = (I_0N_mN.^2)./3;
I_if_NN_square = (I_NN_m0.^2)./3;
I_if_10_square = (I_10_mN.^2)./3;
I_if_N0_square = (I_N0_m1.^2)./3;
I_if_0N_10 = (I_0N_mN .* I_10_mN )./3;
I_if_01_N0 = (I_01_m1 .* I_N0_m1 )./3;
I_if_11_NN = (I_11_m0 .* I_NN_m0 )./3;

% Store in cell array
I_if_square_cell{1, 1} = I_if_01_square;
I_if_square_cell{1, 2} = I_if_11_square;
I_if_square_cell{1, 3} = I_if_0N_square;
I_if_square_cell{1, 4} = I_if_NN_square;
I_if_square_cell{1, 5} = I_if_10_square;
I_if_square_cell{1, 6} = I_if_N0_square;
I_if_square_cell{1, 7} = I_if_0N_10;
I_if_square_cell{1, 8} = I_if_01_N0;
I_if_square_cell{1, 9} = I_if_11_NN;


% Return E as a single scalar if consistent
DeltaE = E(1);

end


% p1 -> d
function [I_if_square_cell, DeltaE, M] = P1toDtrans(Nin, Nfin, Min, Mfin, length)

k = length;

% Dipion mass row
nM = (Mfin - Min) / k;
M = Min:nM:Mfin;

% Initialize variables
I_if_square_cell = cell(1, 9);

%m inical=1
I_10_m1 = zeros(1, k + 1);
I_11_m1 = zeros(1, k + 1);
I_01_m1 = zeros(1, k + 1);
I_N0_m1 = zeros(1, k + 1);
%m inical=0
I_0N_m0 = zeros(1, k + 1);
I_11_m0 = zeros(1, k + 1);

%fins aquí son rellevants de calcular la resta son repetides

E = zeros(1, k + 1);

% Precompute ExpValFunc

ExpValFunc = ExpValFunctions('HQ');
transition10 = Hibrid_ItoF('HQp1m1tod2I10');
transition11 = Hibrid_ItoF('HQp1m1tod1I11');
transition0N = Hibrid_ItoF('HQp1m0tod1I0N');
transition01 = Hibrid_ItoF('HQp1m1tod0I01');
transitionN0 = Hibrid_ItoF('HQp1m1tod0IN0');
transition110 = Hibrid_ItoF('HQp1m0tod0I11');

% Compute first E value and check if consistent
[I_10_m1(1), E(1)] = ExpValFunc(Nin, Nfin, M(1), transition10, 1, 2, false);
[I_11_m1(1), E(1)] = ExpValFunc(Nin, Nfin, M(1), transition11, 1, 2, false);
[I_0N_m0(1), E(1)] = ExpValFunc(Nin, Nfin, M(1), transition0N, 1, 2, false);
[I_01_m1(1), E(1)] = ExpValFunc(Nin, Nfin, M(1), transition01, 1, 2, false);
[I_N0_m1(1), E(1)] = ExpValFunc(Nin, Nfin, M(1), transitionN0, 1, 2, false);
[I_11_m0(1), E(1)] = ExpValFunc(Nin, Nfin, M(1), transition110, 1, 2, false);

% Fill the rest of the arrays
for n = 2:k+1
    [I_10_m1(n), E(n)] = ExpValFunc(Nin, Nfin, M(n), transition10, 1, 2, false);
    [I_11_m1(n), ~] = ExpValFunc(Nin, Nfin, M(n), transition11, 1, 2, false);
    [I_0N_m0(n), ~] = ExpValFunc(Nin, Nfin, M(n), transition0N, 1, 2, false);
    [I_01_m1(n), ~] = ExpValFunc(Nin, Nfin, M(n), transition01, 1, 2, false);
    [I_N0_m1(n), ~] = ExpValFunc(Nin, Nfin, M(n), transitionN0, 1, 2, false);
    [I_11_m0(n), ~] = ExpValFunc(Nin, Nfin, M(n), transition110, 1, 2, false);

    % Check if E is constant
    if E(n) ~= E(1)
        error('E values are not consistent across iterations.');
    end
end

% un cop calculades les transicions rellevants calculem les repetides
%m inical=0
I_NN_m0 = I_11_m0;
I_01_m0 = - I_0N_m0;
%m inicial =-1
I_10_mN = I_N0_m1;
I_0N_mN = - I_01_m1;
I_N0_mN = I_10_m1;
I_NN_mN = I_11_m1;


% Compute spin average de cada terme = per a cada IXX de diferents m inical
% fem spin average en Jini=1
I_if_01_square = (I_01_m1.^2 + I_01_m0.^2)./3;
I_if_11_square = (I_11_m1.^2 + I_11_m0.^2)./3;
I_if_0N_square = (I_0N_m0.^2 + I_0N_mN.^2)./3;
I_if_NN_square = (I_NN_m0.^2 + I_NN_mN.^2)./3;
I_if_10_square = (I_10_m1.^2 + I_10_mN.^2)./3;
I_if_N0_square = (I_N0_m1.^2 + I_N0_mN.^2)./3;
I_if_0N_10 = (I_0N_mN .* I_10_mN + I_10_m1.*0 + I_0N_m0.*0)./3;
I_if_01_N0 = (I_01_m1 .* I_N0_m1 )./3;
I_if_11_NN = (I_11_m0 .* I_NN_m0 )./3;


% Store in cell array
I_if_square_cell{1, 1} = I_if_01_square;
I_if_square_cell{1, 2} = I_if_11_square;
I_if_square_cell{1, 3} = I_if_0N_square;
I_if_square_cell{1, 4} = I_if_NN_square;
I_if_square_cell{1, 5} = I_if_10_square;
I_if_square_cell{1, 6} = I_if_N0_square;
I_if_square_cell{1, 7} = I_if_0N_10;
I_if_square_cell{1, 8} = I_if_01_N0;
I_if_square_cell{1, 9} = I_if_11_NN;


% Return E as a single scalar if consistent
DeltaE = E(1);

end


% (s/d)1 -> p
function [I_if_square_cell, DeltaE, M] = SD1toPtrans(Nin, Nfin, Min, Mfin, length)

k = length;

% Dipion mass row
nM = (Mfin - Min) / k;
M = Min:nM:Mfin;

% Initialize variables
I_if_square_cell = cell(1, 9);

%m inical=1
I_11_m1 = zeros(1, k + 1);
I_NN_m1 = zeros(1, k + 1);
I_01_m1 = zeros(1, k + 1);
I_N0_m1 = zeros(1, k + 1);

%m inical=0
I_10_m0 = zeros(1, k + 1);
I_0N_m0 = zeros(1, k + 1);
I_11_m0 = zeros(1, k + 1);


%fins aquí son rellevants de calcular la resta son repetides

E = zeros(1, k + 1);

% Precompute ExpValFunc

ExpValFunc = ExpValFunctions('HQ');
transition11 = Hibrid_ItoF('HQsd1m1tp1I11');
transitionNN = Hibrid_ItoF('HQsd1m1tp1INN');
transition01 = Hibrid_ItoF('HQsd1m1tp0I01');
transitionN0 = Hibrid_ItoF('HQsd1m1tp0IN0');
transition10 = Hibrid_ItoF('HQsd1m0tp1I10');
transition0N = Hibrid_ItoF('HQsd1m0tp1I0N');
transition110 = Hibrid_ItoF('HQsd1m0tp0I11');

% Compute first E value and check if consistent
[I_11_m1(1), E(1)] = ExpValFunc(Nin, Nfin, M(1), transition11, 1, 1, true);
[I_NN_m1(1), E(1)] = ExpValFunc(Nin, Nfin, M(1), transitionNN, 1, 1, true);
[I_01_m1(1), E(1)] = ExpValFunc(Nin, Nfin, M(1), transition01, 1, 1, true);
[I_N0_m1(1), E(1)] = ExpValFunc(Nin, Nfin, M(1), transitionN0, 1, 1, true);
[I_10_m0(1), E(1)] = ExpValFunc(Nin, Nfin, M(1), transition10, 1, 1, true);
[I_0N_m0(1), E(1)] = ExpValFunc(Nin, Nfin, M(1), transition0N, 1, 1, true);
[I_11_m0(1), E(1)] = ExpValFunc(Nin, Nfin, M(1), transition110, 1, 1, true);


% Fill the rest of the arrays
for n = 2:k+1
    [I_11_m1(n), E(n)] = ExpValFunc(Nin, Nfin, M(n), transition11, 1, 1, true);
    [I_NN_m1(n), ~] = ExpValFunc(Nin, Nfin, M(n), transitionNN, 1, 1, true);
    [I_01_m1(n), ~] = ExpValFunc(Nin, Nfin, M(n), transition01, 1, 1, true);
    [I_N0_m1(n), ~] = ExpValFunc(Nin, Nfin, M(n), transitionN0, 1, 1, true);
    [I_10_m0(n), ~] = ExpValFunc(Nin, Nfin, M(n), transition10, 1, 1, true);
    [I_0N_m0(n), ~] = ExpValFunc(Nin, Nfin, M(n), transition0N, 1, 1, true);
    [I_11_m0(n), ~] = ExpValFunc(Nin, Nfin, M(n), transition110, 1, 1, true);

    % Check if E is constant
    if E(n) ~= E(1)
        error('E values are not consistent across iterations.');
    end
end

% un cop calculades les transicions rellevants calculem les repetides
%m inical=0
I_NN_m0 = I_11_m0;
I_N0_m0 = - I_10_m0;
I_01_m0 = I_0N_m0;
%m inicial =-1
I_0N_mN = I_01_m1;
I_10_mN = - I_N0_m1;
I_NN_mN = - I_11_m1;
I_11_mN = - I_NN_m1;


% Compute spin average de cada terme = per a cada IXX de diferents m inical
% fem spin average en Jini=1
I_if_01_square = (I_01_m1.^2 + I_01_m0.^2)./3;
I_if_11_square = (I_11_m1.^2 + I_11_m0.^2 +  I_11_mN.^2)./3;
I_if_0N_square = (I_0N_mN.^2 + I_0N_m0.^2)./3;
I_if_NN_square = (I_NN_m1.^2 + I_NN_m0.^2 + I_NN_mN.^2)./3;
I_if_10_square = (I_10_m0.^2 + I_10_mN.^2)./3;
I_if_N0_square = (I_N0_m1.^2 + I_N0_m0.^2)./3;
I_if_0N_10 = (I_0N_mN .* I_10_mN + I_0N_m0 .* I_10_m0)./3;
I_if_01_N0 = (I_01_m1 .* I_N0_m1 + I_01_m0 .* I_N0_m0 )./3;
I_if_11_NN = (I_11_m0 .* I_NN_m0 + I_11_m1 .* I_NN_m1 + I_11_mN .* I_NN_mN)./3;


% Store in cell array
I_if_square_cell{1, 1} = I_if_01_square;
I_if_square_cell{1, 2} = I_if_11_square;
I_if_square_cell{1, 3} = I_if_0N_square;
I_if_square_cell{1, 4} = I_if_NN_square;
I_if_square_cell{1, 5} = I_if_10_square;
I_if_square_cell{1, 6} = I_if_N0_square;
I_if_square_cell{1, 7} = I_if_0N_10;
I_if_square_cell{1, 8} = I_if_01_N0;
I_if_square_cell{1, 9} = I_if_11_NN;


% Return E as a single scalar if consistent
DeltaE = E(1);

end
