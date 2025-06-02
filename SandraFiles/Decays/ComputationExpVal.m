
% Programa per calcular les integrals dels factors de forma I_if com a
% vectors en funcio de la dipion mass. Tambe les guarda com a .txt per
% poder plotejar-les

%Valor r0
setr0(3.964)
%massa
load("dades.mat","m_c","m_b")
setm_q(m_b)

%Change this depending on the transition

ComputationIE = TransitionsAdded('QQStoD');

[I_if_square_cell, DeltaE, M] = ComputationIE(5, 1, 0.28, 0.73, 50);

%Obtain the I_if^2 from the cell
I_if_square = I_if_square_cell{1,1};
I_if_c_square = I_if_square_cell{1,2};
I_if_0c_square = I_if_square_cell{1,3};
I_if_s_square = I_if_square_cell{1,4};


%Aquetses variables es passen a fitxers notebook
Mtrans=M';
Itrans=I_if_square';
Itrans_c=I_if_c_square';
Itrans_0c=I_if_0c_square';
Itrans_s=I_if_s_square';

% The following code is to create .txt files where the variables are stored
baseName = 'QQbottom_5s-1d_expandit';
folderPath = '/Users/sandra/Documents/Doctorat/Projectes PhD/Transicions a 2 pions/Lower order Lagrangian/TransitionData';

fileName = fullfile(folderPath, sprintf('%s.txt', baseName));
fileName_c = fullfile(folderPath, sprintf('%s_c.txt', baseName));
fileName_0c = fullfile(folderPath, sprintf('%s_0c.txt', baseName));
fileName_s = fullfile(folderPath, sprintf('%s_s.txt', baseName));

%Make sure the path exists
if ~exist(folderPath, 'dir')
    fprintf('Error: The folder "%s" does not exist.\n', folderPath);
    return; 
end

% Open to write
file = fopen(fileName, 'w');
file_c = fopen(fileName_c, 'w');
file_0c = fopen(fileName_0c, 'w');

% Writhe the heading
fprintf(file,    '%6s %12s\n', 'M', 'I_if^2');
fprintf(file_c,  '%6s %12s\n', 'M', 'I_if_c^2');
fprintf(file_0c, '%6s %12s\n', 'M', 'I_if*I_if_c');

% Wirthe the two columns
for i = 1:length(M)
    fprintf(file,    '%.4f %.15f\n', M(i), I_if_square(i));
    fprintf(file_c,  '%.4f %.15f\n', M(i), I_if_c_square(i));
    fprintf(file_0c, '%.4f %.15f\n', M(i), I_if_0c_square(i));
end


fclose(file);
fclose(file_c);
fclose(file_0c);

%Not all transitions need a _s file so only created for thoes which need it
%(p->p and d->s)
if any(I_if_s_square)
    file_s = fopen(fileName_s, 'w');
    fprintf(file_s, '%6s %12s\n', 'M', 'I_if_s^2');
    for i = 1:length(M)
        fprintf(file_s, '%.4f %.15f\n', M(i), I_if_s_square(i));
    end
    fclose(file_s);
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

