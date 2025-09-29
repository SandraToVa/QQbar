%Create wf plots from a selected state in this case the 4x4 matrix of
%Jcal=1 (p box)

%%Valor r0
setr0(3.964)
setL1(0.059)
setL3(-0.230)
load("dades.mat","m_c","m_b")
setm_q(m_b)
setspin(1)
% l= interpolació (0), llagures distncies (1), bad long distances (2), Vhf
% partial A (3), Vhf2 partial B (4)
setl(0)

[ES1J1p,W,x_wf]=QuarkoniumS0J1(m_q,spin);

len=length(x_wf);
x = 0:len-1;
for k = 1:5  % Loop over the 3rd dimension (5 slices)
    figure; % Create a new figure for each plot
    hold on; % Keep multiple scatter plots in the same figure
    
    for i = 1:3  % Loop over the 3 rows
        y = squeeze(W(i, :, k));  % Extract the i-th row for the k-th slice
        scatter(x, y, 'filled');  % Scatter plot with filled markers
    end
    
    hold off;
    title(['Plot for W(:,:,', num2str(k), ')']);
    xlabel('Index (0 to 164)');
    ylabel('Values');
    legend({'Row 1', 'Row 2', 'Row 3'});
    grid on;
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