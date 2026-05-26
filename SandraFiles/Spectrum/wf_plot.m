%Create wf plots from a selected state in this case the 4x4 matrix of
%Jcal=1 (p box)

%%Valor r0
setr0(3.964)
%setL1(0.059)
%setL3(-0.230)
load("dades.mat","m_c","m_b")
setm_q(m_b)
setspin(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the inifial and final wave functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initial and final states
[ES,Wi,x_wf_i]=QuarkoniumS0J1(m_q,spin);
setspin(0)
[ED,Wf,x_wf_f]=QuarkoniumS0J1(m_q,spin);

if length(x_wf_i)~=length(x_wf_f)
     disp('Change tolerance:');
     disp([length(x_wf_i), length(x_wf_f)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writhing the wf at the same x cordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Extract the first 5 slices into 3D arrays (Shape: 1 x Nx x 5)
wfs = Wi(1, :, 1:5);
wfd = Wf(1, :, 1:5);

% We use 'spline' because wavefunctions are smooth. 
% You could also use 'linear' or 'pchip' depending on your exact needs.

%inital state
wf1s = Wi(1,:,1);
wf2s = Wi(1,:,2);
wf3s = Wi(1,:,3);
wf4s = Wi(1,:,4);
wf5s = Wi(1,:,5);

%final state in terms of the inital state vector
wf1d = interp1(x_wf_f, Wf(1,:,1), x_wf_i, 'spline');
wf2d = interp1(x_wf_f, Wf(1,:,2), x_wf_i, 'spline');
wf3d = interp1(x_wf_f, Wf(1,:,3), x_wf_i, 'spline');
wf4d = interp1(x_wf_f, Wf(1,:,4), x_wf_i, 'spline');
wf5d = interp1(x_wf_f, Wf(1,:,5), x_wf_i, 'spline');



% 2. Calculate the norms using trapz along the 2nd dimension
norm1s = sqrt(trapz(x_wf_i, abs(wf1s).^2));
norm2s = sqrt(trapz(x_wf_i, abs(wf2s).^2));
norm3s = sqrt(trapz(x_wf_i, abs(wf3s).^2));
norm4s = sqrt(trapz(x_wf_i, abs(wf4s).^2));
norm5s = sqrt(trapz(x_wf_i, abs(wf5s).^2));

norm1d = sqrt(trapz(x_wf_i, abs(wf1d).^2));
norm2d = sqrt(trapz(x_wf_i, abs(wf2d).^2));
norm3d = sqrt(trapz(x_wf_i, abs(wf3d).^2));
norm4d = sqrt(trapz(x_wf_i, abs(wf4d).^2));
norm5d = sqrt(trapz(x_wf_i, abs(wf5d).^2));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for k = 1:5  % Loop over the 3rd dimension (5 slices)
    figure; % Create a new figure for each plot
    hold on; % Keep multiple scatter plots in the same figure
    
    y = squeeze(Wi(1, :, k));  % Extract the i-th row for the k-th slice
    %scatter(x, y, 'filled');  % Scatter plot with filled markers
    plot(x_wf_i, y, 'LineWidth', 1.5);  % Solid line
   
    
    hold off;
    title(['Plot for W(:,:,', num2str(k), ')']);
    xlabel('Index (0 to 164)');
    ylabel('Values');
    grid on;
end


Wf_new = cat(3, wf1d, wf2d, wf3d, wf4d, wf5d);
for k = 1:5  % Loop over the 3rd dimension (5 slices)
    figure; % Create a new figure for each plot
    hold on; % Keep multiple scatter plots in the same figure
    
    y = squeeze(Wf_new(:, :, k));  % Extract the i-th row for the k-th slice
    %scatter(x, y, 'filled');  % Scatter plot with filled markers
    plot(x_wf_i, y, 'LineWidth', 1.5);  % Solid line

    
    hold off;
    title(['Plot for W(:,:,', num2str(k), ')']);
    xlabel('Index (0 to 164)');
    ylabel('Values');
    grid on;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Noramlization matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%ev1s = wf1s * wf1s';
evnorm1s = trapz(x_wf_i, conj(wf1s).*wf1s) / (norm1s^2);
%ev2s = wf2s * wf2s';
evnorm2s = trapz(x_wf_i, conj(wf2s).*wf2s) / (norm2s^2);
%ev3s = wf3s * wf3s';
evnorm3s = trapz(x_wf_i, conj(wf3s).*wf3s) / (norm3s^2);
%ev4s = wf4s * wf4s';
evnorm4s = trapz(x_wf_i, conj(wf4s).*wf4s) / (norm4s^2);
%ev5s = wf5s * wf5s';
evnorm5s = trapz(x_wf_i, conj(wf5s).*wf5s) / (norm5s^2);

% 1s initial state
ev12 = trapz(x_wf_i, conj(wf2s).*wf1s)/ (norm1s * norm2s);
ev13 = trapz(x_wf_i, conj(wf3s).*wf1s)/ (norm1s * norm3s);
ev14 = trapz(x_wf_i, conj(wf4s).*wf1s)/ (norm1s * norm4s);
ev15 = trapz(x_wf_i, conj(wf5s).*wf1s)/ (norm1s * norm5s);

% 2s initial state
ev21 = trapz(x_wf_i, conj(wf1s).*wf2s)/ (norm2s * norm1s);
ev23 = trapz(x_wf_i, conj(wf3s).*wf2s)/ (norm2s * norm3s);
ev24 = trapz(x_wf_i, conj(wf4s).*wf2s)/ (norm2s * norm4s);
ev25 = trapz(x_wf_i, conj(wf5s).*wf2s)/ (norm2s * norm5s);

% 3s initial state
ev31 = trapz(x_wf_i, conj(wf1s).*wf3s)/ (norm3s * norm1s);
ev32 = trapz(x_wf_i, conj(wf2s).*wf3s)/ (norm3s * norm2s);
ev34 = trapz(x_wf_i, conj(wf4s).*wf3s)/ (norm3s * norm4s);
ev35 = trapz(x_wf_i, conj(wf5s).*wf3s)/ (norm3s * norm5s);

% 4s initial state
ev41 = trapz(x_wf_i, conj(wf1s).*wf4s)/ (norm4s * norm1s);
ev42 = trapz(x_wf_i, conj(wf2s).*wf4s)/ (norm4s * norm2s);
ev43 = trapz(x_wf_i, conj(wf3s).*wf4s)/ (norm4s * norm3s);
ev45 = trapz(x_wf_i, conj(wf5s).*wf4s)/ (norm4s * norm5s);

% 5s initial state
ev51 = trapz(x_wf_i, conj(wf1s).*wf5s)/ (norm5s * norm1s);
ev52 = trapz(x_wf_i, conj(wf2s).*wf5s)/ (norm5s * norm2s);
ev53 = trapz(x_wf_i, conj(wf3s).*wf5s)/ (norm5s * norm3s);
ev54 = trapz(x_wf_i, conj(wf4s).*wf5s)/ (norm5s * norm4s);

matrix = zeros(5,5);
matrix(1,:)=[evnorm1s,ev12,ev13,ev14,ev15];
matrix(2,:)=[ev21,evnorm2s,ev23,ev24,ev25];
matrix(3,:)=[ev31,ev32,evnorm3s,ev34,ev35];
matrix(4,:)=[ev41,ev42,ev43,evnorm4s,ev45];
matrix(5,:)=[ev51,ev52,ev53,ev54,evnorm5s];

len=length(x_wf_i);
res=0;
for i=1:len
    prod = Wi(1,i,2) .* Wi(1,i,2);
    res = res + prod;
end
res=sqrt(res);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Operator expectation value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Matriu r
O = 1./ x_wf_i;
%Ov2 = diag(eye(len,len))';

u1s = wf1s ./ x_wf_i;
u2s = wf2s ./ x_wf_i;
u3s = wf3s ./ x_wf_i;
u4s = wf4s ./ x_wf_i;
u5s = wf5s ./ x_wf_i;

% 1s initial state
r11 = trapz(x_wf_i, conj(wf1s) .* (O .* wf1s))/ (norm1s * norm1s);
r12 = trapz(x_wf_i, conj(wf2s) .* (O .* wf1s))/ (norm1s * norm2s);
r13 = trapz(x_wf_i, conj(wf3s) .* (O .* wf1s))/ (norm1s * norm3s);
r14 = trapz(x_wf_i, conj(wf4s) .* (O .* wf1s))/ (norm1s * norm4s);
r15 = trapz(x_wf_i, conj(wf5s) .* (O .* wf1s))/ (norm1s * norm5s);

% 1s initial state and d final state
r11sd = trapz(x_wf_i, conj(wf1d) .* (O .* wf1s))/ (norm1s * norm1d);
r12sd = trapz(x_wf_i, conj(wf2d) .* (O .* wf1s))/ (norm1s * norm2d);
r13sd = trapz(x_wf_i, conj(wf3d) .* (O .* wf1s))/ (norm1s * norm3d);
r14sd = trapz(x_wf_i, conj(wf4d) .* (O .* wf1s))/ (norm1s * norm4d);
r15sd = trapz(x_wf_i, conj(wf5d) .* (O .* wf1s))/ (norm1s * norm5d);

% 2s initial state
r21 = trapz(x_wf_i, conj(wf1s) .* (O .* wf2s))/ (norm2s * norm1s);
r22 = trapz(x_wf_i, conj(wf2s) .* (O .* wf2s))/ (norm2s * norm2s);
r23 = trapz(x_wf_i, conj(wf3s) .* (O .* wf2s))/ (norm2s * norm3s);
r24 = trapz(x_wf_i, conj(wf4s) .* (O .* wf2s))/ (norm2s * norm4s);
r25 = trapz(x_wf_i, conj(wf5s) .* (O .* wf2s))/ (norm2s * norm5s);

% 2s initial state and d final state
r21sd = trapz(x_wf_i, conj(wf1d) .* (O .* wf2s))/ (norm2s * norm1d);
r22sd = trapz(x_wf_i, conj(wf2d) .* (O .* wf2s))/ (norm2s * norm2d);
r23sd = trapz(x_wf_i, conj(wf3d) .* (O .* wf2s))/ (norm2s * norm3d);
r24sd = trapz(x_wf_i, conj(wf4d) .* (O .* wf2s))/ (norm2s * norm4d);
r25sd = trapz(x_wf_i, conj(wf5d) .* (O .* wf2s))/ (norm2s * norm5d);

% 3s initial state
r31 = trapz(x_wf_i, conj(wf1s) .* (O .* wf3s))/ (norm3s * norm1s);
r32 = trapz(x_wf_i, conj(wf2s) .* (O .* wf3s))/ (norm3s * norm2s);
r33 = trapz(x_wf_i, conj(wf3s) .* (O .* wf3s))/ (norm3s * norm3s);
r34 = trapz(x_wf_i, conj(wf4s) .* (O .* wf3s))/ (norm3s * norm4s);
r35 = trapz(x_wf_i, conj(wf5s) .* (O .* wf3s))/ (norm3s * norm5s);

% 3s initial state and d final state
r31sd = trapz(x_wf_i, conj(wf1d) .* (O .* wf3s))/ (norm3s * norm1d);
r32sd = trapz(x_wf_i, conj(wf2d) .* (O .* wf3s))/ (norm3s * norm2d);
r33sd = trapz(x_wf_i, conj(wf3d) .* (O .* wf3s))/ (norm3s * norm3d);
r34sd = trapz(x_wf_i, conj(wf4d) .* (O .* wf3s))/ (norm3s * norm4d);
r35sd = trapz(x_wf_i, conj(wf5d) .* (O .* wf3s))/ (norm3s * norm5d);

% 4s initial state
r41 = trapz(x_wf_i, conj(wf1s) .* (O .* wf4s))/ (norm4s * norm1s);
r42 = trapz(x_wf_i, conj(wf2s) .* (O .* wf4s))/ (norm4s * norm2s);
r43 = trapz(x_wf_i, conj(wf3s) .* (O .* wf4s))/ (norm4s * norm3s);
r44 = trapz(x_wf_i, conj(wf4s) .* (O .* wf4s))/ (norm4s * norm4s);
r45 = trapz(x_wf_i, conj(wf5s) .* (O .* wf4s))/ (norm4s * norm5s);

% 4s initial state and d final state
r41sd = trapz(x_wf_i, conj(wf1d) .* (O .* wf4s))/ (norm4s * norm1d);
r42sd = trapz(x_wf_i, conj(wf2d) .* (O .* wf4s))/ (norm4s * norm2d);
r43sd = trapz(x_wf_i, conj(wf3d) .* (O .* wf4s))/ (norm4s * norm3d);
r44sd = trapz(x_wf_i, conj(wf4d) .* (O .* wf4s))/ (norm4s * norm4d);
r45sd = trapz(x_wf_i, conj(wf5d) .* (O .* wf4s))/ (norm4s * norm5d);

% 5s initial state
r51 = trapz(x_wf_i, conj(wf1s) .* (O .* wf5s))/ (norm5s * norm1s);
r52 = trapz(x_wf_i, conj(wf2s) .* (O .* wf5s))/ (norm5s * norm2s);
r53 = trapz(x_wf_i, conj(wf3s) .* (O .* wf5s))/ (norm5s * norm3s);
r54 = trapz(x_wf_i, conj(wf4s) .* (O .* wf5s))/ (norm5s * norm4s);
r55 = trapz(x_wf_i, conj(wf5s) .* (O .* wf5s))/ (norm5s * norm5s);

% 1s initial state and d final state
r51sd = trapz(x_wf_i, conj(wf1d) .* (O .* wf5s))/ (norm5s * norm1d);
r52sd = trapz(x_wf_i, conj(wf2d) .* (O .* wf5s))/ (norm5s * norm2d);
r53sd = trapz(x_wf_i, conj(wf3d) .* (O .* wf5s))/ (norm5s * norm3d);
r54sd = trapz(x_wf_i, conj(wf4d) .* (O .* wf5s))/ (norm5s * norm4d);
r55sd = trapz(x_wf_i, conj(wf5d) .* (O .* wf5s))/ (norm5s * norm5d);

matrixR = zeros(5,5);
matrixR(1,:)=[r11,r12,r13,r14,r15];
matrixR(2,:)=[r21,r22,r23,r24,r25];
matrixR(3,:)=[r31,r32,r33,r34,r35];
matrixR(4,:)=[r41,r42,r43,r44,r45];
matrixR(5,:)=[r51,r52,r53,r54,r55];

% Matriz final s and d
matrixR_sd = zeros(5,5);
matrixR_sd(1,:) = [r11sd,r12sd,r13sd,r14sd,r15sd];
matrixR_sd(2,:) = [r21sd,r22sd,r23sd,r24sd,r25sd];
matrixR_sd(3,:) = [r31sd,r32sd,r33sd,r34sd,r35sd];
matrixR_sd(4,:) = [r41sd,r42sd,r43sd,r44sd,r45sd];
matrixR_sd(5,:) = [r51sd,r52sd,r53sd,r54sd,r55sd];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% second operator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Matriu r
O = 1./ (x_wf_i .^2);


% 1s initial state
r11 = trapz(x_wf_i, conj(wf1s) .* (O .* wf1s))/ (norm1s * norm1s);
r12 = trapz(x_wf_i, conj(wf2s) .* (O .* wf1s))/ (norm1s * norm2s);
r13 = trapz(x_wf_i, conj(wf3s) .* (O .* wf1s))/ (norm1s * norm3s);
r14 = trapz(x_wf_i, conj(wf4s) .* (O .* wf1s))/ (norm1s * norm4s);
r15 = trapz(x_wf_i, conj(wf5s) .* (O .* wf1s))/ (norm1s * norm5s);

% 1s initial state and d final state
r11sd = trapz(x_wf_i, conj(wf1d) .* (O .* wf1s))/ (norm1s * norm1d);
r12sd = trapz(x_wf_i, conj(wf2d) .* (O .* wf1s))/ (norm1s * norm2d);
r13sd = trapz(x_wf_i, conj(wf3d) .* (O .* wf1s))/ (norm1s * norm3d);
r14sd = trapz(x_wf_i, conj(wf4d) .* (O .* wf1s))/ (norm1s * norm4d);
r15sd = trapz(x_wf_i, conj(wf5d) .* (O .* wf1s))/ (norm1s * norm5d);

% 2s initial state
r21 = trapz(x_wf_i, conj(wf1s) .* (O .* wf2s))/ (norm2s * norm1s);
r22 = trapz(x_wf_i, conj(wf2s) .* (O .* wf2s))/ (norm2s * norm2s);
r23 = trapz(x_wf_i, conj(wf3s) .* (O .* wf2s))/ (norm2s * norm3s);
r24 = trapz(x_wf_i, conj(wf4s) .* (O .* wf2s))/ (norm2s * norm4s);
r25 = trapz(x_wf_i, conj(wf5s) .* (O .* wf2s))/ (norm2s * norm5s);

% 2s initial state and d final state
r21sd = trapz(x_wf_i, conj(wf1d) .* (O .* wf2s))/ (norm2s * norm1d);
r22sd = trapz(x_wf_i, conj(wf2d) .* (O .* wf2s))/ (norm2s * norm2d);
r23sd = trapz(x_wf_i, conj(wf3d) .* (O .* wf2s))/ (norm2s * norm3d);
r24sd = trapz(x_wf_i, conj(wf4d) .* (O .* wf2s))/ (norm2s * norm4d);
r25sd = trapz(x_wf_i, conj(wf5d) .* (O .* wf2s))/ (norm2s * norm5d);

% 3s initial state
r31 = trapz(x_wf_i, conj(wf1s) .* (O .* wf3s))/ (norm3s * norm1s);
r32 = trapz(x_wf_i, conj(wf2s) .* (O .* wf3s))/ (norm3s * norm2s);
r33 = trapz(x_wf_i, conj(wf3s) .* (O .* wf3s))/ (norm3s * norm3s);
r34 = trapz(x_wf_i, conj(wf4s) .* (O .* wf3s))/ (norm3s * norm4s);
r35 = trapz(x_wf_i, conj(wf5s) .* (O .* wf3s))/ (norm3s * norm5s);

% 3s initial state and d final state
r31sd = trapz(x_wf_i, conj(wf1d) .* (O .* wf3s))/ (norm3s * norm1d);
r32sd = trapz(x_wf_i, conj(wf2d) .* (O .* wf3s))/ (norm3s * norm2d);
r33sd = trapz(x_wf_i, conj(wf3d) .* (O .* wf3s))/ (norm3s * norm3d);
r34sd = trapz(x_wf_i, conj(wf4d) .* (O .* wf3s))/ (norm3s * norm4d);
r35sd = trapz(x_wf_i, conj(wf5d) .* (O .* wf3s))/ (norm3s * norm5d);

% 4s initial state
r41 = trapz(x_wf_i, conj(wf1s) .* (O .* wf4s))/ (norm4s * norm1s);
r42 = trapz(x_wf_i, conj(wf2s) .* (O .* wf4s))/ (norm4s * norm2s);
r43 = trapz(x_wf_i, conj(wf3s) .* (O .* wf4s))/ (norm4s * norm3s);
r44 = trapz(x_wf_i, conj(wf4s) .* (O .* wf4s))/ (norm4s * norm4s);
r45 = trapz(x_wf_i, conj(wf5s) .* (O .* wf4s))/ (norm4s * norm5s);

% 4s initial state and d final state
r41sd = trapz(x_wf_i, conj(wf1d) .* (O .* wf4s))/ (norm4s * norm1d);
r42sd = trapz(x_wf_i, conj(wf2d) .* (O .* wf4s))/ (norm4s * norm2d);
r43sd = trapz(x_wf_i, conj(wf3d) .* (O .* wf4s))/ (norm4s * norm3d);
r44sd = trapz(x_wf_i, conj(wf4d) .* (O .* wf4s))/ (norm4s * norm4d);
r45sd = trapz(x_wf_i, conj(wf5d) .* (O .* wf4s))/ (norm4s * norm5d);

% 5s initial state
r51 = trapz(x_wf_i, conj(wf1s) .* (O .* wf5s))/ (norm5s * norm1s);
r52 = trapz(x_wf_i, conj(wf2s) .* (O .* wf5s))/ (norm5s * norm2s);
r53 = trapz(x_wf_i, conj(wf3s) .* (O .* wf5s))/ (norm5s * norm3s);
r54 = trapz(x_wf_i, conj(wf4s) .* (O .* wf5s))/ (norm5s * norm4s);
r55 = trapz(x_wf_i, conj(wf5s) .* (O .* wf5s))/ (norm5s * norm5s);

% 1s initial state and d final state
r51sd = trapz(x_wf_i, conj(wf1d) .* (O .* wf5s))/ (norm5s * norm1d);
r52sd = trapz(x_wf_i, conj(wf2d) .* (O .* wf5s))/ (norm5s * norm2d);
r53sd = trapz(x_wf_i, conj(wf3d) .* (O .* wf5s))/ (norm5s * norm3d);
r54sd = trapz(x_wf_i, conj(wf4d) .* (O .* wf5s))/ (norm5s * norm4d);
r55sd = trapz(x_wf_i, conj(wf5d) .* (O .* wf5s))/ (norm5s * norm5d);

matrixR2 = zeros(5,5);
matrixR2(1,:)=[r11,r12,r13,r14,r15];
matrixR2(2,:)=[r21,r22,r23,r24,r25];
matrixR2(3,:)=[r31,r32,r33,r34,r35];
matrixR2(4,:)=[r41,r42,r43,r44,r45];
matrixR2(5,:)=[r51,r52,r53,r54,r55];

% Matriz final s and d
matrixR_sd2 = zeros(5,5);
matrixR_sd2(1,:) = [r11sd,r12sd,r13sd,r14sd,r15sd];
matrixR_sd2(2,:) = [r21sd,r22sd,r23sd,r24sd,r25sd];
matrixR_sd2(3,:) = [r31sd,r32sd,r33sd,r34sd,r35sd];
matrixR_sd2(4,:) = [r41sd,r42sd,r43sd,r44sd,r45sd];
matrixR_sd2(5,:) = [r51sd,r52sd,r53sd,r54sd,r55sd];


%%

Wf_new(1,:,1) = interp1(x_wf_f, Wf(1,:,1), x_wf_i, 'spline');
Wf_new(2,:,1) = interp1(x_wf_f, Wf(2,:,1), x_wf_i, 'spline');
Wf_new(3,:,1) = interp1(x_wf_f, Wf(3,:,1), x_wf_i, 'spline');
Wf_new(1,:,2) = interp1(x_wf_f, Wf(1,:,2), x_wf_i, 'spline');
Wf_new(2,:,2) = interp1(x_wf_f, Wf(2,:,2), x_wf_i, 'spline');
Wf_new(3,:,2) = interp1(x_wf_f, Wf(3,:,2), x_wf_i, 'spline');
Wf_new(1,:,3) = interp1(x_wf_f, Wf(1,:,3), x_wf_i, 'spline');
Wf_new(2,:,3) = interp1(x_wf_f, Wf(2,:,3), x_wf_i, 'spline');
Wf_new(3,:,3) = interp1(x_wf_f, Wf(3,:,3), x_wf_i, 'spline');
Wf_new(1,:,4) = interp1(x_wf_f, Wf(1,:,4), x_wf_i, 'spline');
Wf_new(2,:,4) = interp1(x_wf_f, Wf(2,:,4), x_wf_i, 'spline');
Wf_new(3,:,4) = interp1(x_wf_f, Wf(3,:,4), x_wf_i, 'spline');
Wf_new(1,:,5) = interp1(x_wf_f, Wf(1,:,5), x_wf_i, 'spline');
Wf_new(2,:,5) = interp1(x_wf_f, Wf(2,:,5), x_wf_i, 'spline');
Wf_new(3,:,5) = interp1(x_wf_f, Wf(3,:,5), x_wf_i, 'spline');

for k = 1:5  % Loop over the 3rd dimension (5 slices)
    figure; % Create a new figure for each plot
    hold on; % Keep multiple scatter plots in the same figure
    
    for i = 1:3  % Loop over the 3 rows
        y = squeeze(Wi(i, :, k));  % Extract the i-th row for the k-th slice
        %scatter(x, y, 'filled');  % Scatter plot with filled markers
        plot(x_wf_i, y, 'LineWidth', 1.5);  % Solid line
    end
    
    hold off;
    title(['Plot for W(:,:,', num2str(k), ')']);
    xlabel('Index (0 to 164)');
    ylabel('Values');
    legend({'Row 1', 'Row 2', 'Row 3'});
    grid on;
end

x = 0:len-1;
for k = 1:5  % Loop over the 3rd dimension (5 slices)
    figure; % Create a new figure for each plot
    hold on; % Keep multiple scatter plots in the same figure
    
    for i = 1:3  % Loop over the 3 rows
        y = squeeze(Wf_new(i, :, k));  % Extract the i-th row for the k-th slice
        %scatter(x, y, 'filled');  % Scatter plot with filled markers
        plot(x_wf_i, y, 'LineWidth', 1.5);  % Solid line
    end
    
    hold off;
    title(['Plot for W(:,:,', num2str(k), ')']);
    xlabel('Index (0 to 164)');
    ylabel('Values');
    legend({'Row 1', 'Row 2', 'Row 3'});
    grid on;
end

%-----------------------------------------------------------------------

%%

% x_wf      : vector de posiciones (x > 0)
% psi    : función de onda evaluada en x

psi=wf1s./norm1s;

% 1. Seleccionar región de x pequeños
x_max_small = 0.5;                % ajusta según tu problema
mask = x_wf_i < x_max_small;

% 2. Ajuste lineal psi(x) = a*x + b en la región pequeña
p = polyfit(x_wf_i(mask), psi(mask), 1);
a = p(1);
b = p(2);

fprintf('Pendiente a = %.6e\n', a);
fprintf('Intercepto b = %.6e\n', b);

% 3. Comprobación cuantitativa de linealidad
criterio = abs(b) / abs(a * max(x_wf_i(mask)));
fprintf('b / (a*x_max) = %.6e\n', criterio);

if criterio < 1e-2
    disp('La funcion de onda es lineal cerca del origen')
else
    disp('La funcion de onda NO es lineal cerca del origen')
end

% 4. Gráfica: datos + ajuste lineal
figure
plot(x_wf_i, psi, 'bo-', 'DisplayName', '\psi(x)'); hold on
plot(x_wf_i(mask), a*x_wf_i(mask) + b, 'r--', 'LineWidth', 2, ...
     'DisplayName', 'ajuste lineal')
xlabel('x')
ylabel('\psi(x)')
legend
grid on

% 5. Gráfica del cociente psi(x)/x (solo x pequeños)
figure
plot(x_wf_i(mask), psi(mask) ./ x_wf_i(mask), 'ko-')
xlabel('x')
ylabel('\psi(x) / x')
title('Comprobacion de linealidad cerca del origen')
grid on


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