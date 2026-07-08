% Quarkonium spectrum with new potentials
% Different energies for states S, P, and D

%massa
mass = load("dades.mat","m_c","m_b");
m_q = 1.496;

% arrays for the spectrum


% Quarkonium spin 0 spectrum (orederd by J=L)
[s,ws,~] = QQbarS0J0(m_q);
[p,wp,~] = QQbarS0J1(m_q);
[d,wd,~] = QQbarS0J2(m_q);

% Quarkonium spin 1 spectrum (ordered by J \neq L)
% Without mixing this gives the same as spin 0 states because we don't have
% hyperfine splitting of quarkonium at this order.
% We can see which state is it by the shape of the wave function
[j0,wj0,~] = QQbarS1J0(m_q);
[j1,wj1,~] = QQbarS1J1(m_q);
[j2,wj2,~] = QQbarS1J2(m_q);

% Hybrids spin 0 spectrum (ordered by Jcal=J \neq L)
% We can see which state is it by the shape of the wave function
[h0,wh0,~] = GQQbarS0Jcal0(m_q);
[h1,wh1,~] = GQQbarS0Jcal1(m_q);
[h2,wh2,~] = GQQbarS0Jcal2(m_q);

% Hybrids spin 1 spectrum (ordered by Jcal \neq J \neq L)
% Without mixing this gives the different as spin 0 states because we have
% hyperfine splitting of hybrids at this order.
% We can see which state is it by the shape of the wave function
%[jh0,~,~] = GQQbarS1Jcal0(m_q);
%[jh1,~,~] = GQQbarS1Jcal1(m_q);
%[jh2,~,~] = GQQbarS1Jcal2(m_q);
