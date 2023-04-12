function pa=parLic(pa)
% % physical constants
% pa.ec = 1.602176634e-19; % C, elementary_charge
% pa.bc = 1.380649e-23; % J.K^-1, Boltzmann_constant
% pa.pc = 6.62607015e-34; % J.s, Planck_constant
% pa.es = 3.9*8.8541878128e-12;% F.m^-1, epsilon_sio2
% pa.fe = 9.1093837015e-31;% kg, free_electron_mass
% 
% % % dimensional constants
% % par.L0 = 1e-6;% m
% % par.mu = 1;% 1m^2.V^-1.s^-1
% % par.T0 = 3000;% K
% 
% % environment parameters
% pa.TT = 300;% K
% 
% % geometrical parameters
% pa.ch = 1e-6;% m
% pa.wl = 1;% W/L ratio
% pa.ox = 2e-9;% m
% 
% % material parameters
% pa.eg = 0.76; % eV
% % par.eg = 1.60; % eV
% pa.pn = 0.46; %eV
% pa.pp = 0.30; % eV
% pa.me = 0.15; % m0
% pa.mh = 0.14; % m0
% pa.ue = 100e-4;% m^2.V^-1.s^-1
% pa.uh = 100e-4;% m^2.V^-1.s^-1
% % SRH equivalent mobility
% pa.t0=25e-14;% s
% % par.ueS = par.tau/(par.mh*par.fe/par.ec);% m^2.V^-1.s^-1
% % par.uhS = par.tau/(par.me*par.fe/par.ec);% m^2.V^-1.s^-1
% pa.tn=pa.t0;
% pa.tp=pa.t0;


%% derived parameters
pa.tn=pa.t0;
pa.tp=pa.t0;
pa.D0=4*pi*pa.fe*pa.pc^-2;
pa.De=pa.me*pa.D0;
pa.Dh=pa.mh*pa.D0;
pa.VT=pa.bc*pa.TT*pa.ec^-1;
pa.h1=0.5*pa.eg/pa.VT;
pa.h2=pa.h1*2;
pa.N0=pa.bc*pa.TT*pa.D0; % kTD0
pa.Nc=pa.N0*pa.me;% kTDe
pa.Nv=pa.N0*pa.mh;% kTDh
pa.re=pa.ch^2*(pa.ue*pa.VT*pa.Nc)^-1;% reduced coffee
pa.rh=pa.ch^2*(pa.uh*pa.VT*pa.Nv)^-1;% reduced coffee
pa.te=pa.tn*pa.Nc^-1;% reduced tau_n
pa.th=pa.tp*pa.Nv^-1;% reduced tau_p
pa.Cox=pa.es*pa.ox^-1;
pa.Cq=pa.ec^2*pa.D0;
pa.Cqe=pa.me*pa.Cq;
pa.Cqh=pa.mh*pa.Cq;
pa.ka=pa.Cq/pa.Cox;
pa.ga=pa.Cox/pa.Cq;
pa.kae=pa.ka*pa.me;
pa.kah=pa.ka*pa.mh;

pa.sI=pa.ch^2*pa.re^-1*pa.wl*pa.ec;
pa.tI=pa.ch^2*pa.rh^-1*pa.wl*pa.ec;


% K0
pa.K0=pa.kae*pa.kah*(1+pa.kae+pa.kah).^-1*pa.Cox*(pa.ue+pa.uh)*pa.wl;
end