function pa=parLib1(eta)
% bvp4c acc
pa.acc_a=1e-6;
pa.acc_r=1e-6;
pa.nmax=1e3;

% physical constants
pa.ec = 1.602176634e-19; % C, elementary_charge
pa.bc = 1.380649e-23; % J.K^-1, Boltzmann_constant
pa.pc = 6.62607015e-34; % J.s, Planck_constant
pa.es = 9.1*8.8541878128e-12;% F.m^-1, epsilon_Al2O3
pa.fe = 9.1093837015e-31;% kg, free_electron_mass

% environment parameters
pa.TT = eta.*300;% K

% geometrical parameters
pa.ch = 8.4e-6;% m
pa.wl = (2.0./5.4 + 2.1./8.5 + 1.9./9.4 + 1.4./11.2).^-1;% W/L ratio
pa.ox = 20e-9;% m

% material parameters
pa.eg = 2.1463; % eV
pa.pn = 2.4155e+00; %eV
pa.pp = -2.6917e-01; % eV
pa.me = 0.15; % m0
pa.mh = 0.14; % m0
pa.ue = 3.6366e-04;% m^2.V^-1.s^-1
pa.uh = 1.6449e-03;% m^2.V^-1.s^-1

% SRH equivalent mobility
pa.t0 = 25e-6;% s
% pa.t0 = 25e-15;% s

end