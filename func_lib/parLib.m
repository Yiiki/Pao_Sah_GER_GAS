function par=parLib(eta)
% physical constants
par.ec = 1.602176634e-19; % C, elementary_charge
par.bc = 1.380649e-23; % J.K^-1, Boltzmann_constant
par.pc = 6.62607015e-34; % J.s, Planck_constant
par.es = 9.1*8.8541878128e-12;% F.m^-1, epsilon_Al2O3
par.fe = 9.1093837015e-31;% kg, free_electron_mass

% dimensional constants
par.L0 = 1e-6;% m
par.mu = 1;% 1m^2.V^-1.s^-1
par.T0 = 3000;% K

% environment parameters
par.TT = 300.*eta;% K

% geometrical parameters
par.ch = 8.4e-6;% m
par.wl = (2.0./5.4 + 2.1./8.5 + 1.9./9.4 + 1.4./11.2).^-1;% W/L ratio
par.ox = 20e-9;% m

% material parameters
par.eg = 2.1463; % eV
par.pn = 2.4155e+00; %eV
par.pp = -2.6917e-01; % eV
par.me = 0.15; % m0
par.mh = 0.14; % m0
par.ue = 3.6366e-04;% m^2.V^-1.s^-1
par.uh = 1.6449e-03;% m^2.V^-1.s^-1
end