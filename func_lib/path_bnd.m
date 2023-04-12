function out=path_bnd(vds0,tag)
% test parameters
% Vgs=vgs0; % V
Vds=vds0; % V

%% defined parameters
par=parLib();
% physical constants
elementary_charge = par.ec; % C
Boltzmann_constant = par.bc; % J.K^-1
Planck_constant = par.pc; % J.s
epsilon_sio2=par.es;% F.m^-1
free_electron_mass =par.fe;% kg

% dimensional constants
L0=par.L0;% m
mu0=par.mu;% 1m^2.V^-1.s^-1
T0=par.T0;% K

% environment parameters
T_0=par.TT;% K

% geometrical parameters
ch_l=par.ch;% m
wlr=par.wl;% W/L ratio
t_ox=par.ox;% m

% material parameters
eg=par.eg; % eV
phi_n=par.pn; %eV
phi_p=par.pp; % eV
me=par.me; % m0
mh=par.mh; % m0
mu_e=par.ue;% m^2.V^-1.s^-1
mu_h=par.uh;% m^2.V^-1.s^-1

%% derived parameters, independent with bias
V_T0=Boltzmann_constant*T0/elementary_charge;
VT=Boltzmann_constant*T_0/elementary_charge;
Vt=VT/V_T0;
R0=(L0/mu0)^2*(elementary_charge/free_electron_mass)*V_T0^-1;

chl=ch_l/L0;

mue=mu_e/mu0;
muh=mu_h/mu0;

Cox=epsilon_sio2/t_ox;% F.m^-2
Cq=elementary_charge^2*(4*pi*free_electron_mass*Planck_constant^-2);
gma0=Cox/Cq;

kapa_e=me/gma0;
kapa_h=mh/gma0;
iota_e=R0*chl^2*Vt^-1*mue^-1*me^-1;
iota_h=R0*chl^2*Vt^-1*muh^-1*mh^-1;

ubar=sqrt(mh/me)*exp(-0.5*eg/VT);
vbar=sqrt(me/mh)*exp(-0.5*eg/VT);

% s to In conversion factor
sI=wlr.*Cox.*mu_e.*kapa_e.*VT.^2;
tI=wlr.*Cox.*mu_h.*kapa_h.*VT.^2;

% K0
K0=wlr.*Cox.*(mu_e+mu_h);
switch tag
    case 'con'
        out=K0.*kapa_e.*kapa_h.*VT.^2.*SS((Vds-eg)./VT);
    case 'opt'
        out=K0.*kapa_e.*kapa_h.*VT.^2.*SS((Vds-eg)./VT).*(1+kapa_e+kapa_h).^-1;
end
    function S=SS(x)
        S=(x+1).*(x>=0)+exp(x).*(x<0);
    end
end