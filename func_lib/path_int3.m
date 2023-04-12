function out=path_int3(vgs0,vds0,acc_a,acc_r)
% test parameters
Vgs=vgs0; % V
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
%% boundary conditions, denpendent with bias
%      |   x=0  |   x=1   |
xi_mat=([0+phi_n, Vds+phi_n;...          % xi_n
         0-phi_p, Vds-phi_p]-Vgs)./VT;   % xi_p
         
xi_vec=[sol_bnd(kapa_e,kapa_h,xi_mat(1,1),xi_mat(2,1)),...% x=0
        sol_bnd(kapa_e,kapa_h,xi_mat(1,2),xi_mat(2,2))];  % x=1
uv_mat=[log_exp_plus(xi_vec(1)-xi_mat(1,1)),log_exp_plus(xi_vec(2)-xi_mat(1,2));... % u
        log_exp_plus(xi_mat(2,1)-xi_vec(1)),log_exp_plus(xi_mat(2,2)-xi_vec(2))];   % v
f=@(x)udiag(x);
g=@(x)vdiag(x);

i_n=sI.*integral(f,xi_mat(1,1),xi_mat(1,2),'RelTol',acc_r,'AbsTol',acc_a);
i_p=tI.*integral(g,xi_mat(2,1),xi_mat(2,2),'RelTol',acc_r,'AbsTol',acc_a);
out=i_n+i_p;
% xi find in serial
function u=udiag(xi_n0)
    u=xi_n0;
    parfor i=1:length(xi_n0)
        xi_n=xi_n0(i);
        xi=fzero(@(x)bdeq(kapa_e,kapa_h,xi_n,xi_n-eg/VT,x),ap_sol(kapa_e,kapa_h,xi_n,xi_n-eg/VT));
        u(i)=log_exp_plus(xi-xi_n);
    end
end
function v=vdiag(xi_p0)
    v=xi_p0;
    parfor i=1:length(xi_p0)
        xi_p=xi_p0(i);
        xi=fzero(@(x)bdeq(kapa_e,kapa_h,xi_p+eg/VT,xi_p,x),ap_sol(kapa_e,kapa_h,xi_p+eg/VT,xi_p));
        v(i)=log_exp_plus(xi_p-xi);
    end
end

end