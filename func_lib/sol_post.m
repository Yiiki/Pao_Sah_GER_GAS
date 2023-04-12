function out=sol_post(eta,sol2,vgs0,vds0)
% test parameters
Vgs=vgs0; % V
Vds=vds0; % V

%% defined parameters
par=parLib(eta);
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
         
xi_vec=[sol_bnd(kapa_e,kapa_h,xi_mat(1,1),xi_mat(2,1)),...
        sol_bnd(kapa_e,kapa_h,xi_mat(1,2),xi_mat(2,2))];
uv_mat=[log_exp_plus(xi_vec(1)-xi_mat(1,1)),log_exp_plus(xi_vec(2)-xi_mat(1,2));... % u
        log_exp_plus(xi_mat(2,1)-xi_vec(1)),log_exp_plus(xi_mat(2,2)-xi_vec(2))];   % v


% (psi, vn, vp, jn, jp)-x
xi=-(kapa_e.*sol2.y(1,:)-kapa_h.*sol2.y(2,:));
xin=xi-log_exp_minus(sol2.y(1,:));
xip=xi+log_exp_minus(sol2.y(2,:));

psi=VT.*xi+Vgs;
%  vn=-Vt.*log(exp(g0./me.*sol2.y(1,:))-1)-psb+psix;
%  vp=+Vt.*log(exp(g0./mh.*sol2.y(3,:))-1)+psc+psix;
vn=VT.*xin+Vgs-phi_n;
vp=VT.*xip+Vgs+phi_p;
% W/L=1, unit A
In=-sI.*sol2.y(3,:);
Ip=tI.*sol2.y(4,:);

% recombination rate
rr=4.*pi.*elementary_charge.^2.*Planck_constant.^-2.*VT.*...
    r(sol2.y(1,:),sol2.y(2,:));
% out
out.mesh=sol2.x;
out.u=sol2.y(1,:);
out.v=sol2.y(2,:);
out.bd=[-vn;-vp;phi_n-psi;-phi_p-psi];
out.vn=vn;
out.vp=vp;
out.in=In;
out.ip=Ip;
out.ids=In+Ip;
out.gs=Vgs;
out.ds=Vds;
out.rr=rr;
% output
% subplot(2,2,1)
% semilogy(sol2.x, sol2.y([1 2],:),'LineWidth',1.75)
% xlabel('$x/L$','FontSize',12,'Interpreter','latex')
% ylabel('$n/(kTD_{e,h})$','FontSize',12,'Interpreter','latex')
% subplot(2,2,2)
% plot(sol2.x, out.bd,'LineWidth',1.75)
% xlabel('$x/L$','FontSize',12,'Interpreter','latex')
% ylabel('$E(eV)$','FontSize',12,'Interpreter','latex')
% subplot(2,2,3)
% plot(vn,vp,'LineWidth',1.75)
% xlim([min(0,Vds),max(0,Vds)])
% ylim([min(0,Vds),max(0,Vds)])
% xlabel('$V_n(V)$','FontSize',12,'Interpreter','latex')
% ylabel('$V_p(V)$','FontSize',12,'Interpreter','latex')
% subplot(2,2,4)
% xlabel('$x/L$','FontSize',12,'Interpreter','latex')
% yyaxis left
% plot(sol2.x,[In;In+Ip],'LineWidth',1.75)
% ylabel('$I_e(A)$','FontSize',12,'Interpreter','latex')
% yyaxis right
% plot(sol2.x,Ip,'LineWidth',1.75)
% ylabel('$I_h(A)$','FontSize',12,'Interpreter','latex')

% (s',t')=(iota_e,iota_h)^T*r(u,v)
function y=r(u,v)
y=(u.*v-ubar.*vbar)./(muh.*(u+ubar)+mue.*(v+vbar));
end
end