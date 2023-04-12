function sol2=main_ode(acc,acc1,mag_R,eta,vgs0,vds0,init)
% % bvp4c RelTol
% acc=1e-6;
% acc1=1e-6;
% mag_R=1000000e-6;

% test parameters
Vgs=vgs0; % V
Vds=vds0; % V

%% defined parameters
par=parLib(eta);
% physical constants
elementary_charge = par.ec; % C
Boltzmann_constant = par.bc; % J.K^-1
Planck_constant = par.pc; % J.s
epsilon_Al2O3=par.es;% F.m^-1
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
R0=mag_R.*(L0/mu0)^2*(elementary_charge/free_electron_mass)*V_T0^-1;

chl=ch_l/L0;

mue=mu_e/mu0;
muh=mu_h/mu0;

Cox=epsilon_Al2O3/t_ox;% F.m^-2
Cq=elementary_charge^2*(4*pi*free_electron_mass*Planck_constant^-2);
gma0=Cox/Cq;

kapa_e=me/gma0;
kapa_h=mh/gma0;
iota_e=R0*chl^2*Vt^-1*mue^-1*me^-1;
iota_h=R0*chl^2*Vt^-1*muh^-1*mh^-1;

ubar=sqrt(mh/me)*exp(-0.5*eg/VT);
vbar=sqrt(me/mh)*exp(-0.5*eg/VT);

%% boundary conditions, denpendent with bias
%      |   x=0  |   x=1   |
xi_mat=([0+phi_n, Vds+phi_n;...          % xi_n
         0-phi_p, Vds-phi_p]-Vgs)./VT;   % xi_p
         
xi_vec=[sol_bnd(kapa_e,kapa_h,xi_mat(1,1),xi_mat(2,1)),...
        sol_bnd(kapa_e,kapa_h,xi_mat(1,2),xi_mat(2,2))];
uv_mat=[log_exp_plus(xi_vec(1)-xi_mat(1,1)),log_exp_plus(xi_vec(2)-xi_mat(1,2));... % u
        log_exp_plus(xi_mat(2,1)-xi_vec(1)),log_exp_plus(xi_mat(2,2)-xi_vec(2))];   % v

%% bvp4c
xmesh = linspace(0,1,50);
if nargin==6
solinit2 = bvpinit(xmesh, @guess4);
else 
    solinit2=init;
end
try
options=bvpset('FJacobian',@jac,'RelTol',acc,'AbsTol',acc1,'Stats','on');
sol2 = bvp4c(@bvpfcn4, @bcfcn4, solinit2,options);
catch
options = bvpset('RelTol',acc,'AbsTol',acc1,'Stats','on');
sol2 = bvp4c(@bvpfcn4, @bcfcn4, solinit2,options);
end

%% function definition
function g = guess4(x)
g = [uq(x)
     vq(x)
     sq(x)
     tq(x)];
end
% approximate $\xi_n$
function y=xinq(x)
y=(1-x).*xi_mat(1,1)+x.*xi_mat(1,2);
end
% approximate $\xi_p$
function y=xipq(x)
y=(1-x).*xi_mat(2,1)+x.*xi_mat(2,2);
end
% approximate $\xi$
% $\xi=\alpha_n * \xi_n + \alpha_p * \xi_p$
function a=aln(x)
a=(max(0,xipq(x))>xinq(x)).*kapa_e...
  ./(1+(max(0,xipq(x))>xinq(x)).*kapa_e+(min(0,xinq(x))<xipq(x)).*kapa_h);
end
function a=alp(x)
a=(min(0,xinq(x))<xipq(x)).*kapa_h...
  ./(1+(max(0,xipq(x))>xinq(x)).*kapa_e+(min(0,xinq(x))<xipq(x)).*kapa_h);
end
function xi=xiq(x)
xi=aln(x).*xinq(x)+alp(x).*xipq(x);
end
% approximate u
function y=uq(x)
y=log_exp_plus(xiq(x)-xinq(x));
end
% approximate v
function y=vq(x)
y=log_exp_plus(xipq(x)-xiq(x));
end
% approximate d$\xi$/dx
function y=dxindx(~)
y=xi_mat(1,2)-xi_mat(1,1);
end
function y=dxipdx(~)
y=xi_mat(2,2)-xi_mat(2,1);
end
% approximate du/dx
function y=dudx(x)
y=((aln(x)-1).*dxindx(x)+(alp(x)-0).*dxipdx(x))./(1+exp(-(xiq(x)-xinq(x))));
end
% approximate dv/dx
function y=dvdx(x)
y=((1-alp(x)).*dxipdx(x)+(0-aln(x)).*dxindx(x))./(1+exp(-(xipq(x)-xiq(x))));
end
% approximate ds/dx
function y=sq(x)
y=(f_ext(uq(x))+kapa_e.*uq(x)).*dudx(x)-kapa_h.*uq(x).*dvdx(x);
end
% approximate dt/dx
function y=tq(x)
y=(f_ext(vq(x))+kapa_h.*vq(x)).*dvdx(x)-kapa_e.*vq(x).*dudx(x);
end
% (u',v')^T=1/Det(M)*M*(s,t)^T
function y=DET(u,v,kapa_u,kapa_v)
    y=f_ext(u).*f_ext(v)+kapa_u.*f_ext(v).*u+kapa_v.*f_ext(u).*v;
end
% M-matrix function
function y=mfc(u,v,s,t,kapa)
y=(f_ext(v)+kapa.*v).*s+kapa.*u.*t;
end
% (s',t')=(iota_e,iota_h)^T*r(u,v)
function y=r(u,v)
y=(u.*v-ubar.*vbar)./(muh.*(u+ubar)+mue.*(v+vbar));
end
function dydx=bvpfcn4(~,y)
dydx = [mfc(y(1),y(2),y(3),y(4),kapa_h)./DET(y(1),y(2),kapa_e,kapa_h);
        mfc(y(2),y(1),y(4),y(3),kapa_e)./DET(y(2),y(1),kapa_h,kapa_e);
        iota_e.*r(y(1),y(2));
        iota_h.*r(y(1),y(2))];
end

function y=PDET(u,v,kapa_u,kapa_v)
y=-((f_ext(v)+kapa_h.*v).*Df_ext(u)+kapa_u.*f_ext(v))./DET(u,v,kapa_u,kapa_v).^2;
end

function dfdy=jac(~,y)
dfdy=[j11(y) j12(y) j13(y) j14(y)
      j21(y) j22(y) j23(y) j24(y)
      j31(y) j32(y) 0      0
      j41(y) j42(y) 0      0    ];
end
% du'/du
function dj=j11(y)
dj=PDET(y(1),y(2),kapa_e,kapa_h).*mfc(y(1),y(2),y(3),y(4),kapa_h)...
    +kapa_h.*y(4)./DET(y(1),y(2),kapa_e,kapa_h);
end
% dv'/du
function dj=j21(y)
dj=PDET(y(1),y(2),kapa_e,kapa_h).*mfc(y(2),y(1),y(4),y(3),kapa_e)...
    +kapa_e.*y(4)./DET(y(1),y(2),kapa_e,kapa_h);
end
% du'/dv
function dj=j12(y)
dj=PDET(y(2),y(1),kapa_h,kapa_e).*mfc(y(1),y(2),y(3),y(4),kapa_h)...
    +kapa_h.*y(3)./DET(y(1),y(2),kapa_e,kapa_h);
end
% dv'/dv
function dj=j22(y)
dj=PDET(y(2),y(1),kapa_h,kapa_e).*mfc(y(2),y(1),y(4),y(3),kapa_e)...
    +kapa_e.*y(3)./DET(y(1),y(2),kapa_e,kapa_h);
end
% du'/ds
function dj=j13(y)
dj=(f_ext(y(2))+kapa_h.*y(2))./DET(y(1),y(2),kapa_e,kapa_h);
end
% dv'/dt
function dj=j24(y)
dj=(f_ext(y(1))+kapa_e.*y(1))./DET(y(1),y(2),kapa_e,kapa_h);
end
% du'/dt
function dj=j14(y)
dj=kapa_h.*y(1)./DET(y(1),y(2),kapa_e,kapa_h);
end
% dv'/ds
function dj=j23(y)
dj=kapa_e.*y(2)./DET(y(1),y(2),kapa_e,kapa_h);
end
% dR, denominator of R
function dR=DR(u,v,ubar,vbar,mu_u,mu_v)
dR=(v+vbar).*(v.*mu_u+ubar.*mu_v).*(mu_v.*(u+ubar)+mu_u.*(v+vbar)).^-2;
end
% ds'/du
function dj=j31(y)
dj=iota_e.*DR(y(1),y(2),ubar,vbar,mue,muh);
end
% dt'/du
function dj=j41(y)
dj=iota_h.*DR(y(1),y(2),ubar,vbar,mue,muh);
end
% ds'/dv
function dj=j32(y)
dj=iota_e.*DR(y(2),y(1),vbar,ubar,muh,mue);
end
% dt'/dv
function dj=j42(y)
dj=iota_h.*DR(y(2),y(1),vbar,ubar,muh,mue);
end
% boundary conditions
function res=bcfcn4(ya,yb)
res=[ya(1)-uv_mat(1,1);
     ya(2)-uv_mat(2,1);
     yb(1)-uv_mat(1,2);
     yb(2)-uv_mat(2,2);
    ];
end
end