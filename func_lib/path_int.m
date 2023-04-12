function out=path_int(vgs0,vds0,n,pt)
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
%% vec (vn,vp) path
% path lenght vec
plv=linspace(0,1,n+2);
if abs(Vds)>eg
    alpha_ve=1-eg./Vds;
    na=floor(alpha_ve*(n+1))-1;
    nb=n-na-1;
    pla=linspace(0,alpha_ve,na+2);
    plb=linspace(alpha_ve,1,nb+2);
    plu=[pla(1:na+1) plb];
end
if Vds>0
switch pt
    case 'upper'
        pvec=[zeros(1,n+1) plv
              plv(1:n+1) ones(1,n+2)];
    case 'lower'
        pvec=[plv(1:n+1) ones(1,n+2)
              zeros(1,n+1) plv];
    case 'diaga'
        pvec=[plv
              plv];
    case 'diamd'
        if Vds<=eg
            pvec=[zeros(1,n+1) plv
              plv(1:n+1) ones(1,n+2)];
        else
        pvec=[zeros(1,nb+1) plu 
              plu ones(1,nb+1)];
        end
end
elseif Vds<0
   switch pt
    case 'upper'
        pvec=[plv(1:n+1) ones(1,n+2)
              zeros(1,n+1) plv];
    case 'lower'
        pvec=[zeros(1,n+1) plv
              plv(1:n+1) ones(1,n+2)];
    case 'diaga'
        pvec=[plv
              plv];
    case 'diamd'
        if (-Vds)<=eg
            pvec=[plu(1:n+1) ones(1,n+2)
              zeros(1,n+1) plu];
        else
        pvec=[plv ones(1,nb+1) 
              zeros(1,nb+1) plv];
        end
   end
end
len_pv=size(pvec,2);
% (xi_n,xi_p)-path
xi_pth=kron(ones(1,len_pv),xi_mat(:,2)).*pvec+...
    kron(ones(1,len_pv),xi_mat(:,1)).*(1-pvec);
% psi-path
xi_iii=serial_xi(kapa_e,kapa_h,xi_pth);
% (Kn,Kp)-path
uv_pth=[log_exp_plus(xi_iii-xi_pth(1,:))
        log_exp_plus(xi_pth(2,:)-xi_iii)];
% total current
Itot=trace(uv_pth(:,1:(len_pv-1))*(diff(xi_pth,1,2))'*diag([sI tI]));
% (vn,vp)-path
vnvp=kron(ones(1,len_pv),[Vds;Vds]).*pvec+...
    kron(ones(1,len_pv),[0;0]).*(1-pvec);

% out
out.cur=Itot;
% path used to debug
out.path=vnvp;

% xi find in serial
function xi=serial_xi(a,b,xi_np)
    xi_n=xi_np(1,:);
    xi_p=xi_np(2,:);
    temp=size(xi_np,2);
    xi=zeros(1,temp);
    xi(1)=fzero(@(x)bdeq(a,b,xi_n(1),xi_p(1),x),ap_sol(a,b,xi_n(1),xi_p(1)));
    for i=2:temp
        xi(i)=fzero(@(x)bdeq(a,b,xi_n(i-1),xi_p(i-1),x),xi(i-1));
    end
end

% (s',t')=(iota_e,iota_h)^T*r(u,v)
function y=r(u,v)
y=(u.*v-ubar.*vbar)./(muh.*(u+ubar)+mue.*(v+vbar));
end
end