acc=1e-6;
acc1=1e-6;
mag_R=1e-6;

%% sol0
vgs=3;
vds=0.1;

% raw solution structure
sol0=main_ode(acc,acc1,mag_R,gausf(vgs),vgs,vds);

% ids value post-process
solu=mean(sol_post(gausf(vgs),sol0,vgs,vds).ids);

%% sol1
vgs=3;
vds=0.3;

% raw solution structure
sol1=main_ode(acc,acc1,mag_R,gausf(vgs),vgs,vds);

% ids value post-process
solv=mean(sol_post(gausf(vgs),sol1,vgs,vds).ids);


%% solx, (a,b) --> x
a_vgs=3;
a_vds=0.1;
a_sol=sol0;

b_vgs=3;
b_vds=0.3;
b_sol=sol1;

x_vgs=3;
x_vds=0.5;
x_sol=interp_sol(a_vgs,a_vds,a_sol, ...
                 b_vgs,b_vds,b_sol, ...
                 x_vgs,x_vds);

%% faster

tic
% raw solution structure
solx=main_ode(acc,acc1,mag_R,gausf(x_vgs),x_vgs,x_vds,x_sol);

% ids value post-process
solw=mean(sol_post(gausf(x_vgs),solx,x_vgs,x_vds).ids);
toc

%% slower

tic
% raw solution structure
solx=main_ode(acc,acc1,mag_R,gausf(x_vgs),x_vgs,x_vds);

% ids value post-process
solw=mean(sol_post(gausf(x_vgs),solx,x_vgs,x_vds).ids);
toc

