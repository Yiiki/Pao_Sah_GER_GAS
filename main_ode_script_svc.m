acc=1e-6;
acc1=1e-6;
mag_R=1e-6;

%
Vds=5;
vgslis=(3:0.3:6)';
idsljs=vgslis.*nan;
timljs=vgslis;

%% sol0
vds=Vds;
vgs=vgslis(1);

tic
% raw solution structure
sol0=main_ode(acc,acc1,mag_R,gausf(vgs),vgs,vds);
tim1=toc;

% ids value post-process
solu=mean(sol_post(gausf(vgs),sol0,vgs,vds).ids);

%% sol1
vds=Vds;
vgs=vgslis(2);

tic
% raw solution structure
sol1=main_ode(acc,acc1,mag_R,gausf(vgs),vgs,vds);
tim2=toc;

% ids value post-process
solv=mean(sol_post(gausf(vgs),sol1,vgs,vds).ids);

% filling blank
idsljs(1:2)=[solu,solv]';
timljs(1:2)=[tim1,tim2]';
solljs(length(vgslis))=sol0;
solljs(1)=sol0;
solljs(2)=sol1;

for i=3:length(vgslis)
%% solx, (a,b) --> x
a_vds=Vds;
a_vgs=vgslis(i-2);
a_sol=solljs(i-2);

b_vds=Vds;
b_vgs=vgslis(i-1);
b_sol=solljs(i-1);

x_vds=Vds;
x_vgs=vgslis(i);
x_sol=interp_sol(a_vgs,a_vds,a_sol, ...
                 b_vgs,b_vds,b_sol, ...
                 x_vgs,x_vds);

%% faster
solx=x_sol;
tic
% raw solution structure
for z=1:1
solx=main_ode(acc,acc1,mag_R,gausf(x_vgs),x_vgs,x_vds,solx);
end
timljs(i)=toc;

% ids value post-process
solw=mean(sol_post(gausf(x_vgs),solx,x_vgs,x_vds).ids);
idsljs(i)=solw;
solljs(i)=solx;

end

idx=length(idsljs)-sum(isnan(idsljs));
figure
plot(vgslis(1:idx),idsljs(1:idx))
figure
plot(vgslis(1:idx),timljs(1:idx))
