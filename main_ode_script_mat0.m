acc=1e-6;
acc1=1e-6;
mag_R=1e-8;% test a highest value that all (vgs,vds)-points can stand with

% Set a couple of warnings to temporarily issue errors (exceptions)
s = warning('error','MATLAB:bvp4c:RelTolNotMet');
steta=1;

vgsswp=(3:0.3:6);
vgslen=length(vgsswp);
vdslis=(0.1:0.1:4.9)';
vdslen=length(vdslis);

idsmat=vdslis.*vgsswp.*nan;
timmat=idsmat.^0-1;
% initial shot
sol00=main_ode(acc,acc1,mag_R,gausf(0.1),0.1,0.1);
solmat(vdslen,vgslen)=sol00;

figure

for j=1:vgslen

Vgs=vgsswp(j);
idslis=nan.*vdslis;
timlis=vdslis;

%% sol0
vgs=Vgs;
vds=vdslis(1);

tic
% raw solution structure
sol0=main_ode(acc,acc1,mag_R,gausf(vgs),vgs,vds);
tim1=toc;

% ids value post-process
solu=mean(sol_post(gausf(vgs),sol0,vgs,vds).ids);

%% sol1
vgs=Vgs;
vds=vdslis(2);

tic
% raw solution structure
sol1=main_ode(acc,acc1,mag_R,gausf(vgs),vgs,vds);
tim2=toc;

% ids value post-process
solv=mean(sol_post(gausf(vgs),sol1,vgs,vds).ids);

% filling blank
idslis(1:2)=[solu,solv]';
timlis(1:2)=[tim1,tim2]';
sollis(length(vdslis))=sol0;
sollis(1)=sol0;
sollis(2)=sol1;

for i=3:length(vdslis)
%% solx, (a,b) --> x
a_vgs=Vgs;
a_vds=vdslis(i-2);
a_sol=sollis(i-2);

b_vgs=Vgs;
b_vds=vdslis(i-1);
b_sol=sollis(i-1);

x_vgs=Vgs;
x_vds=vdslis(i);
x_sol=interp_sol(a_vgs,a_vds,a_sol, ...
                 b_vgs,b_vds,b_sol, ...
                 x_vgs,x_vds);

try
   % Regular processing part

%% faster
solx=x_sol;
tic
% raw solution structure
for z=1:1
solx=main_ode(acc,acc1,mag_R,gausf(x_vgs),x_vgs,x_vds,solx);
end
timlis(i)=toc;

catch ME
        % Exception-handling part
        switch ME.identifier
            case 'MATLAB:bvp4c:RelTolNotMet'
                fprintf('warning, RelTol Can''t Met at initial solution\n');
                steta=0;
            case 'MATLAB:bvp4c:SingJac'
                fprintf('error, Singular Jac at initial solution\n');
                steta=0;
            otherwise
                rethrow(ME)
        end
        % exception-handling
        if steta==0
            break
        end
end

% ids value post-process
solw=mean(sol_post(gausf(x_vgs),solx,x_vgs,x_vds).ids);
idslis(i)=solw;
sollis(i)=solx;

end

idx=length(idslis)-sum(isnan(idslis));

subplot(1,2,1)
plot(vdslis(1:idx),idslis(1:idx))
hold on
subplot(1,2,2)
plot(vdslis(1:idx),timlis(1:idx))
hold on

idsmat(1:idx,j)=idslis(1:idx);
timmat(1:idx,j)=timlis(1:idx);
solmat(1:idx,j)=sollis(1:idx);

end

hold off

idsvec=reshape(idsmat,[],1);
timvec=reshape(timmat,[],1);
solvec=reshape(solmat,[],1);