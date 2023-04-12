% main_ode_script_mat0.m [required]
acc=1e-6;
acc1=1e-6;
mag_Rlen=20;
mag_Rlis=10.^linspace(-8,0,mag_Rlen);
% bias -- def1
vgsswp=(3:0.3:6);
vgslen=length(vgsswp);
vdslis=(0.1:0.1:4.9)';
vdslen=length(vdslis);

% figure
% plot(mag_Rlis,'o-')

idscub=idsvec.*mag_Rlis.*nan;
timcub=timvec.*mag_Rlis.*0;

cub_rows=size(idscub,1);
cub_cols=size(idscub,2);


parfor i=1:max(14,cub_rows)
% bias -- def2
vgsswp=(3:0.3:6);
vgslen=length(vgsswp);
vdslis=(0.1:0.1:4.9)';
vdslen=length(vdslis);
% assign
[idx_vds,idx_vgs]=vec2squ(vdslen,vgslen,i);
vds=vdslis(idx_vds);
vgs=vgsswp(idx_vgs);

for j=1:mag_Rlen
mag_Rlis=10.^linspace(-8,0,mag_Rlen);
mag_R=mag_Rlis(j);
if j==1
    sol_init=solvec(i);
else
    sol_init=sol1;
end
tic
sol1=main_ode(acc,acc1,mag_R,gausf(vgs),vgs,vds,sol_init);
timcub(i,j)=toc;
x_vgs=vgs;x_vds=vds;
solw=mean(sol_post(gausf(x_vgs),sol1,x_vgs,x_vds).ids);
idscub(i,j)=solw;
end

end

% bias -- def3
vgsswp=(3:0.3:6);
vgslen=length(vgsswp);
vdslis=(0.1:0.1:4.9)';
vdslen=length(vdslis);

idsedg=idscub(:,end);
idsmat_edg=idsmat;

for z=1:cub_rows
    [zi,zj]=vec2squ(vdslen,vgslen,z);
    idsmat_edg(zi,zj)=idscub(z,end);
end

timtri=zeros(vdslen,vgslen,mag_Rlen);
for z=1:cub_rows
    [zi,zj]=vec2squ(vdslen,vgslen,z);
    for zz=1:mag_Rlen
    timtri(zi,zj,zz)=timcub(z,zz);
    end
end

timpst=mean(sum(timtri,3),2);

figure
plot(vdslis,idsmat_edg)

figure
semilogy(vdslis,timpst)
