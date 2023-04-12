% main_ode_script_matR.m [required]

% integral acc
idsmat0=idsmat;
for j=1:vgslen
    vgs=vgsswp(j);
    for i=1:vdslen
        vds=vdslis(i);
        idsmat0(i,j)=AFET(parLic(parLib1(gausf(vgs))),vgs,vds);
    end
end

err_mat=abs(idsmat_edg-idsmat0).*idsmat_edg.^-1;

figure
semilogy(vdslis,err_mat)

figure
plot(vdslis,err_mat)