function vnvp_plot(sol_pos)
par=parLib();
eg=par.eg; % eV
Vds=sol_pos.ds;
% Vgs=sol_pos.gs;
% vn=sol_pos.vn;
% vp=sol_pos.vp;
rr=sol_pos.rr;
% cur=mean(sol_pos.in+sol_pos.ip);
% output mode 2
subplot(2,2,1)
semilogy(sol_pos.mesh, [sol_pos.u;sol_pos.v],'LineWidth',1.75)
xlabel('$x/L$','FontSize',12,'Interpreter','latex')
ylabel('$n/(kTD_{e,h})$','FontSize',12,'Interpreter','latex')
subplot(2,2,2)
plot(sol_pos.mesh, sol_pos.bd,'LineWidth',1.75)
xlabel('$x/L$','FontSize',12,'Interpreter','latex')
ylabel('$E(eV)$','FontSize',12,'Interpreter','latex')
subplot(2,2,3)
% plot(sol_pos.vn,sol_pos.vp,sol_pos.vn,sol_pos.vn+0.5.*Vds,sol_pos.vn,sol_pos.vn+eg,'LineWidth',1.75)
plot(sol_pos.vn,sol_pos.vp-sol_pos.vn,sol_pos.vn,0.5.*Vds,sol_pos.vn,eg,'LineWidth',1.75)
xlim([min(0,Vds),max(0,Vds)])
ylim([min(0,Vds),max(0,Vds)])
xlabel('$V_n(V)$','FontSize',12,'Interpreter','latex')
ylabel('$V_p(V)$','FontSize',12,'Interpreter','latex')
subplot(2,2,4)
semilogy(sol_pos.mesh,rr,'LineWidth',1.75)
xlabel('$x/L$','FontSize',12,'Interpreter','latex')
ylabel('$R(m^{-2}\cdot s^{-1})$','FontSize',12,'Interpreter','latex')
end