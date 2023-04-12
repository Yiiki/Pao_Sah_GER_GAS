function sol_plot(sol_pos,daset)
Vds=sol_pos.ds;
Vgs=sol_pos.gs;
cur=mean(sol_pos.in+sol_pos.ip);
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
plot(sol_pos.vn,sol_pos.vp,'LineWidth',1.75)
xlim([min(0,Vds),max(0,Vds)])
ylim([min(0,Vds),max(0,Vds)])
xlabel('$V_n(V)$','FontSize',12,'Interpreter','latex')
ylabel('$V_p(V)$','FontSize',12,'Interpreter','latex')
subplot(2,2,4)
semilogy(daset(1,:),daset(2,:),'LineWidth',1.75)
xlabel('$V_{gs}$(V)','FontSize',12,'Interpreter','latex')
ylabel('$I_{ds}$(A)','FontSize',12,'Interpreter','latex')
hold on
semilogy(Vgs,cur,'o')
hold off
end