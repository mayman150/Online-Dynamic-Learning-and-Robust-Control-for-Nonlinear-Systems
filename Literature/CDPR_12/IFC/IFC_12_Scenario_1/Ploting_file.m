function Ploting_file(X, Xd, E, Phi_hat_EE, D, U0, t)
close all;
figure(1)
% -------------------------------------------------------------------------
ax111 = subplot(2,1,1);
hold on
plot(t,E(1,:),'-r','LineWidth',2);
plot(t,E(2,:),'--b','LineWidth',2);
plot(t,E(3,:),'-.','color',[0,0.55,0],'LineWidth',2);
hold
legend('$e_{x}$','$e_{y}$','$e_{z}$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
% title('Position errors','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
ylabel('m','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
% xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')

grid on;
ax111.GridLineStyle = '--';
ax111.GridAlpha = 0.5;
ax111.GridColor = 'k';
ax111.FontSize = 12;
ax111.FontName = 'Times New Roman';
% ax111.FontAngle = 'italic';
% ax111.FontAngle = 'italic';
box on
ax112 = subplot(2,1,2);
hold on
plot(t,E(4,:),'-r','LineWidth',2);
plot(t,E(5,:),'--b','LineWidth',2);
plot(t,E(6,:),'-.','color',[0,0.55,0],'LineWidth',2);
hold on
legend('$e_{\alpha}$','$e_{\beta}$','$e_{\gamma}$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
% title('Orientation errors','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
ylabel('rad','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
% title('Position errors','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
grid on;
ax112.GridLineStyle = '--';
ax112.GridAlpha = 0.5;
ax112.GridColor = 'k';
ax112.FontSize = 12;
ax112.FontName = 'Times New Roman';
% ax112.FontAngle = 'italic';
box on
%%
% -------------------------------------------------------------------------
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% postion and orientation tracking %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
% -------------------------------------------------------------------------
ax211=subplot(3,1,1);
plot(t(:),Xd(1,:),'-r','LineWidth',2);
hold on
plot(t(:),X(1,:),'-.','color',[0,0.55,0],'LineWidth',2);
legend('$x_{d}$','$x$' ,'fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
% xlabel('Time [sec]')
% title('Desired and actual positions and orientation','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
ylabel('m','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
grid on;
ax211.GridLineStyle = '--';
ax211.GridAlpha =0.5;
ax211.GridColor = 'k';
ax211.FontSize= 12;
ax211.FontName= 'Times New Roman';
% ax211.FontAngle= 'italic';
hold off
ax212=subplot(3,1,2);
plot(t(:),Xd(2,:),'-r','LineWidth',2);
hold on
plot(t(:),X(2,:),'-.','color',[0,0.55,0],'LineWidth',2);
hold on
legend('$y_{d}$','$y$' ,'fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
% xlabel('Time [sec]')
% title('Desired and actual positions and orientation','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
ylabel('m','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
grid on;
ax212.GridLineStyle = '--';
ax212.GridAlpha =0.5;
ax212.GridColor = 'k';
ax212.FontSize= 12;
ax212.FontName= 'Times New Roman';
% ax212.FontAngle= 'italic';
hold off
box on
ax213=subplot(3,1,3);
plot(t(:),Xd(3,:),'-r','LineWidth',2);
hold on
plot(t(:),X(3,:),'-.','color',[0,0.55,0],'LineWidth',2);
legend('$z_{d}$','$z$' ,'fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
% xlabel('Time [sec]')
% title('Desired and actual positions and orientation','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
ylabel('m','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
grid on;
ax213.GridLineStyle = '--';
ax213.GridAlpha =0.5;
ax213.GridColor = 'k';
ax213.FontSize= 12;
ax213.FontName= 'Times New Roman';
% ax213.FontAngle= 'italic';
hold off
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(10)
%--------------------------------------------------------------------------
ax241 = subplot(3,1,1);
plot(t(:),Xd(4,:),'-r','LineWidth',2);
hold on
plot(t(:),X(4,:),'-.','color',[0,0.55,0],'LineWidth',2);
hold on 
legend('${\alpha}_{d}$','$\alpha$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
% xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
ylabel('rad','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
grid on;
ax241.GridLineStyle = '--';
ax241.GridAlpha =0.5;
ax241.GridColor = 'k';
ax241.FontSize= 12;
ax241.FontName= 'Times New Roman';
% ax241.FontAngle= 'italic';
hold off
box on

%--------------------------------------------------------------------------
ax25 = subplot(3,1,2);
plot(t(:),Xd(5,:),'-r','LineWidth',2);
hold on
plot(t(:),X(5,:),'-.','color',[0,0.55,0],'LineWidth',2);
hold on
legend('${\beta}_{d}$','$\beta$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
% xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
ylabel('rad','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
grid on;
ax25.GridLineStyle = '--';
ax25.GridAlpha =0.5;
ax25.GridColor = 'k';
ax25.FontSize= 12;
ax25.FontName= 'Times New Roman';
% ax25.FontAngle= 'italic';
hold off
box on
%--------------------------------------------------------------------------
ax26 = subplot(3,1,3);
plot(t(:),Xd(6,:),'-r','LineWidth',2);
hold on
plot(t(:),X(6,:),'-.','color',[0,0.55,0],'LineWidth',2);
hold on
legend('${\gamma}_{d}$','$\gamma$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
ylabel('rad','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
grid on;
ax26.GridLineStyle = '--';
ax26.GridAlpha =0.5;
ax26.GridColor = 'k';
ax26.FontSize= 12;
ax26.FontName= 'Times New Roman';
% ax26.FontAngle= 'italic';
hold off
box on

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% SATURATED CONTROL INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(11)
% -------------------------------------------------------------------------
ax41 = subplot(2,1,1);
hold on
plot(t(:),U0(1,:),'-r','LineWidth',2);
plot(t(:),U0(2,:),'-k','LineWidth',2);
plot(t(:),U0(3,:),'-b','LineWidth',2);
plot(t(:),U0(4,:),'-','Color',[0.64,0.08,0.18],'LineWidth',2);
plot(t(:),U0(5,:),'-','Color',[0.8,0.8,0.8],'LineWidth',2);
plot(t(:),U0(6,:),'-','Color',[0.47,0.67,0.19],'LineWidth',2);
plot(t(:),U0(7,:),'-','Color',[0.49,0.18,0.56],'LineWidth',2);
plot(t(:),U0(8,:),'-','Color',[1,0.41,0.16],'LineWidth',2);
hold off
legend('$\tau_1$','$\tau_2$','$\tau_3$','$\tau_4$','$\tau_5$','$\tau_6$','$\tau_7$','$\tau_8$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
ylabel('N','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
% title('Positive tensions','fontsize',10,'FontName','Times New Roman','Interpreter','latex')
grid on;
ax41.GridLineStyle = '--';
ax41.GridAlpha = 0.5;
ax41.GridColor = 'k';
ax41.FontSize = 12;
ax41.FontName = 'Times New Roman';
% ax41.FontAngle = 'italic';
box on
% figure(12)
% -------------------------------------------------------------------------
% ax51 = subplot(1,1,1);
% hold on
% plot(t(:),Tau(1,:),'-r','LineWidth',2);
% plot(t(:),Tau(2,:),'-k','LineWidth',2);
% plot(t(:),Tau(3,:),'-b','LineWidth',2);
% plot(t(:),Tau(4,:),'-','Color',[0.64,0.08,0.18],'LineWidth',2);
% plot(t(:),Tau(5,:),'-','Color',[0.27,0.25,0.25],'LineWidth',2);
% plot(t(:),Tau(6,:),'-','Color',[0.47,0.67,0.19],'LineWidth',2);
% plot(t(:),Tau(7,:),'-','Color',[0.49,0.18,0.56],'LineWidth',2);
% plot(t(:),Tau(8,:),'-','Color',[1,0.41,0.16],'LineWidth',2);
% hold off
% legend('$\tau_1$','$\tau_2$','$\tau_3$','$\tau_4$','$\tau_5$','$\tau_6$','$\tau_7$','$\tau_8$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
% ylabel('N','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
% xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
% % title('Positive tensions','fontsize',10,'FontName','Times New Roman','Interpreter','latex')
% grid on;
% ax51.GridLineStyle = '--';
% ax51.GridAlpha = 0.5;
% ax51.GridColor = 'k';
% ax51.FontSize = 12;
% ax51.FontName = 'Times New Roman';
% ax51.FontAngle = 'italic';
% box on
figure(13)

ax61 = subplot(3,1,1);
hold on
plot(t(:),D(1,:),'-r','LineWidth',2);
hold on
plot(t(:),Phi_hat_EE(1,:),'--','color',[0,0.55,0],'LineWidth',2);
ylabel('$m/s^2$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
legend('${f}_{a_{1_x}}$','$\hat{f}_{a_{1_x}}$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
grid on;
ax61.GridLineStyle = '--';
ax61.GridAlpha = 0.5;
ax61.GridColor = 'k';
ax61.FontSize = 12;
ax61.FontName = 'Times New Roman';
% ax61.FontAngle = 'italic';
box on

ax62 = subplot(3,1,2);
hold on
hold on
plot(t(:),D(2,:),'-r','LineWidth',2);
hold on
plot(t(:),Phi_hat_EE(2,:),'--','color',[0,0.55,0],'LineWidth',2);
ylabel('$m/s^2$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
legend('${f}_{a_{1_y}}$','$\hat{f}_{a_{1_y}}$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
grid on;
ax62.GridLineStyle = '--';
ax62.GridAlpha = 0.5;
ax62.GridColor = 'k';
ax62.FontSize = 12;
ax62.FontName = 'Times New Roman';
% ax62.FontAngle = 'italic';
box on
hold off


ax63 = subplot(3,1,3);
hold on
hold on
plot(t(:),D(3,:),'-r','LineWidth',2);
hold on
plot(t(:),Phi_hat_EE(3,:),'--','color',[0,0.55,0],'LineWidth',2);
ylabel('$m/s^2$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
legend('${f}_{a_{1_z}}$','$\hat{f}_{a_{1_z}}$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
grid on;
ax63.GridLineStyle = '--';
ax63.GridAlpha = 0.5;
ax63.GridColor = 'k';
ax63.FontSize = 12;
ax63.FontName = 'Times New Roman';
% ax63.FontAngle = 'italic';
box on
hold off












figure(14)

ax71 = subplot(3,1,1);
hold on
plot(t(:),D(4,:),'-r','LineWidth',2);
hold on
plot(t(:),Phi_hat_EE(4,:),'--','color',[0,0.55,0],'LineWidth',2);
ylabel('$rad/s^2$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
legend('${f}_{a_{1_\alpha}}$','$\hat{f}_{a_{1_\alpha}}$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
grid on;
ax71.GridLineStyle = '--';
ax71.GridAlpha = 0.5;
ax71.GridColor = 'k';
ax71.FontSize = 12;
ax71.FontName = 'Times New Roman';
ax71.FontAngle = 'italic';
box on

ax72 = subplot(3,1,2);
hold on
hold on
plot(t(:),D(5,:),'-r','LineWidth',2);
hold on
plot(t(:),Phi_hat_EE(5,:),'--','color',[0,0.55,0],'LineWidth',2);
ylabel('$rad/s^2$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
legend('${f}_{a_{1_\beta}}$','$\hat{f}_{a_{1_\beta}}$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
grid on;
ax72.GridLineStyle = '--';
ax72.GridAlpha = 0.5;
ax72.GridColor = 'k';
ax72.FontSize = 12;
ax72.FontName = 'Times New Roman';
ax72.FontAngle = 'italic';
box on
hold off


ax73 = subplot(3,1,3);
hold on
hold on
plot(t(:),D(6,:),'-r','LineWidth',2);
hold on
plot(t(:),Phi_hat_EE(6,:),'--','color',[0,0.55,0],'LineWidth',2);
ylabel('$rad/s^2$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
legend('${f}_{a_{1_\gamma}}$','$\hat{f}_{a_{1_\gamma}}$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
grid on;
ax73.GridLineStyle = '--';
ax73.GridAlpha = 0.5;
ax73.GridColor = 'k';
ax73.FontSize = 12;
ax73.FontName = 'Times New Roman';
ax73.FontAngle = 'italic';
box on
hold off



figure(20);
% -------------------------------------------------------------------------
ax255=subplot(1,1,1);
plot(Xd(1,:),Xd(2,:),'-r','LineWidth',2);
hold on
plot(X(1,:),X(2,:),'--b','LineWidth',2);
 legend('$X_{d}$','$X$' ,'fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
% xlabel('Time [sec]')
% title('Desired and actual positions and orientation','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
ylabel('y(m)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
xlabel('x(m)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
grid on;
ax255.GridLineStyle = '--';
ax255.GridAlpha =0.5;
ax255.GridColor = 'k';
ax255.FontSize= 12;
ax255.FontName= 'Times New Roman';
% ax255.FontAngle= 'italic';
axis(ax255,'equal')

hold off
