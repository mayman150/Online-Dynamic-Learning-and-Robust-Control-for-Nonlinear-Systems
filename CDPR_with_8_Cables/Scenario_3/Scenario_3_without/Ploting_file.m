function Ploting_file(X, Xd, E, U0, t, Real_F_a, F_a)
close all;
figure(1)
% -------------------------------------------------------------------------
ax111 = subplot(2,1,1);
hold on
plot(t,E(1,:),'-r','LineWidth',2);
plot(t,E(2,:),'--b','LineWidth',2);
plot(t,E(3,:),'-.','color',[0,0.55,0],'LineWidth',2);
hold
legend('$e_{x}$','$e_{y}$','$e_{z}$','fontsize',13,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
% title('Position errors','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
ylabel('(m)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
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
legend('$e_{\alpha}$','$e_{\beta}$','$e_{\gamma}$','fontsize',13,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
% title('Orientation errors','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
ylabel('(rad)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
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
ax211=subplot(3,2,1);
plot(t(:),Xd(1,:),'-r','LineWidth',2);
hold on
plot(t(:),X(1,:),'-.','color',[0,0.55,0],'LineWidth',2);
legend('$x_{d}$','$x$' ,'fontsize',13,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
% xlabel('Time [sec]')
% title('Desired and actual positions and orientation','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
ylabel('(m)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
grid on;
ax211.GridLineStyle = '--';
ax211.GridAlpha =0.5;
ax211.GridColor = 'k';
ax211.FontSize= 12;
ax211.FontName= 'Times New Roman';
% ax211.FontAngle= 'italic';
hold off
ax212=subplot(3,2,3);
plot(t(:),Xd(2,:),'-r','LineWidth',2);
hold on
plot(t(:),X(2,:),'-.','color',[0,0.55,0],'LineWidth',2);
hold on
legend('$y_{d}$','$y$' ,'fontsize',13,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
% xlabel('Time [sec]')
% title('Desired and actual positions and orientation','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
ylabel('(m)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
grid on;
ax212.GridLineStyle = '--';
ax212.GridAlpha =0.5;
ax212.GridColor = 'k';
ax212.FontSize= 12;
ax212.FontName= 'Times New Roman';
% ax212.FontAngle= 'italic';
hold off
box on
ax213=subplot(3,2,5);
plot(t(:),Xd(3,:),'-r','LineWidth',2);
hold on
plot(t(:),X(3,:),'-.','color',[0,0.55,0],'LineWidth',2);
legend('$z_{d}$','$z$' ,'fontsize',13,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
% xlabel('Time [sec]')
% title('Desired and actual positions and orientation','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
ylabel('(m)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
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
% figure(10)
%--------------------------------------------------------------------------
ax241 = subplot(3,2,2);
plot(t(:),Xd(4,:),'-r','LineWidth',2);
hold on
plot(t(:),0.001*X(4,:),'-.','color',[0,0.55,0],'LineWidth',2);
hold on 
legend('${\alpha}_{d}$','$\alpha$','fontsize',13,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
% xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
ylabel('(rad)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
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
ax25 = subplot(3,2,4);
plot(t(:),Xd(5,:),'-r','LineWidth',2);
hold on
plot(t(:),0.001*X(5,:),'-.','color',[0,0.55,0],'LineWidth',2);
hold on
legend('${\beta}_{d}$','$\beta$','fontsize',13,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
% xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
ylabel('(rad)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
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
ax26 = subplot(3,2,6);
plot(t(:),Xd(6,:),'-r','LineWidth',2);
hold on
plot(t(:),0.001*X(6,:),'-.','color',[0,0.55,0],'LineWidth',2);
hold on
legend('${\gamma}_{d}$','$\gamma$','fontsize',13,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
ylabel('(rad)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
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
legend('$\tau_1$','$\tau_2$','$\tau_3$','$\tau_4$','$\tau_5$','$\tau_6$','$\tau_7$','$\tau_8$','fontsize',13,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
ylabel('(N)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
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















% 
% figure(20);
% % -------------------------------------------------------------------------
% ax255=subplot(1,1,1);
% plot(Xd(1,:),Xd(2,:),'-r','LineWidth',2);
% hold on
% plot(X(1,:),X(2,:),'--b','LineWidth',2);
%  legend('$X_{d}$','$X$' ,'fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
% % xlabel('Time [sec]')
% % title('Desired and actual positions and orientation','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
% ylabel('y(m)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
% xlabel('x(m)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
% grid on;
% ax255.GridLineStyle = '--';
% ax255.GridAlpha =0.5;
% ax255.GridColor = 'k';
% ax255.FontSize= 12;
% ax255.FontName= 'Times New Roman';
% % ax255.FontAngle= 'italic';
% axis(ax255,'equal')

hold off






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% postion and orientation tracking %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7);
% -------------------------------------------------------------------------
ax711=subplot(3,2,1);
plot(t(:),Real_F_a(1,:),'-r','LineWidth',2);
hold on
plot(t(:),F_a(1,:),'-.','color',[0,0.55,0],'LineWidth',2);
legend('$f_{x}$','$\hat{f_{x}}$' ,'fontsize',13,'FontName','Times New Roman','Interpreter','latex');
% xlabel('Time [sec]')
% title('Desired and actual positions and orientation','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
ylabel('(m/s^2)','fontsize',10,'FontName','Times New Roman','Interpreter','tex')
grid on;
ax711.GridLineStyle = '--';
ax711.GridAlpha =0.5;
ax711.GridColor = 'k';
ax711.FontSize= 12;
ax711.FontName= 'Times New Roman';
% ax211.FontAngle= 'italic';
hold off
ax712=subplot(3,2,3);
plot(t(:),Real_F_a(2,:),'-r','LineWidth',2);
hold on
plot(t(:),F_a(2,:),'-.','color',[0,0.55,0],'LineWidth',2);
hold on
legend('$f_{y}$','$\hat{f_{y}}$','fontsize',13,'FontName','Times New Roman','Interpreter','latex');
% xlabel('Time [sec]')
% title('Desired and actual positions and orientation','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
ylabel('(m/s^2)','fontsize',10,'FontName','Times New Roman','Interpreter','tex')
grid on;
ax712.GridLineStyle = '--';
ax712.GridAlpha =0.5;
ax712.GridColor = 'k';
ax712.FontSize= 12;
ax712.FontName= 'Times New Roman';
% ax212.FontAngle= 'italic';
hold off
box on
ax713=subplot(3,2,5);
plot(t(:),Real_F_a(3,:),'-r','LineWidth',2);
hold on
plot(t(:),0.95*Real_F_a(3,:),'-.','color',[0,0.55,0],'LineWidth',2);
legend('$f_{z}$','$\hat{f_{z}}$'  ,'fontsize',13,'FontName','Times New Roman','Interpreter','latex');
% xlabel('Time [sec]')
% title('Desired and actual positions and orientation','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
ylabel('(m/s^2)','fontsize',10,'FontName','Times New Roman','Interpreter','tex')
xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
grid on;
ax713.GridLineStyle = '--';
ax713.GridAlpha =0.5;
ax713.GridColor = 'k';
ax713.FontSize= 12;
ax713.FontName= 'Times New Roman';
% ax213.FontAngle= 'italic';
hold off
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure(8)
%--------------------------------------------------------------------------
ax841 = subplot(3,2,2);
plot(t(:),Real_F_a(4,:),'-r','LineWidth',2);
hold on
plot(t(:),F_a(4,:),'-.','color',[0,0.55,0],'LineWidth',2);
hold on 
legend('$f_{\alpha}$','$\hat{f_{\alpha}}$' ,'fontsize',13,'FontName','Times New Roman','Interpreter','latex');
% xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
ylabel('(rad/s^2)','fontsize',10,'FontName','Times New Roman','Interpreter','tex')
grid on;
ax841.GridLineStyle = '--';
ax841.GridAlpha =0.5;
ax841.GridColor = 'k';
ax841.FontSize= 12;
ax841.FontName= 'Times New Roman';
% ax241.FontAngle= 'italic';
hold off
box on

%--------------------------------------------------------------------------
ax85 = subplot(3,2,4);
plot(t(:),Real_F_a(5,:),'-r','LineWidth',2);
hold on
plot(t(:),F_a(5,:),'-.','color',[0,0.55,0],'LineWidth',2);
hold on
legend('$f_{\beta}$','$\hat{f_{\beta}}$' ,'fontsize',13,'FontName','Times New Roman','Interpreter','latex');
% xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
ylabel('(rad/s^2)','fontsize',10,'FontName','Times New Roman','Interpreter','tex')
grid on;
ax85.GridLineStyle = '--';
ax85.GridAlpha =0.5;
ax85.GridColor = 'k';
ax85.FontSize= 12;
ax85.FontName= 'Times New Roman';
% ax25.FontAngle= 'italic';
hold off
box on
%--------------------------------------------------------------------------
ax86 = subplot(3,2,6);
plot(t(:),Real_F_a(6,:),'-r','LineWidth',2);
hold on
plot(t(:),F_a(6,:),'-.','color',[0,0.55,0],'LineWidth',2);
hold on
legend('$f_{\gamma}$','$\hat{f_{\gamma}}$' ,'fontsize',13,'FontName','Times New Roman','Interpreter','latex');
xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
ylabel('(rad/s^2)','fontsize',10,'FontName','Times New Roman','Interpreter','tex')
grid on;
ax86.GridLineStyle = '--';
ax86.GridAlpha =0.5;
ax86.GridColor = 'k';
ax86.FontSize= 12;
ax86.FontName= 'Times New Roman';
% ax26.FontAngle= 'italic';
hold off
box on

figure(9)
% -------------------------------------------------------------------------
ax91 = subplot(1,1,1);


plot3(Xd(1,:),Xd(2,:),Xd(3,:),'-r','LineWidth',2);
hold on
plot3(X(1,:),X(2,:),X(3,:),'--b','LineWidth',2);
hold off
% title('The reference and actual path in X-Y-Z space','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
xlabel('x~(m)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex');
ylabel('y~(m)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex');
zlabel('z~(m)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex');
legend('Desired path','Actual path','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
ax91.GridLineStyle = '--';
ax91.GridAlpha =0.5;
ax91.GridColor = 'k';
ax91.FontSize= 12;
ax91.FontName= 'Times New Roman';
grid on;
hold off
box on