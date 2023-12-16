function Ploting_file(X, Xd, E, U0, t,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% error ploting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
% -------------------------------------------------------------------------
ax11 = subplot(2,1,1);
hold on
plot(t(1:n),E(1,:),'-r','LineWidth',2);
plot(t(1:n),E(2,:),'--b','LineWidth',2);
plot(t(1:n),E(3,:),'-.','color',[0,0.5,0],'LineWidth',2);
hold off
legend('$e_x$','$e_y$','$e_z$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
% title('Position errors','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
ylabel('Position Error (m)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
grid on;
ax11.GridLineStyle = '--';
ax11.GridAlpha = 0.5;
ax11.GridColor = 'k';
ax11.FontSize = 12;
ax11.FontName = 'Times New Roman';
ax11.FontAngle = 'italic';
box on
% -------------------------------------------------------------------------
ax12=subplot(2,1,2);
hold on
plot(t(1:n),E(4,:),'-r','LineWidth',2);
plot(t(1:n),E(5,:),'--b','LineWidth',2);
plot(t(1:n),E(6,:),'-.','color',[0,0.5,0],'LineWidth',2);
hold off
legend('$e_{\alpha}$','$e_\beta$','$e_\gamma$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
% title('Orientation errors','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
ylabel('Orientation Error (rad)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
grid on;
ax12.GridLineStyle = '--';
ax12.GridAlpha = 0.5;
ax12.GridColor = 'k';
ax12.FontSize = 12;
ax12.FontName = 'Times New Roman';
ax12.FontAngle = 'italic';
box on
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Postion tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
% -------------------------------------------------------------------------
ax21=subplot(3,1,1);
plot(t(1:n),Xd(1,:),'-r','LineWidth',2);
hold on
plot(t(1:n),X(1,:),'--b','LineWidth',2);
legend('$x_{d}$','x','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
% xlabel('Time [sec]')
%title('Desired and actual position','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
ylabel('m','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
grid on;
ax21.GridLineStyle = '--';
ax21.GridAlpha =0.5;
ax21.GridColor = 'k';
ax21.FontSize= 12;
ax21.FontName= 'Times New Roman';
ax21.FontAngle= 'italic';
hold off
box on
%--------------------------------------------------------------------------
ax22=subplot(3,1,2);
plot(t(1:n),Xd(2,:),'-r','LineWidth',2);
hold on
plot(t(1:n),X(2,:),'--b','LineWidth',2);
legend('$y_{d}$','y','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
ylabel('m','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
grid on;
ax22.GridLineStyle = '--';
ax22.GridAlpha = 0.5;
ax22.GridColor = 'k';
ax22.FontSize= 12;
ax22.FontName= 'Times New Roman';
ax22.FontAngle= 'italic';
hold off
%--------------------------------------------------------------------------
ax23 = subplot(3,1,3);
plot(t(1:n),Xd(3,:),'-r','LineWidth',2);
hold on
plot(t(1:n),X(3,:),'--b','LineWidth',2);
legend('$z_{d}$','$z$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
ylabel('m','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
grid on;
ax23.GridLineStyle = '--';
ax23.GridAlpha = 0.5;
ax23.GridColor = 'k';
ax23.FontSize = 12;
ax23.FontName = 'Times New Roman';
ax23.FontAngle= 'italic';
hold off
box on
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% orientation tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
%--------------------------------------------------------------------------
ax24 = subplot(3,1,1);
plot(t(1:n),Xd(4,:),'-r','LineWidth',2);
hold on
plot(t(1:n),X(4,:),'--b','LineWidth',2);
legend('${\alpha}_{d}$','$\alpha$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
% xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
%title('Desired and actual orientation','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')
ylabel('rad','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
grid on;
ax24.GridLineStyle = '--';
ax24.GridAlpha =0.5;
ax24.GridColor = 'k';
ax24.FontSize= 12;
ax24.FontName= 'Times New Roman';
ax24.FontAngle= 'italic';
hold off
box on
%--------------------------------------------------------------------------
ax25 = subplot(3,1,2);
plot(t(1:n),Xd(5,:),'-r','LineWidth',2);
hold on
plot(t(1:n),X(5,:),'--b','LineWidth',2);
legend('${\beta}_{d}$','$\beta$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
% xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
ylabel('rad','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
grid on;
ax25.GridLineStyle = '--';
ax25.GridAlpha =0.5;
ax25.GridColor = 'k';
ax25.FontSize= 12;
ax25.FontName= 'Times New Roman';
ax25.FontAngle= 'italic';
hold off
box on
%--------------------------------------------------------------------------
ax26 = subplot(3,1,3);
plot(t(1:n),Xd(6,:),'-r','LineWidth',2);
hold on
plot(t(1:n),X(6,:),'--b','LineWidth',2);
legend('${\gamma}_{d}$','$\gamma$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
ylabel('rad','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
grid on;
ax26.GridLineStyle = '--';
ax26.GridAlpha =0.5;
ax26.GridColor = 'k';
ax26.FontSize= 12;
ax26.FontName= 'Times New Roman';
ax26.FontAngle= 'italic';
hold off
Pref = 1.3*[rms(E(1,:)); rms(E(2,:)); rms(E(3,:)); rms(E(4,:)); rms(E(5,:)); rms(E(6,:))]
box on
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  control input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
% -------------------------------------------------------------------------
ax31 = subplot(2,1,1);
hold on
plot(t(1:n),U0(1,:),'--r','LineWidth',2);
plot(t(1:n),U0(2,:),'-.k','LineWidth',2);
plot(t(1:n),U0(3,:),':b','LineWidth',2);
plot(t(1:n),U0(4,:),'-.','Color',[0.64,0.08,0.18],'LineWidth',2);
plot(t(1:n),U0(5,:),'-','Color',[0.8,0.8,0.8],'LineWidth',2);
plot(t(1:n),U0(6,:),'-','Color',[0.47,0.67,0.19],'LineWidth',2);
plot(t(1:n),U0(7,:),'-','Color',[0.49,0.18,0.56],'LineWidth',2);
plot(t(1:n),U0(8,:),'-','Color',[1,0.41,0.16],'LineWidth',2);
hold off
legend('$\tau_{1}$','$\tau_{2}$','$\tau_{3}$','$\tau_{4}$','$\tau_{5}$','$\tau_{6}$','$\tau_{7}$','$\tau_{8}$','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex');
ylabel('Tension (N)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
xlabel('Time (sec)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','latex')
%title('Cable tensions','fontsize',10,'FontName','Times New Roman','Interpreter','latex')
grid on;
ax31.GridLineStyle = '--';
ax31.GridAlpha = 0.5;
ax31.GridColor = 'k';
ax31.FontSize = 12;
ax31.FontName = 'Times New Roman';
ax31.FontAngle = 'italic';
box on
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Consumed time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure(5)
% -------------------------------------------------------------------------
%ax71 = subplot(1,1,1);

%plot(t(1:n),H3(:),'-r','LineWidth',2);
% title('Consumed time for each sample (T = 0.001)','fontsize',10,'FontName','Times New Roman','FontWeight','Bold','Interpreter','Latex')

%grid on;
%ax71.GridLineStyle = '--';
%ax71.GridAlpha = 0.5;
%ax71.GridColor = 'k';
%ax71.FontSize = 12;
%ax71.FontName = 'Times New Roman';
%ax71.FontAngle = 'italic';


