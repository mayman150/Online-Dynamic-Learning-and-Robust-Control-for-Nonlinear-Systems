%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Simulation of the proposed controller for the 6DOF CDPR %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------- ==> Sampling time == T == 0.001 => 1000 Hz
%% ---------------- ==> Simulation time == tfinal == 20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear
clc
%%
tfinal = 20;
T = 0.001;
t = 0:T:tfinal;
n = length(t);

%% system initalization
% d_x = [(1.8/1.2)-1.5; (1.8/1.2)-1.5; 0; 0; 0; 0]; 
dd_x = zeros(6,1);

d_x = zeros(6,1);

x = zeros(6,1);
X = x;
u = zeros(12,1);
u_sat = zeros(12,1);
%%
f_hat_init = [0; 0; -9.8; 0; 0; 0];
gi_hat_init = [0.144028334071706,0.180501788378900,-0.144028334071706,-0.180501788378900,-0.144028334071706,-0.136168015794609,0.144028334071706,0.136168015794609;-0.323226377626038,0.305586361027437,0.323226377626038,-0.305586361027437,0.323226377626038,-0.327753247319582,-0.323226377626038,0.327753247319582;-0.0837374035300617,-0.0791674510433774,-0.0837374035300617,-0.0791674510433774,0.0837374035300618,0.0791674510433774,0.0837374035300618,0.0791674510433774;2.17243263875160,-2.05387255065366,-2.17243263875160,2.05387255065366,2.17243263875160,-2.05387255065366,-2.17243263875160,2.05387255065366;1.02345715425631,0.967602179419056,-1.02345715425631,-0.967602179419056,1.02345715425631,0.967602179419057,-1.02345715425631,-0.967602179419057;-3.22389003590738,-5.07991144195005,-3.22389003590738,-5.07991144195005,-3.22389003590738,5.07991144195005,-3.22389003590738,5.07991144195005];
I_e = [0;0;0;0;0;0];
I_switching_term = [0; 0; 0; 0; 0; 0];
%%
Xd = [0; 0; 0; 0; 0; 0];
E = [0; 0; 0; 0; 0; 0];
%%
R_1 = [cos(0)*cos(-pi), cos(-pi)*sin(0)*sin(0)-cos(0)*sin(-pi), cos(0)*cos(-pi)*sin(0)+sin(0)*sin(-pi);
      cos(0)*sin(-pi), cos(0)*cos(-pi)+sin(0)*sin(0)*sin(-pi), -cos(-pi)*sin(0)+cos(0)*sin(0)*sin(-pi);
      -sin(0), cos(0)*sin(0), cos(0)*cos(0)];
R_2 = [cos(0)*cos(+pi), cos(+pi)*sin(0)*sin(0)-cos(0)*sin(+pi), cos(0)*cos(+pi)*sin(0)+sin(0)*sin(+pi);
      cos(0)*sin(+pi), cos(0)*cos(+pi)+sin(0)*sin(0)*sin(+pi), -cos(+pi)*sin(0)+cos(0)*sin(0)*sin(+pi);
      -sin(0), cos(0)*sin(0), cos(0)*cos(0)];
%%
J_hat = zeros(12,6);
%% ------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
o_t = tic;
for r = 1:n
H3tic = tic;
xG = x(1);
y = x(2);
z = x(3);
rx = x(4);
ry = x(5);
rz = x(6);

nn = 12; %number of cables 
a = 0.14; 
b = 0.07;   
h = 0.5;   
fa = 1; 
fb = 2; 
fh = 0;  



%% B1
Eini(1,1)=+a;
Eini(1,2)=-b;
%% B2
Eini(2,1)=-a;
Eini(2,2)=+b;
%% B3
Eini(3,1)=-a;
Eini(3,2)=+b;
%% B4
Eini(4,1)=+a;
Eini(4,2)=-b;
%% B5
Eini(5,1)=-a;
Eini(5,2)=+b;
%% B6
Eini(6,1)=-a;
Eini(6,2)=+b;
%% B7
Eini(7,1)=+a;
Eini(7,2)=-b;
%% B8
Eini(8,1)=+a;
Eini(8,2)=-b;
%% B9
Eini(9,1)=-a;
Eini(9,2)=-b;
%% B10
Eini(10,1)=+a;
Eini(10,2)=+b;
%% B11
Eini(11,1)=+a;
Eini(11,2)=+b;
%% B12
Eini(12,1)=-a;
Eini(12,2)=-b;



Eini(1:4,3) = h;
Eini(9,3) = h;
Eini(11,3) = h;
Eini(5:8,3) = -h;
Eini(10,3) = -h;
Eini(12,3) = -h;
%% A1
Ai(1,1) = -fa;
Ai(1,2) = +fb;
%% A2
Ai(2,1) = -fa;
Ai(2,2) = -fb;
%% A3
Ai(3,1) = +fa;
Ai(3,2) = -fb;
%% A4
Ai(4,1) = +fa;
Ai(4,2) = +fb;
%% A5
Ai(5,1) = +fa;
Ai(5,2) = -fb;
%% A6
Ai(6,1) = +fa;
Ai(6,2) = +fb;
%% A7
Ai(7,1) = -fa;
Ai(7,2) = +fb;
%% A8
Ai(8,1) = -fa;
Ai(8,2) = -fb;
%% A9
Ai(9,1) = +fa;
Ai(9,2) = -fb;
%% A10
Ai(10,1) = +fa;
Ai(10,2) = +fb;
%% A11
Ai(11,1) = -fa;
Ai(11,2) = +fb;
%% A12
Ai(12,1) = -fa;
Ai(12,2) = -fb;

%%
Ai(1:4,3) = fh;
Ai(5:8,3) = -fh;
Ai(9,3) = fh;
Ai(11,3) = fh;
Ai(10,3) = -fh;
Ai(12,3) = -fh;

%%
for i=1:4
    Ai(i,:) = R_1*Ai(i,:)';
end
for i=5:8
    Ai(i,:) = R_2*Ai(i,:)';
end

%% 
G = repmat([xG y z rx ry rz],nn,1);
Rotation = [cos(G(1,5))*cos(G(1,6)), cos(G(1,6))*sin(G(1,4))*sin(G(1,5))-cos(G(1,4))*sin(G(1,6)), cos(G(1,4))*cos(G(1,6))*sin(G(1,5))+sin(G(1,4))*sin(G(1,6));
            cos(G(1,5))*sin(G(1,6)), cos(G(1,4))*cos(G(1,6))+sin(G(1,4))*sin(G(1,5))*sin(G(1,6)), -cos(G(1,6))*sin(G(1,4))+cos(G(1,4))*sin(G(1,5))*sin(G(1,6));
            -sin(G(1,5)), cos(G(1,5))*sin(G(1,4)), cos(G(1,4))*cos(G(1,5))]; 
Ei(:,1:3) = (Rotation*Eini(:,1:3).').';
Bi = Ei + G(:,1:3);

AiBi = Bi - Ai;
Li = sqrt(dot(AiBi.',AiBi.'));
for ii=1:nn
    Si(ii,:) = AiBi(ii,:)/Li(ii);
    EixSi(ii,:) = cross(Ei(ii,:),Si(ii,:));
    J(ii,:) = [Si(ii,:) EixSi(ii,:)] ;
end
%% Kinematic uncertainty
nn = 12; %number of cables 
a = 1.1*0.14; 
b = 1.1*0.07;   
h = 1.1*0.5;   
fa = 1.1*1; 
fb = 1.1*2; 
fh = 1.1*0;   
%%%%%%%%%%%%
%% B1
Eini(1,1)=+a;
Eini(1,2)=-b;
%% B2
Eini(2,1)=-a;
Eini(2,2)=+b;
%% B3
Eini(3,1)=-a;
Eini(3,2)=+b;
%% B4
Eini(4,1)=+a;
Eini(4,2)=-b;
%% B5
Eini(5,1)=-a;
Eini(5,2)=+b;
%% B6
Eini(6,1)=-a;
Eini(6,2)=+b;
%% B7
Eini(7,1)=+a;
Eini(7,2)=-b;
%% B8
Eini(8,1)=+a;
Eini(8,2)=-b;
%% B9
Eini(9,1)=-a;
Eini(9,2)=-b;
%% B10
Eini(10,1)=+a;
Eini(10,2)=+b;
%% B11
Eini(11,1)=+a;
Eini(11,2)=+b;
%% B12
Eini(12,1)=-a;
Eini(12,2)=-b;



Eini(1:4,3)=h;
Eini(9,3)=h;
Eini(11,3)=h;
Eini(5:8,3)=-h;
Eini(10,3)=-h;
Eini(12,3)=-h;
%% A1
Ai(1,1)=-fa;
Ai(1,2)=+fb;
%% A2
Ai(2,1)=-fa;
Ai(2,2)=-fb;
%% A3
Ai(3,1)=+fa;
Ai(3,2)=-fb;
%% A4
Ai(4,1)=+fa;
Ai(4,2)=+fb;
%% A5
Ai(5,1)=+fa;
Ai(5,2)=-fb;
%% A6
Ai(6,1)=+fa;
Ai(6,2)=+fb;
%% A7
Ai(7,1)=-fa;
Ai(7,2)=+fb;
%% A8
Ai(8,1)=-fa;
Ai(8,2)=-fb;
%% A9
Ai(9,1)=+fa;
Ai(9,2)=-fb;
%% A10
Ai(10,1)=+fa;
Ai(10,2)=+fb;
%% A11
Ai(11,1)=-fa;
Ai(11,2)=+fb;
%% A12
Ai(12,1)=-fa;
Ai(12,2)=-fb;

%%
Ai(1:4,3)=fh;
Ai(5:8,3)=-fh;
Ai(9,3)=fh;
Ai(11,3)=fh;
Ai(10,3)=-fh;
Ai(12,3)=-fh;

%%
 for i=1:4
    Ai(i,:)=R_1*Ai(i,:)';
end
for i=5:8
    Ai(i,:)=R_2*Ai(i,:)';
end


%%%%%%%%%%%%%%%%%%%%%%%%%
Ei(:,1:3) = (Rotation*Eini(:,1:3).').';
Bi = Ei + G(:,1:3);

AiBi = Bi - Ai;
Li = sqrt(dot(AiBi.',AiBi.'));
for ii=1:nn
    Si(ii,:)=AiBi(ii,:)/Li(ii);
    EixSi(ii,:)= cross(Ei(ii,:),Si(ii,:));
    J_hat(ii,:)=[Si(ii,:) EixSi(ii,:)] ;
end
delta_J = J' - J_hat';
% h2 = toc(H2tic);
% H2(r) = h2;


% H4tic = tic;
%% CDPR parameters --------------------------------------------------------
g = 9.8;
m = 2.5;                                                      % The Mass of moving platform M
I_p = [0.212 0 0;0 0.225 0;0 0 0.03];                         % The moment of inertia of the MP
%% 
E_matrix = [cos(x(5))*cos(x(6)), -sin(x(6)), 0; cos(x(5))*sin(x(6)), cos(x(6)), 0; -sin(x(5)), 0, 1];
d_E = [ -d_x(5)*sin(x(5))*cos(x(6))-d_x(6)*cos(x(5))*sin(x(6)), -d_x(6)*cos(x(6)), 0; 
        -d_x(5)*sin(x(5))*sin(x(6))+d_x(6)*cos(x(5))*cos(x(6)), -d_x(6)*sin(x(6)), 0;
        -d_x(5)*cos(x(5))                                     ,  0               , 0 ];
%% matrices of the model --------------------------------------------------
M = [m*eye(3) zeros(3); zeros(3) E_matrix'*I_p*E_matrix];
M_hat = 1.1*M;
del_M = M_hat - M;
d_Angles = [d_x(4); d_x(5); d_x(6)];
C = [ zeros(3,1); E_matrix'*(I_p*d_E*d_Angles + cross(E_matrix*d_Angles,I_p*E_matrix*d_Angles)) ];
C_matrix = [ zeros(3,3),           zeros(3,3); 
             zeros(3,3),           E_matrix'*(I_p*d_E + cross(E_matrix,I_p*E_matrix)) ];
C_matrix_hat = 1.1*C_matrix;

C_hat = 1.1*C;
del_C = C_hat - C;
G = [0;0;m*g;0;0;0];
G_hat = 1.1*G;
del_G = G_hat - G;
%% desired trajectory -----------------------------------------------------
xd = [0.3+1.5*exp(-r*T)- 1.8*exp(-(r*T)/1.2);
      0.3+1.5*exp(-r*T)- 1.8*exp(-(r*T)/1.2);
      0;
       0.2+exp(-r*T)-1.2*exp(-(r*T)/1.2);
       0.2+exp(-r*T)-1.2*exp(-(r*T)/1.2);
       0];
Xd(:,r) = xd;
d_xd = [ 1.5*exp(-r*T)-1.8*(1/1.2)*(1/1.2)*exp(-(r*T)/1.2); 
          1.5*exp(-r*T)-1.8*(1/1.2)*(1/1.2)*exp(-(r*T)/1.2);
          0;
          exp(-r*T)-(1/1.2)*exp(-(r*T)/1.2);
          exp(-r*T)-(1/1.2)*exp(-(r*T)/1.2);
          0];
dd_Xd = [ (+ exp(-r*T) - 1/1.2*exp(-r*T/1.2));
          (+ exp(-r*T) - 1/1.2*exp(-r*T/1.2));
          0;
          0;
          0;
          0];
ddd_Xd = [ (- exp(-r*T) + 1/(1.2*1.2)*exp(-r*T/1.2));
           (- exp(-r*T) + 1/(1.2*1.2)*exp(-r*T/1.2));
           0;
           0;
           0;
           0];
e = xd - x;
E(:,r) = e;
I_e = I_e + e*T; 
d_e = d_xd - d_x;
dd_e = dd_Xd - dd_x;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%% Adaptive sliding control %%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% --------- defining parameters
Lambda = 10*diag([1,1,1,1,1,1]);
K_d = 50*diag([1,1,1,1,1,1]);
%% --------- defining errors
x_tilde = x - xd;
d_x_tilde = d_x - d_xd;
%I_x_tilde = I_x_tilde + x_tilde*T;
%x_r = xd - Lambda*I_x_tilde;
d_x_r = d_xd - Lambda*x_tilde;
dd_x_r = dd_Xd - Lambda*d_x_tilde;
%% ----------- sliding surface
s = d_x_tilde + Lambda*x_tilde;

%% ----------- control law
F = M_hat*dd_x_r + C_matrix_hat*d_x_r + G_hat - K_d*tanh(s);
tau = - pinv(J_hat')*F;
U0(:,r) = tau;

%% ------- Null space calculation
% Xn = fmincon(@(xn) xn(1)^2+xn(2)^2+xn(3)^2+xn(4)^2+xn(5)^2+xn(6)^2+xn(7)^2+xn(8)^2+xn(9)^2+xn(10)^2+xn(11)^2+xn(12)^2,zeros(12,1),[],[],J_hat',zeros(6,1),0.1*ones(12,1),[],[],options);
options = optimset('display','off');
[Xn,fval,exitflag,output]= fmincon(@(xn) xn(1)^2+xn(2)^2+xn(3)^2+xn(4)^2+xn(5)^2+xn(6)^2+xn(7)^2+xn(8)^2+xn(9)^2+xn(10)^2+xn(11)^2+xn(12)^2,zeros(12,1),[],[],J_hat',zeros(6,1),1e-15*ones(12,1),[],[],options);
if exitflag == 0
    break
end
%% -------- IFC Block
umin = 1;
for i=1:8
if tau(i)<umin
    
     R = (umin-tau(i))/abs(Xn(i));
     if(Xn(i)>0)
     tau = tau+R*Xn;
     else 
     tau = tau-R*Xn; 
     end   
end
end


u_sat = tau;
U_sat(:,r) = u_sat; %#ok


%% dynamic equation of cable robot ----------------------------------------
dd_x = M \ ( - J'*u_sat - C - G);
% dd_x = M \ ( - J'*u_sat - C_matrix*d_x - G);
d_x = d_x + dd_x*T;
x = x + d_x*T;
X(:,r) = x;
% h = toc(H1tic);
% H(r) = h; 

h3 = toc(H3tic);
H3(r) = h3; 
end
toc

%% plot -------------------------------------------------------------------
Perf1 = [rms(E(1,:)); rms(E(2,:)); rms(E(3,:)); rms(E(4,:)); rms(E(5,:)); rms(E(6,:))]
%% plot -------------------------------------------------------------------
% Ploting_file_6(X,Xd,E,E,E,U0,t);
Ploting_file(X, Xd, E, E, E, U_sat, t)
 
