%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Simulation of the proposed controller for the CDPR %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------- ==> Sampling time == T == 0.0005
%% ---------------- ==> Simulation time == tfinal == 20
%% ---------------- ==> Simulation with cylindrical path without kinematic uncertainty and disturbance for 12-cables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear;
clc;
%%
tfinal = 20;
T = 0.0005;
t = 0:T:tfinal;
n = length(t);
%% system initalization
umin = 1;
umax = 200;
us = (umax-umin)/2; 
up = (umax+umin)/2;
tau = umin*ones(12,1);
u = umin*ones(12,1);
tau_d = tau;
d_x = [0.2, 0, 0.02, 0, 0, 0]';
x = [0, 0.1, 0, 0, 0, 0]'; 
dd_x_d = zeros(6,1);
%% Control Parameters
Ks = 30*diag([1,1,2,1,1,1]);
K_pp = 100*diag([1,1,2,1,1,1]);
K_sw = 2*diag([1,1,1,1,1,1]);
%% RL Parameters
rng(42);
hiddenlayer_a = 3000;  
u_a = 100*randn(hiddenlayer_a,1);%
hiddenlayer_c = 64; % 50
u_c = 10*randn(hiddenlayer_c,1);%
eta_c = 10*ones(hiddenlayer_c,1);%
eta_a = 500*ones(hiddenlayer_a,1);%
K_a = 5*diag([1, 1, 1, 1, 1, 1]);
Q = 10*diag([1, 1, 1, 1, 1, 1]); %
R = 0.01*diag([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]);%
psi = 1000;%
alpha_c = 0.001*diag([1, 1, 1, 0.1, 0.01, 0.001]);%
alpha_a = 0.0001*diag([1, 1, 1, 1, 1, 1]);%
beta_c = 0.01*diag([1, 1, 1, 1, 1, 1]);
beta_a = 0.01*diag([1, 1, 1, 1, 1, 1]);
%% Neural networks initialization
%% Weights of RBF Network of notation a
Wt = zeros(6*hiddenlayer_a + 6*hiddenlayer_c,1);%5*randn(6*hiddenlayer_a + 6*hiddenlayer_c,1);
what_a1 = Wt(1:hiddenlayer_a,1);
what_a2 = Wt(hiddenlayer_a+1:2*hiddenlayer_a,1);
what_a3 = Wt(2*hiddenlayer_a+1:3*hiddenlayer_a,1);
what_a4 = Wt(3*hiddenlayer_a+1:4*hiddenlayer_a,1);
what_a5 = Wt(4*hiddenlayer_a+1:5*hiddenlayer_a,1);
what_a6 = Wt(5*hiddenlayer_a+1:6*hiddenlayer_a,1);
what_a = [what_a1,what_a2,what_a3,what_a4,what_a5,what_a6];
%% Weights of RBF Network of notation c
what_c1 = Wt(6*hiddenlayer_a+1:6*hiddenlayer_a+hiddenlayer_c,1);
what_c2 = Wt(6*hiddenlayer_a+hiddenlayer_c+1:6*hiddenlayer_a+2*hiddenlayer_c,1);
what_c3 = Wt(6*hiddenlayer_a+2*hiddenlayer_c+1:6*hiddenlayer_a+3*hiddenlayer_c,1);
what_c4 = Wt(6*hiddenlayer_a+3*hiddenlayer_c+1:6*hiddenlayer_a+4*hiddenlayer_c,1);
what_c5 = Wt(6*hiddenlayer_a+4*hiddenlayer_c+1:6*hiddenlayer_a+5*hiddenlayer_c,1);
what_c6 = Wt(6*hiddenlayer_a+5*hiddenlayer_c+1:6*hiddenlayer_a+6*hiddenlayer_c,1);
what_c = [what_c1,what_c2,what_c3,what_c4,what_c5,what_c6];
%% ------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
Xd = [0; 0; 0; 0; 0; 0];
E = [0; 0; 0; 0; 0; 0];
%% Rotation matrices
R_1 = [cos(0)*cos(-pi), cos(-pi)*sin(0)*sin(0)-cos(0)*sin(-pi), cos(0)*cos(-pi)*sin(0)+sin(0)*sin(-pi);
      cos(0)*sin(-pi), cos(0)*cos(-pi)+sin(0)*sin(0)*sin(-pi), -cos(-pi)*sin(0)+cos(0)*sin(0)*sin(-pi);
      -sin(0), cos(0)*sin(0), cos(0)*cos(0)];
R_2 = [cos(0)*cos(+pi), cos(+pi)*sin(0)*sin(0)-cos(0)*sin(+pi), cos(0)*cos(+pi)*sin(0)+sin(0)*sin(+pi);
      cos(0)*sin(+pi), cos(0)*cos(+pi)+sin(0)*sin(0)*sin(+pi), -cos(+pi)*sin(0)+cos(0)*sin(0)*sin(+pi);
      -sin(0), cos(0)*sin(0), cos(0)*cos(0)];
%% ------------------------------Main Loop----------------------------------------------> 
tic
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
    J_d(ii,:)=[Si(ii,:) EixSi(ii,:)] ;
end
delta_J = J' - J_d';
%% CDPR parameters --------------------------------------------------------  
g = 9.8;
m = 2.5;
I_p = [0.212 0 0;0 0.225 0;0 0 0.03];            
%% 
E_matrix = [cos(x(5))*cos(x(6)), -sin(x(6)), 0; cos(x(5))*sin(x(6)), cos(x(6)), 0; -sin(x(5)), 0, 1];
d_E = [ -d_x(5)*sin(x(5))*cos(x(6))-d_x(6)*cos(x(5))*sin(x(6)), -d_x(6)*cos(x(6)), 0; 
        -d_x(5)*sin(x(5))*sin(x(6))+d_x(6)*cos(x(5))*cos(x(6)), -d_x(6)*sin(x(6)), 0;
        -d_x(5)*cos(x(5))   ,  0               , 0 ];
%% matrices of the model --------------------------------------------------
M = [m*eye(3) zeros(3); zeros(3) E_matrix'*I_p*E_matrix];
G  = [0;0;m*g;0;0;0];
d_Angles = [d_x(4); d_x(5); d_x(6)];
C_x = [ zeros(3,1); E_matrix'*(I_p*d_E*d_Angles + cross(E_matrix*d_Angles,I_p*E_matrix*d_Angles)) ];
C_matrix = [ zeros(3,3),           zeros(3,3); 
             zeros(3,3),           E_matrix'*(I_p*d_E + cross(E_matrix,I_p*E_matrix)) ];
%% desired trajectory -----------------------------------------------------
xd = 1/2 * [0.2*sin(2*r*T);
      0.2*cos(2*r*T);
      0.2*sin(0.2*r*T);
      0;
      0;
      0];
Xd(:,r) = xd;
d_xd = 1/2 * [ 0.2*2*cos(2*r*T);
         -0.2*2*sin(2*r*T);
         0.2*0.2*cos(0.2*r*T);
         0;
         0;
         0];
dd_xd = 1/2 * [ -0.2*2*2*sin(2*r*T);
          -0.2*2*2*cos(2*r*T);
          -0.2*0.2*0.2*sin(0.2*r*T);
          0;
          0;
          0];
Xd(:,r) = xd;
%% Tracking errors 
e_1 = xd - x;
E(:,r) = e_1;
e_2 = d_xd - d_x; 
%% sliding surface and parameters -----------------------------------------------------
s = e_2 + Ks*e_1;
%% Determine the features of the RBF
z1 = e_1;
d_z1 = e_2;
%% RBF Network of notation c
Z_c = z1;
d_Z_c = d_z1;
eta_c_p2 = eta_c.^2;
h_c = zeros(hiddenlayer_c,1);
for i=1:hiddenlayer_c
h_c(i,1) = exp(-(Z_c - u_c(i,1))'*(Z_c - u_c(i,1))/eta_c_p2(i,1));
end
%% Cost function of RBF network for critic
phi_t = s'*Q*s + tau_d'*R*tau_d;
%% Adaptation law of RBF network for critic   
gradient_h_c = zeros(hiddenlayer_c,6);
for i=1:hiddenlayer_c
gradient_h_c(i,:) = (1/eta_c_p2(i,1))*h_c(i,1)*( -2*(Z_c - u_c(i,1))' );
end
lambda = - (h_c/psi) + gradient_h_c*d_Z_c;     %% Learning rate
d_what_c = - ((phi_t + what_c'*lambda)*lambda')'*alpha_c - what_c*beta_c*alpha_c; 
what_c = what_c + d_what_c*T;
%% RBF Network for actor
Z_a = [x; d_x; xd; d_xd; dd_xd; tau_d];
eta_a_p2 = eta_a.^2;
h_a = zeros(hiddenlayer_a,1);
for i=1:hiddenlayer_a
h_a(i,1) = exp(-(Z_a - u_a(i,1))'*(Z_a - u_a(i,1))/eta_a_p2(i,1));
end
%% Output of RBF Network for critic
I_hat = what_c'*h_c;  %% 6*1
%% Adaptation law of RBF network for actor
d_what_a = - ((what_a'*h_a - alpha_a\s + K_a*I_hat)*h_a')'*alpha_a - what_a*beta_a*alpha_a; %% what_a^T : 2*L 
what_a = what_a + d_what_a*T;
%% Output of RBF Network for actor
f_a = what_a'*h_a;  %% 6*1
F_a(:,r) = f_a;
%% Delta of u
delta_u = tau - u;
%% Disturbance definition
f_1 = -1+0.5^0*cos(3^0*pi*t(r))+0.5^1*cos(3^1*pi*t(r))+0.5^2*cos(3^2*pi*t(r));%+0.5^3*cos(3^3*pi*t(r))+0.5^4*cos(3^4*pi*t(r))+0.5^5*cos(3^5*pi*t(r))+0.5^6*cos(3^6*pi*t(r))+0.5^7*cos(3^7*pi*t(r))+0.5^8*cos(3^8*pi*t(r))+0.5^9*cos(3^9*pi*t(r))+0.5^10*cos(3^10*pi*t(r))+0.5^11*cos(3^11*pi*t(r))+0.5^12*cos(3^12*pi*t(r))+0.5^13*cos(3^13*pi*t(r))+0.5^14*cos(3^14*pi*t(r))+0.5^15*cos(3^15*pi*t(r))+0.5^16*cos(3^16*pi*t(r))+0.5^17*cos(3^17*pi*t(r))+0.5^18*cos(3^18*pi*t(r))+0.5^19*cos(3^19*pi*t(r))+0.5^20*cos(3^20*pi*t(r)); 

if r<50000

T_d = [20*(x(1))^2+15*(d_x(1))^2, 0, 5*cos(d_x(3)), 0, 5*f_1, 0]' ;
 elseif r>=50000&& r<= 70000
         Q = (1.2 + exp(-r*T+5) - 1.2*exp((-r*T+5)/1.2));

T_d = [20*(x(1))^2+15*(d_x(1))^2, 10*sin(x(2)), 10*f_1, 0, 0, 0]' +J1' *(Q*f_w);
else

 T_d = [20*(x(1))^2+15*(d_x(3))^2+2*cos(d_x(1)), 0,0, 0, 10*f_1, -25*x(6)]' + J1' * (Q*f_w);
end
N = C_x + G + T_d + J_d'*delta_u + delta_J*tau;
Real_f_a = N + M*dd_xd + M*Ks*e_2 + C_matrix*s;
Real_F_a(:,r) = Real_f_a;
%% Unconstrained control law ----------------------------------------------
f_w =  f_a + K_pp*s + K_sw*tanh(s);  %
u = - pinv(J_d') * ( f_w);
%% Input saturation via tanh
tau = us*tanh(((u)-[up,up,up,up,up,up,up,up,up,up,up,up]')/us)+up;
tau_d = tau;
Tau(:,r) = tau;%#ok

%% dynamic equation of cable robot ----------------------------------------
dd_x = -M \ (C_x + G + T_d) - (M \ J'*tau) ;%
d_x = d_x + dd_x*T;
x = x + d_x*T;
X(:,r) = x;
dd_x_d = dd_x;
h3 = toc(H3tic);
H3(r) = h3; 
end
toc
perf1 = [rms(E(1,:)); rms(E(2,:)); rms(E(3,:)); rms(E(4,:)); rms(E(5,:)); rms(E(6,:))];
%% plot -------------------------------------------------------------------
Ploting_file(X, Xd, E, Tau, t, Real_F_a, F_a)
 
