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
T = 0.0005;
t = 0:T:tfinal;
n = length(t);

%% system initalization
% d_x = [(1.8/1.2)-1.5; (1.8/1.2)-1.5; 0; 0; 0; 0]; 
dd_x = zeros(6,1);

d_x = [0.2, 0, 0.02, 0, 0, 0]';
x = [0, 0.1, 0, 0, 0, 0]'; 
X = x;
u = zeros(8,1);
u_sat = zeros(8,1);
%%
f_hat_init = [0; 0; -9.8; 0; 0; 0];
gi_hat_init = [0.144028334071706,0.180501788378900,-0.144028334071706,-0.180501788378900,-0.144028334071706,-0.136168015794609,0.144028334071706,0.136168015794609;-0.323226377626038,0.305586361027437,0.323226377626038,-0.305586361027437,0.323226377626038,-0.327753247319582,-0.323226377626038,0.327753247319582;-0.0837374035300617,-0.0791674510433774,-0.0837374035300617,-0.0791674510433774,0.0837374035300618,0.0791674510433774,0.0837374035300618,0.0791674510433774;2.17243263875160,-2.05387255065366,-2.17243263875160,2.05387255065366,2.17243263875160,-2.05387255065366,-2.17243263875160,2.05387255065366;1.02345715425631,0.967602179419056,-1.02345715425631,-0.967602179419056,1.02345715425631,0.967602179419057,-1.02345715425631,-0.967602179419057;-3.22389003590738,-5.07991144195005,-3.22389003590738,-5.07991144195005,-3.22389003590738,5.07991144195005,-3.22389003590738,5.07991144195005];
I_e = [0;0;0;0;0;0];
I_switching_term = [0; 0; 0; 0; 0; 0];
%%
Xd = [0; 0; 0; 0; 0; 0];
E = [0; 0; 0; 0; 0; 0];
%%
R1 = [cos(0)*cos(-pi), cos(-pi)*sin(0)*sin(0)-cos(0)*sin(-pi), cos(0)*cos(-pi)*sin(0)+sin(0)*sin(-pi);
      cos(0)*sin(-pi), cos(0)*cos(-pi)+sin(0)*sin(0)*sin(-pi), -cos(-pi)*sin(0)+cos(0)*sin(0)*sin(-pi);
      -sin(0), cos(0)*sin(0), cos(0)*cos(0)];
R2 = [cos(0)*cos(+pi), cos(+pi)*sin(0)*sin(0)-cos(0)*sin(+pi), cos(0)*cos(+pi)*sin(0)+sin(0)*sin(+pi);
      cos(0)*sin(+pi), cos(0)*cos(+pi)+sin(0)*sin(0)*sin(+pi), -cos(+pi)*sin(0)+cos(0)*sin(0)*sin(+pi);
      -sin(0), cos(0)*sin(0), cos(0)*cos(0)];
%% ------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
for r = 1:n
H3tic = tic;
% H2tic = tic;
%% kinematic
xG=x(1);
y=x(2);
z=x(3);
rx=x(4);
ry=x(5);
rz=x(6);

n = 8; %number of cables 
a = 0.14; 
b = 0.07;   
h = 0.5;   
fa = 1; 
fb = 2; 
fh = 0;  




Eini(1,1)=+a;
Eini(1,2)=-b;
Eini(2,1)=-a;
Eini(2,2)=+b;
Eini(3,1)=-a;
Eini(3,2)=+b;
Eini(4,1)=+a;
Eini(4,2)=-b;
Eini(5,1)=-a;
Eini(5,2)=+b;
Eini(6,1)=-a;
Eini(6,2)=+b;
Eini(7,1)=+a;
Eini(7,2)=-b;
Eini(8,1)=+a;
Eini(8,2)=-b;
Eini(1:4,3)=h;
Eini(5:8,3)=-h;

Ai(1,1)=-fa;
Ai(1,2)=+fb;
Ai(2,1)=-fa;
Ai(2,2)=-fb;
Ai(3,1)=+fa;
Ai(3,2)=-fb;
Ai(4,1)=+fa;
Ai(4,2)=+fb;
Ai(5,1)=+fa;
Ai(5,2)=-fb;
Ai(6,1)=+fa;
Ai(6,2)=+fb;
Ai(7,1)=-fa;
Ai(7,2)=+fb;
Ai(8,1)=-fa;
Ai(8,2)=-fb;
Ai(1:4,3)=fh;
Ai(5:8,3)=-fh;


 for i=1:4
    Ai(i,:)=R1*Ai(i,:)';
end
for i=5:8
    Ai(i,:)=R2*Ai(i,:)';
end

%% 
G = repmat([xG y z rx ry rz],n,1);
Rotation = [cos(G(1,5))*cos(G(1,6)), cos(G(1,6))*sin(G(1,4))*sin(G(1,5))-cos(G(1,4))*sin(G(1,6)), cos(G(1,4))*cos(G(1,6))*sin(G(1,5))+sin(G(1,4))*sin(G(1,6));
            cos(G(1,5))*sin(G(1,6)), cos(G(1,4))*cos(G(1,6))+sin(G(1,4))*sin(G(1,5))*sin(G(1,6)), -cos(G(1,6))*sin(G(1,4))+cos(G(1,4))*sin(G(1,5))*sin(G(1,6));
            -sin(G(1,5)), cos(G(1,5))*sin(G(1,4)), cos(G(1,4))*cos(G(1,5))]; 
Ei(:,1:3) = (Rotation*Eini(:,1:3).').';
Bi = Ei + G(:,1:3);

AiBi = Bi - Ai;
Li = sqrt(dot(AiBi.',AiBi.'));
for ii=1:n
    Si(ii,:)=AiBi(ii,:)/Li(ii);
    EixSi(ii,:)= cross(Ei(ii,:),Si(ii,:));
    J(ii,:)=[Si(ii,:) EixSi(ii,:)] ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Uncertain kinematic
n=8; %number of cables 
a = 1.1*0.14 ; 
b = 1.1*0.07 ;   
h = 1.1*0.5;   
fa = 1.1*1 ; 
fb = 1.1*2  ; 
fh = 1.1*0 ;  
%%%%%%%%%%%%
Eini(1,1)=+a;
Eini(1,2)=-b;
Eini(2,1)=-a;
Eini(2,2)=+b;
Eini(3,1)=-a;
Eini(3,2)=+b;
Eini(4,1)=+a;
Eini(4,2)=-b;
Eini(5,1)=-a;
Eini(5,2)=+b;
Eini(6,1)=-a;
Eini(6,2)=+b;
Eini(7,1)=+a;
Eini(7,2)=-b;
Eini(8,1)=+a;
Eini(8,2)=-b;
Eini(1:4,3)=h;
Eini(5:8,3)=-h;

Ai(1,1)=-fa;
Ai(1,2)=+fb;
Ai(2,1)=-fa;
Ai(2,2)=-fb;
Ai(3,1)=+fa;
Ai(3,2)=-fb;
Ai(4,1)=+fa;
Ai(4,2)=+fb;
Ai(5,1)=+fa;
Ai(5,2)=-fb;
Ai(6,1)=+fa;
Ai(6,2)=+fb;
Ai(7,1)=-fa;
Ai(7,2)=+fb;
Ai(8,1)=-fa;
Ai(8,2)=-fb;
Ai(1:4,3)=fh;
Ai(5:8,3)=-fh;


 for i=1:4
    Ai(i,:)=R1*Ai(i,:)';
end
for i=5:8
    Ai(i,:)=R2*Ai(i,:)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%
Ei(:,1:3) = (Rotation*Eini(:,1:3).').';
Bi = Ei + G(:,1:3);

AiBi = Bi - Ai;
Li = sqrt(dot(AiBi.',AiBi.'));
for ii=1:n
    Si(ii,:)=AiBi(ii,:)/Li(ii);
    EixSi(ii,:)= cross(Ei(ii,:),Si(ii,:));
    J_hat(ii,:)=[Si(ii,:) EixSi(ii,:)] ;
end

del_J = J_hat - J;
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
dd_xd = 0 * [ -0.2*2*2*sin(2*r*T);
          -0.2*2*2*cos(2*r*T);
          -0.2*0.2*0.2*sin(0.2*r*T);
          0;
          0;
          0];
e = xd - x;
E(:,r) = e;
I_e = I_e + e*T; 
d_e = d_xd - d_x;
dd_e = dd_xd - dd_x;
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
dd_x_r = dd_xd - Lambda*d_x_tilde;
%% ----------- sliding surface
s = d_x_tilde + Lambda*x_tilde;

%% ----------- control law
F = M_hat*dd_x_r + C_matrix_hat*d_x_r + G_hat - K_d*tanh(s);
tau = - pinv(J_hat')*F;
U0(:,r) = tau;

%% ------- Null space calculation
% Xn = fmincon(@(xn) xn(1)^2+xn(2)^2+xn(3)^2+xn(4)^2+xn(5)^2+xn(6)^2+xn(7)^2+xn(8)^2+xn(9)^2+xn(10)^2+xn(11)^2+xn(12)^2,zeros(12,1),[],[],J_hat',zeros(6,1),0.1*ones(12,1),[],[],options);
options = optimset('display','off');
[Xn,fval,exitflag,output]= fmincon(@(xn) xn(1)^2+xn(2)^2+xn(3)^2+xn(4)^2+xn(5)^2+xn(6)^2+xn(7)^2+xn(8)^2,zeros(8,1),[],[],J_hat',zeros(6,1),1e-15*ones(8,1),[],[],options);
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

f_1 = -1+0.5^0*cos(3^0*pi*t(r))+0.5^1*cos(3^1*pi*t(r))+0.5^2*cos(3^2*pi*t(r));%+0.5^3*cos(3^3*pi*t(r))+0.5^4*cos(3^4*pi*t(r))+0.5^5*cos(3^5*pi*t(r))+0.5^6*cos(3^6*pi*t(r))+0.5^7*cos(3^7*pi*t(r))+0.5^8*cos(3^8*pi*t(r))+0.5^9*cos(3^9*pi*t(r))+0.5^10*cos(3^10*pi*t(r))+0.5^11*cos(3^11*pi*t(r))+0.5^12*cos(3^12*pi*t(r))+0.5^13*cos(3^13*pi*t(r))+0.5^14*cos(3^14*pi*t(r))+0.5^15*cos(3^15*pi*t(r))+0.5^16*cos(3^16*pi*t(r))+0.5^17*cos(3^17*pi*t(r))+0.5^18*cos(3^18*pi*t(r))+0.5^19*cos(3^19*pi*t(r))+0.5^20*cos(3^20*pi*t(r)); 

 if r<50000

T_d = [20*(x(1))^2+15*(d_x(1))^2, 0, 5*cos(d_x(3)), 0, 5*f_1, 0]' ;
 elseif r>=50000&& r<= 70000
         Q = (1.2 + exp(-r*T+5) - 1.2*exp((-r*T+5)/1.2));

T_d = [20*(x(1))^2+15*(d_x(1))^2, 10*sin(x(2)), 10*f_1, 0, 0, 0]' +J1' *(Q*f_w);
else

 T_d = [20*(x(1))^2+15*(d_x(3))^2+2*cos(d_x(1)), 0,0, 0, 10*f_1, -25*x(6)]' + J1' * (Q*f_w);
end
tau_d = T_d; 
%% dynamic equation of cable robot ----------------------------------------
dd_x = M \ ( - J'*u_sat - C - G-T_d);
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
Perf1 = [rms(E(1,:)); rms(E(2,:)); rms(E(3,:)); rms(E(4,:)); rms(E(5,:)); rms(E(6,:))];
%% plot -------------------------------------------------------------------
% Ploting_file_6(X,Xd,E,E,E,U0,t);
Ploting_file(X, Xd, E, E, E, U_sat, t)
 