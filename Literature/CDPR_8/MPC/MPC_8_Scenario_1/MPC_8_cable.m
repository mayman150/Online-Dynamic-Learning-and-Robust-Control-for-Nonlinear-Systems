%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Simulation of the MPC controller for the 6DOF CDPR %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------- ==> Sampling time == T == 0.001 => 1000 Hz
%% ---------------- ==> Simulation time == tfinal == 20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear
clc
%%
% options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',1500);
% options = optimoptions(@fmincon,'Algorithm','trust-region-reflective','MaxIterations',1500);
% options = optimoptions(@fmincon,'Algorithm','sqp-legacy','MaxIterations',1500);
options = optimoptions(@fmincon,'Algorithm','active-set','MaxIterations',1500);

h_c = 10;
h_p = 14;
tfinal = 20;
T = 0.001;
t = 0:T:tfinal;
n = length(t) - h_p;
%% Initialization of the Robot --------------------------------------------
n_a = 6;              %% number of DoF
n_b = 8;              %% number of actuators
%% system initalization


d_x = [0; 0; 0; 0; 0; 0];

x = [0; 0; 0; 0; 0; 0];
d_X = d_x;
X = x;

U_t = ones(n_b*h_c,1);

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
%%
J_hat = zeros(n_b,n_a);
%% ------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
o_t = tic;
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

nn = 8; %number of cables 
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
G = repmat([xG y z rx ry rz],nn,1);
Rotation = [cos(G(1,5))*cos(G(1,6)), cos(G(1,6))*sin(G(1,4))*sin(G(1,5))-cos(G(1,4))*sin(G(1,6)), cos(G(1,4))*cos(G(1,6))*sin(G(1,5))+sin(G(1,4))*sin(G(1,6));
            cos(G(1,5))*sin(G(1,6)), cos(G(1,4))*cos(G(1,6))+sin(G(1,4))*sin(G(1,5))*sin(G(1,6)), -cos(G(1,6))*sin(G(1,4))+cos(G(1,4))*sin(G(1,5))*sin(G(1,6));
            -sin(G(1,5)), cos(G(1,5))*sin(G(1,4)), cos(G(1,4))*cos(G(1,5))]; 
Ei(:,1:3) = (Rotation*Eini(:,1:3).').';
Bi = Ei + G(:,1:3);

AiBi = Bi - Ai;
Li = sqrt(dot(AiBi.',AiBi.'));
for ii=1:nn
    Si(ii,:)=AiBi(ii,:)/Li(ii);
    EixSi(ii,:)= cross(Ei(ii,:),Si(ii,:));
    J(ii,:)=[Si(ii,:) EixSi(ii,:)] ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Uncertain kinematic
nn=8; %number of cables 
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
for ii=1:nn
    Si(ii,:)=AiBi(ii,:)/Li(ii);
    EixSi(ii,:)= cross(Ei(ii,:),Si(ii,:));
    J_hat(ii,:)=[Si(ii,:) EixSi(ii,:)] ;
end

del_J = J_hat - J;
% H4tic = tic;
options = optimset('display','off');
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
C_hat = 1.1*C;
C_matrix = [ zeros(3,3),           zeros(3,3); 
             zeros(3,3),           E_matrix'*(I_p*d_E + cross(E_matrix,I_p*E_matrix)) ];
del_C = C_hat - C;
G = [0;0;m*g;0;0;0];
G_hat = 1.1*G;
del_G = G_hat - G;
%% desired trajectory -----------------------------------------------------
xd = [0.3 + 1.5*exp(-r*T) - 1.8*exp(-r*T/1.2);
      0.3 + 1.5*exp(-r*T) - 1.8*exp(-r*T/1.2);
      0;
      0.2 + exp(-r*T) - 1.2* exp(-r*T/1.2);
      0.2 + exp(-r*T) - 1.2*exp(-r*T/1.2);
      0];
Xd(:,r) = xd;
d_xd = [ - 1.5*exp(-r*T) + 1.8/1.2*exp(-r*T/1.2);
         - 1.5*exp(-r*T) + 1.8/1.2*exp(-r*T/1.2);
         0;
         - exp(-r*T) + 1*exp(-r*T/1.2);
         - exp(-r*T) + 1*exp(-r*T/1.2);
         0];
dd_Xd = [ + 1.5*exp(-r*T) - 1.8/(1.2*1.2)*exp(-r*T/1.2);
          + 1.5*exp(-r*T) - 1.8/(1.2*1.2)*exp(-r*T/1.2);
          0;
          + exp(-r*T) - 1/1.2*exp(-r*T/1.2);
          + exp(-r*T) - 1/1.2*exp(-r*T/1.2);
          0];
e = xd - x;
E(:,r) = e;
d_e = d_xd - d_x;
%% control law
u0 = U_t(1:8,1);
U0(:,r) = u0;%#ok
%% external disturbance & friction
f_1 = -1+0.5^0*cos(3^0*pi*t(r))+0.5^1*cos(3^1*pi*t(r))+0.5^2*cos(3^2*pi*t(r));%+0.5^3*cos(3^3*pi*t(r))+0.5^4*cos(3^4*pi*t(r))+0.5^5*cos(3^5*pi*t(r))+0.5^6*cos(3^6*pi*t(r))+0.5^7*cos(3^7*pi*t(r))+0.5^8*cos(3^8*pi*t(r))+0.5^9*cos(3^9*pi*t(r))+0.5^10*cos(3^10*pi*t(r))+0.5^11*cos(3^11*pi*t(r))+0.5^12*cos(3^12*pi*t(r))+0.5^13*cos(3^13*pi*t(r))+0.5^14*cos(3^14*pi*t(r))+0.5^15*cos(3^15*pi*t(r))+0.5^16*cos(3^16*pi*t(r))+0.5^17*cos(3^17*pi*t(r))+0.5^18*cos(3^18*pi*t(r))+0.5^19*cos(3^19*pi*t(r))+0.5^20*cos(3^20*pi*t(r)); 

 if r<50000

T_d = [20*(x(1))^2+15*(d_x(1))^2, 0, 5*cos(d_x(3)), 0, 5*f_1, 0]' ;
 elseif r>=50000&& r<= 70000
         Q = (1.2 + exp(-r*T+5) - 1.2*exp((-r*T+5)/1.2));

T_d = [20*(x(1))^2+15*(d_x(1))^2, 10*sin(x(2)), 10*f_1, 0, 0, 0]' +J1' *(Q*f_w);
else

 T_d = [20*(x(1))^2+15*(d_x(3))^2+2*cos(d_x(1)), 0,0, 0, 10*f_1, -25*x(6)]' + J1' * (Q*f_w);
end
tau_d = T_d; %[0.0707200027595757;0.0714024021283772;2.86065083265569;0.00145770971816599;0.00539846283852691;0.404274427208022];
F = 0*sin(d_x);
%% dynamic equation of cable robot ----------------------------------------
dd_x = M \ ( - J'*u0 - C - G - F - tau_d); 
d_x = d_x + dd_x*T;
x = x + d_x*T;
X(:,r) = x;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% MPC approach %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cost function parameters -----------------------------------------------
K_y = diag([5*10^7, 0.7*10^7, 5*10^7, 0.4*10^6, 5*10^6, 0.4*10^6, 0.6*10^(-3), 0.6*10^(-3), 0.6*10^(-3), 0.6*10^(-3), 0.6*10^(-3), 0.6*10^(-3)]); % size y_t

K_Y = zeros(h_p*2*n_a);
for i = 1:h_p
    K_Y(2*n_a*(i-1)+1:2*n_a*i, 2*n_a*(i-1)+1:2*n_a*i) = K_y;
end

K_U = 3.3* 10^(-2)* eye(h_c*n_b); % K_U = 3.3* 10^(-2)* eye(h_c*n_b) % scalar * I(h_c*m)

Q = zeros(h_c*n_b);  %% Its dimention ==> m*hc × m*h_c
for i = 1:h_c
    Q(n_b*(i-1)+1:n_b*i, n_b*(i-1)+1:n_b*i) = eye(n_b);
    if i>1
    Q(n_b*(i-1)+1:n_b*i, n_b*(i-2)+1:n_b*(i-1)) = -eye(n_b);
    end
end

z = [u0;
     zeros((h_c-1)*n_b,1)];   %% Its dimention ==> m*hc × 1
k_Du = 0.1; % scalar
K_D = k_Du* eye(h_c*n_b); %% k_Du * I(h_c*m)
%% Discontinuous system matrices ------------------------------------------
A_t = [eye(n_a)    T*eye(n_a);  
       zeros(n_a)  eye(n_a)-T*M_hat^(-1)*C_matrix]; 
B_t = [zeros(n_a, n_b);
       -T*M_hat^(-1)*J_hat']; 
%% Constructing matrices for Y equation -----------------------------------
D_t = zeros(h_p*2*n_a,2*n_a);
for i=1:h_p
    D_t(2*n_a*(i-1)+1:2*n_a*i,:) = A_t^i;
end
E_t = zeros(h_p*2*n_a,h_c*n_b);
for i=1:h_p
    E_t(2*n_a*(i-1)+1:2*n_a*i, 1:n_b) = A_t^(i-1)*B_t;
    if i>1
    E_t(2*n_a*(i-1)+1:2*n_a*i, n_b+1:h_c*n_b) = E_t(2*n_a*(i-2)+1:2*n_a*(i-1), 1:(h_c-1)*n_b);
    end
end
v_t = [zeros(n_a, 1);
       -T*M_hat^(-1)*G_hat];  % [7] page 6
F_t = zeros(2*n_a*h_p,1);
sum = zeros(2*n_a,1);
for i=1:h_p
    sum = sum + A_t^(i-1)*v_t;
    F_t(2*n_a*(i-1)+1:2*n_a*i,1) = sum;
end
%% Calculating Y ----------------------------------------------------------
y_t = [x; d_x];                   %% [7] page 6
%Y_t = D_t*y_t + E_t*U_t + F_t;    %% equation 4 of reference [15]
%% 
for i =1:h_p
    Y_d(2*n_a*(i-1)+1:2*n_a*i,1) = [0.3 + 1.5*exp(-t(r+1)) - 1.8*exp(-t(r+1)/1.2);              
                                0.3 + 1.5*exp(-t(r+1)) - 1.8*exp(-t(r+1)/1.2);
                                0;
                                0.2 + exp(-t(r+1)) - 1.2*exp(-t(r+1)/1.2);
                                0.2 + exp(-t(r+1)) - 1.2*exp(-t(r+1)/1.2);
                                0;
                                - 1.5*exp(-t(r+1)) + 1.8/1.2*exp(-t(r+1)/1.2);
                                - 1.5*exp(-t(r+1)) + 1.8/1.2*exp(-t(r+1)/1.2);
                                0;
                                - exp(-t(r+1)) + 1*exp(-t(r+1)/1.2);
                                - exp(-t(r+1)) + 1*exp(-t(r+1)/1.2);
                                0]; 
end
%% Matrices of the cost function ------------------------------------------
H_c = E_t'*K_Y*E_t + K_U + Q'*K_D*Q;
dT = (D_t*y_t + F_t - Y_d)'*K_Y*E_t - z'*K_D*Q; %% d transpose

%% Optimization problem ---------------------------------------------------
U_t = fmincon(@(xn) xn'*H_c*xn + 2*dT*xn, zeros(n_b*h_c,1),[],[],[],[],1*ones(n_b*h_c,1),200*ones(n_b*h_c,1),[],options);

%% Time calculation -------------------------------------------------------
h3 = toc(H3tic);
H3(r) = h3; 
end
overall_time = toc(o_t)
%% Ploting ----------------------------------------------------------------
Ploting_file(X,Xd,E,U0,t,n);