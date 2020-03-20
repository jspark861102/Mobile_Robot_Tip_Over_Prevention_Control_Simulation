clear all
close all
clc

T = 0.01; dt = T;
%% design parameters
N=15; % Batch size
Qs = 10;
Rs = 1;

%1:central, 2:backward, 3:forward constraint
zmp_switch = 1;

%1:nonlinear model, 2:linear model
model_switch = 2;

%% system parameters
%mass
m1 = 17;
m2 = 10;
m = m1 + m2;
g = 9.81;

%wheel radius
r = 0.1;
%wheel base
D = 0.3;
%width
L = 0.4;

%z value from global origin to link1 mass center
h1 = 0.20;
%z value from global origin to link2 mass center
h2 = 0.9;
%z value from joint 0 to link1 mass center
h0c1 = h1;
%z value from joint 1 to link2 mass center
h1c2 = 0.7;

%x value from global origin to link2 mass center
a = 0.09;
%y value from global origin to link2 mass center
b = 0.05;
%moment of inertia
I = 1/12*m1*(D^2 + L^2) + m2*(a^2+b^2); %m2�� lumped mass�� ����

%% reference
v_max = 0.3; %0.85;%1;
road_width = 0.6;%2;
circle_r = (road_width - L*cos(pi/4)) / (1-cos(pi/4));

t = 0 : dt : pi *circle_r/(2*v_max) +0.5 ;

%make input reference
num = 20;
v_ref = v_max*[-cos(pi .* [0:1/num:1])/2+0.5 1*ones(1,length(t) - 2*(num+1)) +cos(pi .* [0:1/num:1])/2+0.5];
w_ref = 1*(v_max/circle_r)*[-cos(pi .* [0:1/num:1])/2+0.5 1*ones(1,length(t) -2*(num+1)) +cos(pi .* [0:1/num:1])/2+0.5];

t = 0 : dt : pi *circle_r/(2*v_max) +0.5 +1; %to add zeros in front and back

v_ref = [zeros(1,50) v_ref zeros(1,50)];
w_ref = [zeros(1,50) w_ref zeros(1,50)];

%make state reference
theta_ref(1) = 0;
for i = 1 : length(t)-1   
    theta_ref(i+1) = theta_ref(i) + w_ref(i) * dt; %forward integration    
end

dx_ref = v_ref.*cos(theta_ref);
dy_ref = v_ref.*sin(theta_ref);
x_ref(1) = 0;
y_ref(1) = 0;
for i = 1 : length(t)-1
    x_ref(i+1) = x_ref(i) + dx_ref(i) * dt; %forward integration
    y_ref(i+1) = y_ref(i) + dy_ref(i) * dt; %forward integration
end
ddx_ref(1:length(dx_ref)-1) = (dx_ref(2:end)-dx_ref(1:end-1))/T;  %forward derivative
ddy_ref(1:length(dx_ref)-1) = (dy_ref(2:end)-dy_ref(1:end-1))/T;  %forward derivative
ddtheta_ref(1:length(dx_ref)-1) = (w_ref(2:end)-w_ref(1:end-1))/T;%forward derivative
ddx_ref(length(dx_ref)) = 0;
ddy_ref(length(dx_ref)) = 0;
ddtheta_ref(length(dx_ref)) = 0;

%define state and input reference
xr = [x_ref; y_ref; theta_ref];
ur = [v_ref; w_ref];

ddxr = [ddx_ref; ddy_ref; ddtheta_ref];

%% State space equation
for i = 1 : length(t)
    A(:,:,i) = eye(3) + [0 0 -v_ref(i)*sin(theta_ref(i))*T;
                         0 0  v_ref(i)*cos(theta_ref(i))*T;
                         0 0  0];
    B(:,:,i) = [cos(theta_ref(i))*T 0;
                sin(theta_ref(i))*T 0;
                0                   T];    
end
% linearization_test;

%% MPC parameters
n = size(A,1);
m = size(B,2);
Q = Qs*eye(n);
R = Rs*eye(m);

%% initial state
x0 = [0 0 0]'; %[x y theta]
x0_tilt = x0 - xr(:,1);

%% Optimization
% options = optimoptions('quadprog',...
%      'Algorithm','interior-point-convex','Display','iter');
options = optimoptions('quadprog',...
     'Algorithm','interior-point-convex');

x_tilt_current = x0_tilt; 
for k = 1 : length(t)-N    
    %% MPC formulation
    cur_state = k;
    [Sx,Su,Qb,Rb] = Batch_formulation(A, B, N, Q, R, cur_state);

    H = Su'*Qb*Su + Rb;
    F = Sx'*Qb*Su;
    Y = Sx'*Qb*Sx;

    %% constraints
%     [G0,E0,w0] = QP_constraints_ZMPm(cur_state, A, B, h2, T, N, D, L, xr, ddxr, zmp_switch);
    [G0,E0,w0] = QP_constraints_ZMP(cur_state, A, B, h2, T, N, D, L, ddxr, zmp_switch);
%     [G0,E0,w0] = QP_constraints_basic(cur_state, A, B, h2, T, N, D, L, ddxr);

%     [u_tilt,fval] = quadprog(blkdiag(H, 0.15*eye(n*N,n*N)),[2*F'*x_tilt_current; zeros(n*N,1)],G0,w0+E0*x_tilt_current,[],[],[],[],[],options);  
    [u_tilt,fval] = quadprog(H,2*F'*x_tilt_current,G0,w0+E0*x_tilt_current,[],[],[],[],[],options);  
%     [u_tilt,fval] = quadprog(H,2*F'*x_tilt_current,zeros(size(G0)),zeros(size(w0+E0*x_tilt_current)),[],[],[],[],[],options);  
%     u_tilt = -inv(H)*F'*x_tilt_current;
    
    u_tilt_set(:,k) = u_tilt;
    x_tilt_set(:,k) = x_tilt_current;     
    
    u(:,k) = u_tilt(1:m) + ur(:,k);
    x_state(:,k) = x_tilt_current + xr(:,k);
    
    if model_switch == 1 %nonlinear model
        dx_state(:,k) = [u(1,k)*cos(x_state(3,k)); u(1,k)*sin(x_state(3,k)); u(2,k)];  %forward derivative     
        x_state_next = x_state(:,k) + dx_state(:,k)*T;    %forward derivative
        x_tilt_current = x_state_next - xr(:,k+1);%     
        
    elseif model_switch == 2 %linear model    
        x_tilt_current_next = A(:,:,k)*x_tilt_current + B(:,:,k)*u_tilt(1:m);  
        dx_state(:,k) = ((x_tilt_current_next+xr(:,k+1)) - x_state(:,k))/T; %forward derivative
        x_tilt_current = x_tilt_current_next;
    end    
end

%% zmp        
ddx_state(:,1:length(dx_state)-1) = (dx_state(:,2:end) - dx_state(:,1:end-1))/T;    %forward derivative
ddx_state(:,length(dx_state)) = 0;
zmp = -h2/g*ddx_state;
zmp_tilt = -h2/g*(ddx_state - ddxr(:,1:length(ddx_state)));


%% plot
data_plot