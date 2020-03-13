clear all
close all
clc

T = 0.01; dt = T;
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
v_max = 0.3%0.85;%1;
road_width = 0.6;%2;
circle_r = (road_width - L*cos(pi/4)) / (1-cos(pi/4));

t = 0 : dt : pi *circle_r/(2*v_max) +0.5 ;

%make input reference
num = 20;
v_ref = v_max*[-cos(pi .* [0:1/num:1])/2+0.5 1*ones(1,length(t) - 2*(num+1)) +cos(pi .* [0:1/num:1])/2+0.5];
w_ref = 1*(v_max/circle_r)*[-cos(pi .* [0:1/num:1])/2+0.5 1*ones(1,length(t) -2*(num+1)) +cos(pi .* [0:1/num:1])/2+0.5];

t = 0 : dt : pi *circle_r/(2*v_max) +0.5 +1;

v_ref = [zeros(1,50) v_ref zeros(1,50)];
w_ref = [zeros(1,50) w_ref zeros(1,50)];

%make state reference
theta_ref(1) = 0;
for i = 1 : length(t)-1   
    theta_ref(i+1) = theta_ref(i) + w_ref(i) * dt;    
end

dx_ref = v_ref.*cos(theta_ref);
dy_ref = v_ref.*sin(theta_ref);
x_ref(1) = 0;
y_ref(1) = 0;
for i = 1 : length(t)-1
    x_ref(i+1) = x_ref(i) + dx_ref(i) * dt;
    y_ref(i+1) = y_ref(i) + dy_ref(i) * dt;
end
ddx_ref(1) = 0;
ddy_ref(1) = 0;
ddtheta_ref(1) = 0;
ddx_ref(2:length(dx_ref)) = (dx_ref(2:end)-dx_ref(1:end-1))/T;
ddy_ref(2:length(dx_ref)) = (dy_ref(2:end)-dy_ref(1:end-1))/T;
ddtheta_ref(2:length(dx_ref)) = (w_ref(2:end)-w_ref(1:end-1))/T;

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

%% MPC parameters
N=15; % Batch size
n = size(A,1);
m = size(B,2);
Q = 10*eye(n);
R = 10*eye(m);

%% initial state
x0 = [0 0 0]'; %[x y theta]
u0 = [0 0]'; %[v w]

x0_tilt = x0 - xr(:,1);
u0_tilt = u0 - ur(:,1);

%% Optimization
% options = optimoptions('quadprog',...
%      'Algorithm','interior-point-convex','Display','iter');
options = optimoptions('quadprog',...
     'Algorithm','interior-point-convex');

x_tilt_current = x0_tilt; 
x_state(:,1) = x0;
dx_state(:,1) = [0;0;0];
ddx_state(:,1) = [0;0;0];
for k = 1 : length(t)-N    
    %% MPC formulation
    cur_state = k;
    [Sx,Su,Qb,Rb] = Batch_formulation(A, B, N, Q, R, cur_state);

    H = Su'*Qb*Su + Rb;
%     size(H)
%     rank(H)
    F = Sx'*Qb*Su;
    Y = Sx'*Qb*Sx;

    %% constraints
    [G0,E0,w0] = QP_constraints_ZMP(cur_state, A, B, h2, T, N, D, L, ddxr);
    
    %% QP    
    %make u_tilt
    [u_tilt,fval] = quadprog(H,2*F'*x_tilt_current,G0,w0+E0*x_tilt_current,[],[],[],[],[],options);  
%     [u_tilt,fval] = quadprog(H,2*F'*x_tilt_current,zeros(size(G0)),zeros(size(w0+E0*x_tilt_current)),[],[],[],[],[],options);  
%     u_tilt = -inv(H)*F'*x_tilt_current;
    
    u_tilt_set(:,k) = u_tilt;
    x_tilt_set(:,k) = x_tilt_current;
    
    %make real input to simulate with real plant
    u(:,k) = u_tilt(1:m) + ur(:,k);
    
    %real plant simulation
    dx_state(:,k+1) = [u(1,k)*cos(x_state(3,k)); u(1,k)*sin(x_state(3,k)); u(2,k)];       
    x_state(:,k+1) = x_state(:,k) + dx_state(:,k+1)*T;
    x_tilt_current = x_state(:,k+1) - xr(:,k+1);
    
    %zmp    
    zmp(:,k) = -h2/g*ddx_state(:,k);
    ddx_state(:,k+1) = (dx_state(:,k+1) - dx_state(:,k))/T;    
    
end

%% plot
figure;
plot(xr(1,:),xr(2,:),'r','LineWidth',2)
hold on
plot(x_state(1,:),x_state(2,:),'--b')
title('trajectory')
legend('reference','real')
grid on

figure;
subplot(2,3,1)
plot(t,v_ref,'b','LineWidth',2)
hold on
plot(t(1:end-N),u(1,:),'--k')
plot(t,w_ref,'r','LineWidth',2)
plot(t(1:end-N),u(2,:),'--k')
legend('v ref','v real','w ref','w real')
title('input')
grid on

subplot(2,3,2)
% figure;
plot(t(1:end-N),zmp(1,:),'b')
hold on
% plot([t(1) t(end-N)], [D/2 D/2],'k')
% plot([t(1) t(end-N)], [-D/2 -D/2],'k')
plot(t(1:end), -(h2/g*ddxr(1,:)+D/2),'k')
plot(t(1:end), -(h2/g*ddxr(1,:)-D/2),'k')
title('zmp x')
grid on

subplot(2,3,3)
plot(t(1:end-N),zmp(2,:),'b')
hold on
% plot([t(1) t(end-N)], [L/2 L/2],'k')
% plot([t(1) t(end-N)], [-L/2 -L/2],'k')
plot(t(1:end), -(h2/g*ddxr(2,:)+L/2),'k')
plot(t(1:end), -(h2/g*ddxr(2,:)-L/2),'k')
title('zmp y')
grid on

subplot(2,3,4)
plot(t(1:end-N+1),dx_state(1,:),'b')
hold on
plot(t(1:end-N+1),dx_state(2,:),'r')
plot(t(1:end-N+1),dx_state(3,:),'k')
title('dstate')
legend('x','y','theta')
grid on

subplot(2,3,5)
plot(t(1:end-N),x_tilt_set(1,:),'b')
hold on
plot(t(1:end-N),x_tilt_set(2,:),'r')
plot(t(1:end-N),x_tilt_set(3,:),'k')
title('x tilt')
legend('x','y','theta')
grid on

subplot(2,3,6)
plot(t(1:end-N),u_tilt_set(1,:),'b')
hold on
plot(t(1:end-N),u_tilt_set(2,:),'r')
title('u tilt')
legend('v','w')
grid on

% figure;plot((dx_state(1,2:end)-dx_state(1,1:end-1))/T)