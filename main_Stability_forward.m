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
% ddx_ref(1) = 0;
% ddy_ref(1) = 0;
% ddtheta_ref(1) = 0;
% ddx_ref(2:length(dx_ref)) = (dx_ref(2:end)-dx_ref(1:end-1))/T;  %backward derivative
% ddy_ref(2:length(dx_ref)) = (dy_ref(2:end)-dy_ref(1:end-1))/T;  %backward derivative
% ddtheta_ref(2:length(dx_ref)) = (w_ref(2:end)-w_ref(1:end-1))/T;%backward derivative

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

% x(:,1) = xr(:,1);
% xx(:,1) = xr(:,1);
% u = [sin(2*pi*t);2*sin(2*pi*t)];
% for i = 1 : length(t)
%     x(:,i+1) = A(:,:,i)*x(:,i) + B(:,:,i)*u(:,i);
%     dxx(:,i) = [u(1,i)*cos(xx(3,i)); u(1,i)*sin(xx(3,i)); u(2,i)];
%     xx(:,i+1) = xx(:,i) + dxx(:,i)*T;
% end
% figure;plot(x(2,:),'b')
% hold on
% plot(xx(2,:),'--r')

%% MPC parameters
N=30; % Batch size
n = size(A,1);
m = size(B,2);
Q = 10*eye(n);
R = 1*eye(m);

%1:central, 2:backward, 3:forward
zmp_switch = 1;

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
% dx_state(:,1) = [0;0;0];
% ddx_state(:,1) = [0;0;0];
for k = 1 : length(t)-N    
    %% MPC formulation
    cur_state = k;
    [Sx,Su,Qb,Rb] = Batch_formulation(A, B, N, Q, R, cur_state);

    H = Su'*Qb*Su + Rb;
    F = Sx'*Qb*Su;
    Y = Sx'*Qb*Sx;

    %% constraints
    [G0,E0,w0] = QP_constraints_ZMP(cur_state, A, B, h2, T, N, D, L, ddxr, zmp_switch);
%     [G0,E0,w0] = QP_constraints_basic(cur_state, A, B, h2, T, N, D, L, ddxr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% QP            
%     %make u_tilt
%     [u_tilt,fval] = quadprog(H,2*F'*x_tilt_current,G0,w0+E0*x_tilt_current,[],[],[],[],[],options);  
% %     [u_tilt,fval] = quadprog(H,2*F'*x_tilt_current,zeros(size(G0)),zeros(size(w0+E0*x_tilt_current)),[],[],[],[],[],options);  
% %     u_tilt = -inv(H)*F'*x_tilt_current;
%     
%     u_tilt_set(:,k) = u_tilt;
%     x_tilt_set(:,k) = x_tilt_current;       
%     
%     %make real input to simulate with real plant
%     u(:,k) = u_tilt(1:m) + ur(:,k);
%     
%     %real plant simulation
%     dx_state(:,k) = [u(1,k)*cos(x_state(3,k)); u(1,k)*sin(x_state(3,k)); u(2,k)];  %forward derivative     
%        
%     %x_state
%     x_state(:,k+1) = x_state(:,k) + dx_state(:,k)*T;    %forward derivative
%     
%     %x_tilt for next step qp
%     x_tilt_current = x_state(:,k+1) - xr(:,k+1);
%     
% %     %zmp        
% %     ddx_state(:,k) = (dx_state(:,k+1) - dx_state(:,k))/T;    
% %     zmp(:,k) = -h2/g*ddx_state(:,k);
% %     zmp_tilt(:,k) = -h2/g*(ddx_state(:,k)- ddxr(:,k));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% QP            
    %make u_tilt
    [u_tilt,fval] = quadprog(H,2*F'*x_tilt_current,G0,w0+E0*x_tilt_current,[],[],[],[],[],options);  
%     [u_tilt,fval] = quadprog(H,2*F'*x_tilt_current,zeros(size(G0)),zeros(size(w0+E0*x_tilt_current)),[],[],[],[],[],options);  
%     u_tilt = -inv(H)*F'*x_tilt_current;
    
    u_tilt_set(:,k) = u_tilt;
    x_tilt_set(:,k) = x_tilt_current;     
    
    x_state(:,k) = x_tilt_current + xr(:,k);
    u(:,k) = u_tilt(1:m) + ur(:,k);
    
    x_tilt_current_next = A(:,:,k)*x_tilt_current + B(:,:,k)*u_tilt(1:m);   
    
    dx_state(:,k) = ((x_tilt_current_next+xr(:,k+1)) - x_state(:,k))/T;
    
    x_tilt_current = x_tilt_current_next;
    
end
%zmp        
ddx_state(:,1:length(dx_state)-1) = (dx_state(:,2:end) - dx_state(:,1:end-1))/T;    %forward derivative
ddx_state(:,length(dx_state)) = 0;
zmp = -h2/g*ddx_state;
zmp_tilt = -h2/g*(ddx_state - ddxr(:,1:length(ddx_state)));

%% plot
%trajectory
figure;
plot(xr(1,:),xr(2,:),'r','LineWidth',2)
hold on
plot(x_state(1,:),x_state(2,:),'--b')
title('trajectory')
legend('reference','real')
grid on

%state & input
figure;
subplot(2,2,1)
plot(t,v_ref,'b','LineWidth',2)
hold on
plot(t(1:end-N),u(1,:),'--k')
plot(t,w_ref,'r','LineWidth',2)
plot(t(1:end-N),u(2,:),'--k')
legend('v ref','v real','w ref','w real')
title('input')
grid on

subplot(2,2,2)
plot(t(1:end-N),dx_state(1,:),'b')
hold on
plot(t(1:end-N),dx_state(2,:),'r')
plot(t(1:end-N),dx_state(3,:),'k')
title('dstate')
legend('x','y','theta')
grid on

subplot(2,2,3)
plot(t(1:end-N),x_tilt_set(1,:),'b')
hold on
plot(t(1:end-N),x_tilt_set(2,:),'r')
plot(t(1:end-N),x_tilt_set(3,:),'k')
title('x tilt')
legend('x','y','theta')
grid on

subplot(2,2,4)
plot(t(1:end-N),u_tilt_set(1,:),'b')
hold on
plot(t(1:end-N),u_tilt_set(2,:),'r')
title('u tilt')
legend('v','w')
grid on

%zmp
figure;
subplot(2,2,1)
% figure;
plot(t(1:end-N),zmp(1,:),'b')
hold on
plot([t(1) t(end-N)], [D/2 D/2],'k')
plot([t(1) t(end-N)], [-D/2 -D/2],'k')
title('zmp x')
grid on

subplot(2,2,2)
plot(t(1:end-N),zmp(2,:),'b')
hold on
plot([t(1) t(end-N)], [L/2 L/2],'k')
plot([t(1) t(end-N)], [-L/2 -L/2],'k')
title('zmp y')
grid on

subplot(2,2,3)
plot(t(1:end-N),zmp_tilt(1,:),'b')
hold on
plot(t(1:end), (h2/g*ddxr(1,1:end)+D/2),'k')
plot(t(1:end), (h2/g*ddxr(1,1:end)-D/2),'k')
title('zmp tilt x')
grid on

subplot(2,2,4)
plot(t(1:end-N),zmp_tilt(2,:),'b')
hold on
plot(t(1:end), (h2/g*ddxr(2,1:end)+L/2),'k')
plot(t(1:end), (h2/g*ddxr(2,1:end)-L/2),'k')
title('zmp tilt y')
grid on