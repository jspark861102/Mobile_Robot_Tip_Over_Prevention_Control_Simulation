clear all
close all
clc

T = 0.01; dt = T;
%% design parameters
N = 15; % Batch size
Qs = 10;
Rs = 1;

%1:nonlinear model, 2:linear model
model_switch = 1;

%1:on 0:off zmp constraint
switch_zmp = 1;

%1:on 0:off input constraint
switch_input = 0;
thr_input_plus = [0.1; 0.1]; 
thr_input_minus = [-0.1;-0.1 ];
thr_input = [thr_input_plus; -thr_input_minus];

%1:on 0:off state constraint
switch_state = 0;
thr_state_plus =  [ 0.008;  0.008;  0.008]; 
thr_state_minus = [-0.008; -0.008; -0.008];
thr_state = [thr_state_plus; -thr_state_minus];

%% system parameters
%mass
g = 9.81;
%wheel base
D = 0.3;
%width
L = 0.4;
%z value from global origin to link2 mass center
h2 = 0.9;

%% reference
v_max = 0.3; %0.85;%1;
road_width = 0.6;%2;
circle_r = (road_width - L*cos(pi/4)) / (1-cos(pi/4));

t = 0 : dt : pi *circle_r/(2*abs(v_max)) +0.5 ;

%make input reference
num = 20;
v_ref = v_max*[-cos(pi .* [0:1/num:1])/2+0.5 1*ones(1,length(t) - 2*(num+1)) +cos(pi .* [0:1/num:1])/2+0.5];
w_ref = 1*(v_max/circle_r)*[-cos(pi .* [0:1/num:1])/2+0.5 1*ones(1,length(t) -2*(num+1)) +cos(pi .* [0:1/num:1])/2+0.5];

t = 0 : dt : pi *circle_r/(2*abs(v_max)) +0.5 +3; %to add zeros in front and back

v_ref = [zeros(1,150) v_ref zeros(1,150)];
w_ref = [zeros(1,150) w_ref zeros(1,150)];

%make state reference
theta_ref(1) = 0;
for i = 1 : length(t)-1   
    theta_ref(i+1) = theta_ref(i) + w_ref(i+1) * dt; %backward integration    
end

dx_ref = v_ref.*cos(theta_ref);
dy_ref = v_ref.*sin(theta_ref);
x_ref(1) = 0;
y_ref(1) = 0;
for i = 1 : length(t)-1
    x_ref(i+1) = x_ref(i) + dx_ref(i+1) * dt; %backward integration
    y_ref(i+1) = y_ref(i) + dy_ref(i+1) * dt; %backward integration
end
ddx_ref(1) = 0;
ddy_ref(1) = 0;
ddtheta_ref(1) = 0;% %% input limit
ddx_ref(2:length(dx_ref)) = (dx_ref(2:end)-dx_ref(1:end-1))/T;  %backward derivative
ddy_ref(2:length(dx_ref)) = (dy_ref(2:end)-dy_ref(1:end-1))/T;  %backward derivative
ddtheta_ref(2:length(dx_ref)) = (w_ref(2:end)-w_ref(1:end-1))/T;%backward derivative

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
x0 = [0 0 0]'; %[x y theta]switch_zmp
x0_tilt = x0 - xr(:,1);

%% Optimization
options = optimoptions('quadprog',...
     'Algorithm','interior-point-convex');
 
x_tilt_current = x0_tilt; 
dx_state(:,1) = [0 0 0]';
for k = 1 : length(t)-N    
    %% MPC formulation
    cur_state = k;
    [Sx,Su,Qb,Rb] = Batch_formulation(A, B, N, Q, R, cur_state);

    H = Su'*Qb*Su + Rb;
    F = Sx'*Qb*Su;
    Y = Sx'*Qb*Sx;

    %% constraints
    if k == 1 || k == 2
        x_tilt_m1 = x0_tilt;
    else
        x_tilt_m1 = x_tilt_m1_dummy;
    end        
    [G0,E0,w0] = QP_constraints_ZMP(cur_state, ddxr, x_tilt_m1, N, A, B, h2, D, L, T, switch_zmp, switch_input, thr_input, switch_state, thr_state);
    
    %% QP
    [u_tilt,fval] = quadprog(H,2*F'*x_tilt_current,G0,w0+E0*x_tilt_current,[],[],[],[],[],options);  
    
    %% Simulation
    u_tilt_set(:,k) = u_tilt;
    x_tilt_set(:,k) = x_tilt_current;  
    x_tilt_m1_dummy = x_tilt_current;  
    
    u(:,k) = u_tilt(1:m) + ur(:,k);
    x_state(:,k) = x_tilt_current + xr(:,k);

    if model_switch == 1 %nonlinear model
        if k < length(t)-N  
            dx_state(:,k+1) = [u(1,k)*cos(x_state(3,k)); u(1,k)*sin(x_state(3,k)); u(2,k)];  %backward derivative    
            x_state_next = x_state(:,k) + dx_state(:,k+1)*T;    %backward derivative
            x_tilt_current = x_state_next - xr(:,k+1);
        end
    
    elseif model_switch == 2 %linear model    
        if k < length(t)-N  
            x_tilt_current_next = A(:,:,k)*x_tilt_current + B(:,:,k)*u_tilt(1:m);   
            dx_state(:,k+1) = ((x_tilt_current_next+xr(:,k+1)) - x_state(:,k))/T;  %backward derivative
            x_tilt_current = x_tilt_current_next;  
        end
    end
end

%% zmp        
ddx_state(:,1) = zeros(3,1);
ddx_state(:,2:length(dx_state)) = (dx_state(:,2:end) - dx_state(:,1:end-1))/T;    %backward derivative

zmp = -h2/g*ddx_state;
zmpr = -h2/g*ddxr(:,1:end-N);
zmp_tilt = zmp-zmpr;

%% plot
data_plot