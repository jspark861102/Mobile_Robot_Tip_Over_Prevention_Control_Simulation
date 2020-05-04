clear all
% close all
clc

T = 0.01; dt = T;
%% design parameters
N = 15; % Batch size
Qs = 10;
Rs = 1;

%1:nonlinear model, 2:linear model
model_switch = 2;

%1:on 0:off zmp constraint
switch_zmp = 1;

% initial state
% x0 = [0 0 0]'; %[x y theta]
x0 = [0.1 0.05 pi/12]'; %[x y theta]

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

%% reference trajectory
v_max = 0.3; %0.85;%1;
road_width = 0.6;%2;
%xr=[x_ref,y_ref,theta_ref], ur=[v_ref,w_ref]
[t, xr, ur, ddxr] = ref_trajectory_diff_wv(L, v_max, road_width,dt,T);

%% State space equation
for i = 1 : length(t)
    A(:,:,i) = eye(3) + [0 0 -ur(1,i)*sin(xr(3,i))*T;
                         0 0  ur(1,i)*cos(xr(3,i))*T;
                         0 0  0];
    B(:,:,i) = [cos(xr(3,i))*T 0;
                sin(xr(3,i))*T 0;
                0              T];    
end
% linearization_test;

%% MPC parameters
n = size(A,1);
m = size(B,2);
Q = Qs*eye(n);
R = Rs*eye(m);

%% Optimization
options = optimoptions('quadprog',...
     'Algorithm','interior-point-convex');

x0_tilt = x0 - xr(:,1);
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
    [G0,E0,w0] = QP_constraints_ZMP(cur_state, xr, ddxr, x_tilt_current, x_tilt_m1, N, A, B, h2, D, L, T, switch_zmp, switch_input, thr_input, switch_state, thr_state);
    
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
data_plot_diff