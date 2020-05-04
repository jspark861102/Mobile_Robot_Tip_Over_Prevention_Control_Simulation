%1. in linear model, ZMP is working but, QP is working only for t=1:860
%1-->. So, the time is reduced and simulation is completed. but don't know why
%2. in nonlinear model, simulation not working why?

clear all
close all
clc

T = 0.01; dt = T;

%% design parameters
N = 15; % Batch size
Qs = 5;
Rs = 1;

%actuator dynamics
tau1=0.01; tau2=0.01; tau3=0.01;
Au = [eye(3) - diag([T/tau1, T/tau2, T/tau3])];
Bu = diag([T/tau1, T/tau2, T/tau3]);

%1:nonlinear model, 2:linear model
model_switch = 2;

%reference trajectory 1:wv, 2:xy
ref_traj_switch = 1;

%1:on 0:off zmp constraint
switch_zmp = 1;

%initial state : [x y theta]
x0 = [0 0 0 0 0 0]'; %actuator dynamics
% if ref_traj_switch == 1
%     x0 = [0.1 0.2 pi/8]'; %wv
% elseif ref_traj_switch == 2
%     x0 = [-0.1 0.1 -2.0344]'; %xy
% end

%1:on 0:off input constraint
switch_input = 0;
mag = 0.3;
thr_input_plus = [mag; mag; mag]; 
thr_input_minus = [-mag; -mag; -mag ];
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

%% reference trajectory xr=[x_ref, y_ref, theta_ref], ur=[dx_ref, dy_ref, dtheta_ref]
if ref_traj_switch == 1
    v_max = 0.3; %0.85;%1;
    road_width = 0.6;%2;
    [t, xr, ur, ddxr] = ref_trajectory_mecanum_wv(L, v_max, road_width,dt,T);    
    
    %actuator dynamics
    xr = [xr;ur];
    for  i = 1 : length(ur)-1
        vr(:,i) = inv(Bu)*(ur(:,i+1)-Au*ur(:,i));
    end
    ur = vr;
    ur(:,end+1) = ur(:,end);
    
    t = t(1:end-200);
    xr = xr(:,1:end-200);
    ur = ur(:,1:end-200);
    ddxr = ddxr(:,1:end-200);    
    
elseif ref_traj_switch == 2
    [t, xr, ur, ddxr] = ref_trajectory_mecanum_xy(T);
%     [t, xr, ur, ddxr] = ref_trajectory_mecanum_xy1(T);
end

%% State space equation
for i = 1 : length(t)    
    An(:,:,i) = eye(3) + [0 0  (-ur(1,i)*sin(xr(3,i)) - ur(2,i)*cos(xr(3,i)))*T;
                         0 0  ( ur(1,i)*cos(xr(3,i)) - ur(2,i)*sin(xr(3,i)))*T;
                         0 0  0];
    Bn(:,:,i) = [cos(xr(3,i))*T -sin(xr(3,i))*T 0;
                sin(xr(3,i))*T  cos(xr(3,i))*T 0;
                0               0              T];    
    
    %actuator dynamics
    A(:,:,i) = [An(:,:,i) Bn(:,:,i); zeros(3,3) Au];
    B(:,:,i) = [zeros(3,3); Bu];
end

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
dx_state(:,1) = [0 0 0 0 0 0]';     %actuator dynamics
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
    [G0,E0,w0] = QP_constraints_ZMP_actuator(cur_state, xr, ddxr, x_tilt_current, x_tilt_m1, N, A, B, h2, D, L, T, switch_zmp, switch_input, thr_input, switch_state, thr_state);
    
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
            dx_state(1:3,k+1) = [x_state(4,k)*cos(x_state(3,k)) - x_state(5,k)*sin(x_state(3,k)); x_state(4,k)*sin(x_state(3,k)) + x_state(5,k)*cos(x_state(3,k)); x_state(6,k)];  %backward derivative    
            dx_state(4:6,k+1) = -diag([1/tau1, 1/tau2, 1/tau3])*x_state(4:6,k) + diag([1/tau1, 1/tau2, 1/tau3])*u(:,k);  
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
%actuator dynamics
ddx_state(:,1) = zeros(6,1);

ddx_state(:,2:length(dx_state)) = (dx_state(:,2:end) - dx_state(:,1:end-1))/T;    %backward derivative

zmp = -h2/g*ddx_state;
zmpr = -h2/g*ddxr(:,1:end-N);
zmp_tilt = zmp(1:3,:)-zmpr;

%% plot
data_plot_mecanum