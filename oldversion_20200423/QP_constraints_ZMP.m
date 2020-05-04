%% QP constraints, ZMP constraints
function [G0,E0,w0] = QP_constraints_ZMP(cur_state, xr, ddxr, x_tilt_current, x_tilt_m1, N, A, B, h2, D, L, T, switch_zmp, switch_input, thr_input, switch_state, thr_state)
n = size(A,1);
m = size(B,2);

%% ZMP matrix
h = h2*1.0;
g = 9.81;

%ZMP = x-h/g*ddx --> -h/g*ddx
Z = [(2*h)/(g*T^2) 0               0;
     0             (2*h)/(g*T^2)   0;
     0             0               0];
Zm = [-h/(g*T^2)  0           0;
      0          -h/(g*T^2)   0;
      0           0           0];
  
G = [h/(g*T^2) 0         0;
     0         h/(g*T^2) 0;
     0         0         0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% backward finite difference for zmp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if switch_zmp == 1
    %% G0{2*n*N, m*N}  
    for b_hat_num = 1 : N    
        %in writing Bi(k+j), in code Bj(k+i), 
        eval(['B' num2str(b_hat_num) '= B_hat(A,B,cur_state,cur_state+' num2str(b_hat_num) ');']);
    end    

    G01 = zeros(n*N, m*N);
    for i = 1 : N 
        for j = 1 : N        
            if i-j == 0 %diagonal
                G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['Zm*B' num2str(i) '(:,:,' num2str(j) ');']);
            elseif i-j == 1 %one lower than diagonal 
                G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['Zm*B' num2str(i) '(:,:,' num2str(j) ') - G*(-2*B' num2str(i-1) '(:,:,' num2str(j) '));']);
            elseif i-j < 0 %upper diagonal is 0
                G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = 0;
            else %lower diagonal
                G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['Zm*B' num2str(i) '(:,:,' num2str(j) ') - G*(-2*B' num2str(i-1) '(:,:,' num2str(j) ') + B' num2str(i-2) '(:,:,' num2str(j) '));']);                            
            end
        end
    end
    G0 = [G01;-G01];

    %% E0{2*n*N,n}
    E01 = zeros(n*N,n);
    for i = 1 : N
        if i == 1
            E01((i-1)*n+1 : i*n,:) = eval(['Zm*A_hat(A,cur_state,cur_state+' num2str(i) ') - G*(-2*eye(n));']);
        elseif i == 2
            E01((i-1)*n+1 : i*n,:) = eval(['Zm*A_hat(A,cur_state,cur_state+' num2str(i) ') - G*(-2*A_hat(A,cur_state,cur_state+' num2str(i-1) ') + eye(n));']);
        else
            E01((i-1)*n+1 : i*n,:) = eval(['Zm*A_hat(A,cur_state,cur_state+' num2str(i) ') - G*(-2*A_hat(A,cur_state,cur_state+' num2str(i-1) ') + A_hat(A,cur_state,cur_state+' num2str(i-2) '));']);
        end
    end    
    E0 = [-E01;-(-E01)];

    %% W0{2*n*N,1}
    w0 = zeros(2*n*N,1);
    w01 = zeros(n*N,1);
    w02 = zeros(n*N,1);
%     for i = 1 : N
%         w01((i-1)*n+1 : i*n,:) = [[D/2;L/2;0] + h/g*[ddxr(1,cur_state+i); ddxr(2,cur_state+i); 0]];
%     end
%     for i = 1 : N
%         w02((i-1)*n+1 : i*n,:) = [[D/2;L/2;0] - h/g*[ddxr(1,cur_state+i); ddxr(2,cur_state+i); 0]];
%     end
    for i = 1 : N
        w01((i-1)*n+1 : i*n,:) = [ abs(cos(x_tilt_current(3)+xr(3,cur_state))) abs(sin(x_tilt_current(3)+xr(3,cur_state))) 0; 
                                  abs(sin(x_tilt_current(3)+xr(3,cur_state))) abs(cos(x_tilt_current(3)+xr(3,cur_state))) 0;
                                   0                      0                      1] * [D/2;L/2;0] + h/g*[ddxr(1,cur_state+i); ddxr(2,cur_state+i); 0];
    end
    for i = 1 : N
        w02((i-1)*n+1 : i*n,:) = [ abs(cos(x_tilt_current(3)+xr(3,cur_state))) abs(sin(x_tilt_current(3)+xr(3,cur_state))) 0; 
                                  abs(sin(x_tilt_current(3)+xr(3,cur_state))) abs(cos(x_tilt_current(3)+xr(3,cur_state))) 0;
                                   0                      0                      1] * [D/2;L/2;0] - h/g*[ddxr(1,cur_state+i); ddxr(2,cur_state+i); 0];
    end
       
    w01 = w01 +[G*x_tilt_m1; zeros(n*(N-1),1)];
    w02 = w02 -[G*x_tilt_m1; zeros(n*(N-1),1)];
    w0 = [w01;w02];
else
    G0 = [];
    E0 = [];
    w0 = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% input limit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if switch_input == 1
    %% input limit
    Au = [eye(m);-eye(m)];
    G1_dummy = Au;
    for i = 2 : N    
        G1 = blkdiag(G1_dummy,Au);
        G1_dummy = G1;
    end    
    E1 = zeros(2*m*N,n);
    w1 = thr_input;
    for i = 2 : N
        w1 = [w1 ;thr_input];
    end
else
    G1 = [];
    E1 = [];
    w1 = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% state limit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if switch_state == 1
    %% state limit
    E22 = zeros(n*N,n);
    for i = 1 : N
        E22((i-1)*n+1 : i*n,:) = eval(['mA(A,cur_state,cur_state+' num2str(i) '-1);']);
    end
    E2 = [-E22;-(-E22)];

    G22 = zeros(n*N, (m)*N);
    for i = 1 : N 
        eval(['B' num2str(i) '= B_hat(A,B,cur_state,cur_state+' num2str(i) ');']);
        for j = 1 : i     
            G22((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['B' num2str(i) '(:,:,' num2str(j) ');']);
        end
    end
    G2 = [G22;-G22];

    w2 = zeros(2*n*N,1);
    w21 = zeros(n*N,1);
    w22 = zeros(n*N,1);
    for i = 1 : N
        w21((i-1)*n+1 : i*n,:) =  thr_state(1:3);
    end
    for i = 1 : N
        w22((i-1)*n+1 : i*n,:) =  thr_state(4:6);
    end
    w2 = [w21;w22];
else
    G2 = [];
    E2 = [];
    w2 = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% constraint completion
G0 = [G0; G1; G2];
E0 = [E0; E1; E2];
w0 = [w0; w1; w2];

if switch_zmp == 0 && switch_input == 0 && switch_state == 0
    G0 = zeros(1,m*N);
    E0 = zeros(1,n);
    w0 = zeros(1,1);
end

