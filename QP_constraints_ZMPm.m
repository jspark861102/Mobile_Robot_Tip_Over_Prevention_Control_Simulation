%% QP constraints, ZMP constraints
function [G0,E0,w0] = QP_constraints_ZMPm(cur_state, A, B, h2, T, N, D, L, xr, ddxr,zmp_switch)
n = size(A,1);
m = size(B,2);

%% ZMP matrix
h = h2*1.0;
g = 9.81;

%% E0{nu*N + nx*N + nf, n}
E01 = zeros(n*N,n);
for i = 1 : N
    E01((i-1)*n+1 : i*n,:) = eval(['mA(A,cur_state,cur_state+' num2str(i) '-1);']);
end
E0 = [-E01;-(-E01)];

%% G0{2*n*N, (m+n)*N}  
G01 = zeros(n*N, (m)*N);
for i = 1 : N 
    eval(['B' num2str(i) '= B_hat(A,B,cur_state,cur_state+' num2str(i) ');']);
    for j = 1 : i     
        G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['B' num2str(i) '(:,:,' num2str(j) ');']);
    end
end
G02 = zeros(n*N, (n)*N);
for i = 1 : N    
    G02((i-1)*n+1 : i*n, (i-1)*n+1 : i*n) =  [-h/g 0 0;0 -h/g 0;0 0 0];
end
G0 = [G01 G02;-G01 -G02];

%% W0{2*n*N,1}
w0 = zeros(2*n*N,1);
w01 = zeros(n*N,1);
w02 = zeros(n*N,1);
for i = 1 : N
    w01((i-1)*n+1 : i*n,:) = [D/2;L/2;0] -[[xr(1,cur_state+i);xr(2,cur_state+i);0]-h/g*[ddxr(1,cur_state+i); ddxr(2,cur_state+i); 0]];
end
for i = 1 : N
    w02((i-1)*n+1 : i*n,:) = [D/2;L/2;0] +[[xr(1,cur_state+i);xr(2,cur_state+i);0]-h/g*[ddxr(1,cur_state+i); ddxr(2,cur_state+i); 0]];
end
w0 = [w01;w02];








