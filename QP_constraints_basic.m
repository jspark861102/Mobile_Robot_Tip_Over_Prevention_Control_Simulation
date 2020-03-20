%% QP constraints, ZMP constraints
function [G0,E0,w0] = QP_constraints_basic(cur_state, A, B, h2, T, N, D, L, ddxr)
n = size(A,1);
m = size(B,2);

%% G0{2*m*N, m*N}  
Au = [eye(2);-eye(2)];
thrp = [-0.1; 0.01]; 
thrm = [-0.2;-0.2 ];
thr = [thrp;-thrm];

G0_dummy = Au;
for i = 2 : N    
    G0 = blkdiag(G0_dummy,Au);
    G0_dummy = G0;
end
    
%% E0{2*n*N,n}
E0 = zeros(2*m*N,n);

%% W0{2*n*N,1}
w0 = thr;
for i = 2 : N
    w0 = [w0 ;thr];
end
