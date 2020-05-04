%% Batch formulation of system for QP algorithm
function [Sx,Su,Qb,Rb] = Batch_formulation(A, B, N, Q, R, cur_state)

%% parameter size
n = size(A,1);
m = size(B,2);

%% Sx{n*N,n}
Sx = A(:,:,cur_state);
for i = 1 : N-1    
    Sx = [Sx ; mA(A,cur_state,cur_state+i)];    
end

%% Su{n*N,m*N}
Su = zeros(n*N,m*N);
for i = 1 : N
    for j = 1 : N
        if i-j < 0 %upper than diagonal 
            At = zeros(n,n);
        elseif i-j == 0 %diagonal
            At = eye(n,n);
        else %i-j>0
            At = mA(A,cur_state+1,cur_state+(i-j));
        end        
        Su((i-1)*n+1 : i*n , (j-1)*m+1 : j*m) = [At * B(:,:,cur_state+(j-1))];
    end
end

%% Qb{n*N,n*N}
Q_dummy = Q;
for i = 2 : N    
    Qb = blkdiag(Q_dummy,Q);    
    Q_dummy = Qb;
end

%% Rb{m*N,m*N}
R_dummy = R;
for i = 2 : N
    Rb = blkdiag(R_dummy,R);
    R_dummy = Rb;
end
