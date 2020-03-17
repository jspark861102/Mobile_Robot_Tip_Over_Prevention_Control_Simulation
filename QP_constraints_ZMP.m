%% QP constraints, ZMP constraints
function [G0,E0,w0] = QP_constraints_ZMP(cur_state, A, B, h2, T, N, D, L, ddxr,zmp_switch)
n = size(A,1);
m = size(B,2);

%% ZMP matrix
h = h2;
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

if zmp_switch == 1 %central 
    %%%%%%%%%%%%% i+1 -> i+N-1 central, i+N backward %%%%%%%%%%%%%%
    %% G0{2*n*N, m*N}  
    for b_hat_num = 1 : N    
        %in writing Bi(k+j), in code Bj(k+i), 
        eval(['B' num2str(b_hat_num) '= B_hat(A,B,cur_state,cur_state+' num2str(b_hat_num) ');']);
    end    
    %ZMP with central differentiate method
    G01 = zeros(n*N, m*N);
    for i = 1 : N 
        for j = 1 : N        
            if i == N %backword derivative
                if i-j == 0 %diagonal
                    G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['Zm*B' num2str(i) '(:,:,' num2str(j) ');']);
                elseif i-j == 1 %one lower than diagonal 
                    G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['Zm*B' num2str(i) '(:,:,' num2str(j) ') - G*(-2*B' num2str(i-1) '(:,:,' num2str(j) '));']);
                else %lower diagonal
                    G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['Zm*B' num2str(i) '(:,:,' num2str(j) ') - G*(-2*B' num2str(i-1) '(:,:,' num2str(j) ') + B' num2str(i-2) '(:,:,' num2str(j) '));']);
                end    
            else %central derivative
                if i-j == -1 %one upper than diagonal
                    G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['-G*B' num2str(i+1) '(:,:,' num2str(j) ');']);             
                elseif i-j == 0 %diagonal
                    G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['Z*B' num2str(i) '(:,:,' num2str(j) ') - G*B' num2str(i+1) '(:,:,' num2str(j) ');']);
                elseif i-j < -1 %second upper than diagonal is 0
                    G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = 0;
                else %lower diagonal
                    G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['Z*B' num2str(i) '(:,:,' num2str(j) ') - G*(B' num2str(i+1) '(:,:,' num2str(j) ') + B' num2str(i-1) '(:,:,' num2str(j) '));']);
                end                
            end
        end
    end
    G0 = [G01;-G01];
    % G0 = -G0;

    %% E0{2*n*N,n}
    E01 = zeros(n*N,n);
    for i = 1 : N
        if i == 1
            E01((i-1)*n+1 : i*n,:) = eval(['Z*A_hat(A,cur_state,cur_state+' num2str(i) ') - G*(A_hat(A,cur_state,cur_state+' num2str(i+1) ') + eye(3));']);
        elseif i == N
            E01((i-1)*n+1 : i*n,:) = eval(['Zm*A_hat(A,cur_state,cur_state+' num2str(i) ') - G*(-2*A_hat(A,cur_state,cur_state+' num2str(i-1) ') + A_hat(A,cur_state,cur_state+' num2str(i-2) '));']);
        else
            E01((i-1)*n+1 : i*n,:) = eval(['Z*A_hat(A,cur_state,cur_state+' num2str(i) ') - G*(A_hat(A,cur_state,cur_state+' num2str(i+1) ') + A_hat(A,cur_state,cur_state+' num2str(i-1) '));']);
        end
    end
    E0 = [-E01;-(-E01)];
    % E0=-E0;

elseif zmp_switch == 2 %backward    
    %%%%%%%%%%%% i+1 central, i+2 -> i+N backward %%%%%%%%%%%%%%
    %% G0{2*n*N, m*N}  
    for b_hat_num = 1 : N    
        %in writing Bi(k+j), in code Bj(k+i), 
        eval(['B' num2str(b_hat_num) '= B_hat(A,B,cur_state,cur_state+' num2str(b_hat_num) ');']);
    end    
    %ZMP with central differentiate method
    G01 = zeros(n*N, m*N);
    for i = 1 : N 
        for j = 1 : N        
            if i == 1 %central derivative
                if i-j == -1 %one upper than diagonal
    %                 G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['-G*B' num2str(i+1) '(:,:,' num2str(j) ');']);                             
                    G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['-G*B' num2str(i) '(:,:,' num2str(j-1) ');']);                             
                elseif i-j == 0 %diagonal
    %                 G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['Z*B' num2str(i) '(:,:,' num2str(j) ') - G*B' num2str(i+1) '(:,:,' num2str(j) ');']);                
                    G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['Z*B' num2str(i) '(:,:,' num2str(j) ') - G*B' num2str(i) '(:,:,' num2str(j) ');']);                
                else 
                    G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = 0;            
                end                            
            else %backward derivative            
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
    end
    G0 = [G01;-G01];
    % G0 = -G0;

    %% E0{2*n*N,n}
    E01 = zeros(n*N,n);
    for i = 1 : N
        if i == 1
            E01((i-1)*n+1 : i*n,:) = eval(['Z*A_hat(A,cur_state,cur_state+' num2str(i) ') - G*(A_hat(A,cur_state,cur_state+' num2str(i+1) ') + eye(3));']);
        elseif i == 2
            E01((i-1)*n+1 : i*n,:) = eval(['Zm*A_hat(A,cur_state,cur_state+' num2str(i) ') - G*(-2*A_hat(A,cur_state,cur_state+' num2str(i-1) ') + eye(3));']);
        else
            E01((i-1)*n+1 : i*n,:) = eval(['Zm*A_hat(A,cur_state,cur_state+' num2str(i) ') - G*(-2*A_hat(A,cur_state,cur_state+' num2str(i-1) ') + A_hat(A,cur_state,cur_state+' num2str(i-2) '));']);
        end
    end
    E0 = [-E01;-(-E01)];
    % E0=-E0;
    
elseif zmp_switch == 3 %forward
    %%%%%%%%%%%%% i+1 -> i+N-2 forward, i+N-1 central, i+N backward %%%%%%%%%%%%%%
    %% G0{2*n*N, m*N}  
    for b_hat_num = 1 : N    
        %in writing Bi(k+j), in code Bj(k+i), 
        eval(['B' num2str(b_hat_num) '= B_hat(A,B,cur_state,cur_state+' num2str(b_hat_num) ');']);
    end    
    %ZMP with central differentiate method
    G01 = zeros(n*N, m*N);
    for i = 1 : N 
        for j = 1 : N        
            if i == N %backward
                if i-j == 0 %diagonal
                    G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['Zm*B' num2str(i) '(:,:,' num2str(j) ');']);
                elseif i-j == 1 %one lower than diagonal 
                    G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['Zm*B' num2str(i) '(:,:,' num2str(j) ') - G*(-2*B' num2str(i-1) '(:,:,' num2str(j) '));']);
                else %lower diagonal
                    G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['Zm*B' num2str(i) '(:,:,' num2str(j) ') - G*(-2*B' num2str(i-1) '(:,:,' num2str(j) ') + B' num2str(i-2) '(:,:,' num2str(j) '));']);
                end   
            elseif i == N-1 %central
                if i-j == -1 %one upper than diagonal
                    G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['-G*B' num2str(i+1) '(:,:,' num2str(j) ');']);             
                elseif i-j == 0 %diagonal
                    G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['Z*B' num2str(i) '(:,:,' num2str(j) ') - G*B' num2str(i+1) '(:,:,' num2str(j) ');']);
                else 
                    G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['Z*B' num2str(i) '(:,:,' num2str(j) ') - G*(B' num2str(i+1) '(:,:,' num2str(j) ') + B' num2str(i-1) '(:,:,' num2str(j) '));']);
                end               
            else %forward
                if i-j == -2 %second upper diagonal
                    G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['-G*B' num2str(i+2) '(:,:,' num2str(j) ');']);             
                elseif i-j == -1 %one upper diagonal
                    G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['-G*(B' num2str(i+2) '(:,:,' num2str(j) ')-2*B' num2str(i+1) '(:,:,' num2str(j) '));']);             
                elseif i-j < -2
                    G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = 0;
                else
                    G01((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = eval(['Zm*B' num2str(i) '(:,:,' num2str(j) ') - G*(B' num2str(i+2) '(:,:,' num2str(j) ') -2* B' num2str(i+1) '(:,:,' num2str(j) '));']);
                end                 
            end        
        end
    end
    G0 = [G01;-G01];
    % G0 = -G0;

    %% E0{2*n*N,n}
    E01 = zeros(n*N,n);
    for i = 1 : N
        if i == N
            E01((i-1)*n+1 : i*n,:) = eval(['Zm*A_hat(A,cur_state,cur_state+' num2str(i) ') - G*(-2*A_hat(A,cur_state,cur_state+' num2str(i-1) ') + A_hat(A,cur_state,cur_state+' num2str(i-2) '));']);
        elseif i == N-1
            E01((i-1)*n+1 : i*n,:) = eval(['Z*A_hat(A,cur_state,cur_state+' num2str(i) ') - G*(A_hat(A,cur_state,cur_state+' num2str(i+1) ') + A_hat(A,cur_state,cur_state+' num2str(i-1) '));']);
        else
            E01((i-1)*n+1 : i*n,:) = eval(['Zm*A_hat(A,cur_state,cur_state+' num2str(i) ') - G*(A_hat(A,cur_state,cur_state+' num2str(i+2) ') -2* A_hat(A,cur_state,cur_state+' num2str(i+1) '));']);
        end
    end
    E0 = [-E01;-(-E01)];
    % E0=-E0;
end

%% W0{2*n*N,1}
w0 = zeros(2*n*N,1);
w01 = zeros(n*N,1);
w02 = zeros(n*N,1);
for i = 1 : N
    w01((i-1)*n+1 : i*n,:) = [[D/2;L/2;0]*1.0 + h/g*[ddxr(1,cur_state+i); ddxr(2,cur_state+i); 0]]*1;
%     w01((i-1)*n+1 : i*n,:) = [[D/2;L/2;0] + h/g*[ddxr(1,cur_state+i); ddxr(2,cur_state+i); 0]]*1;
%     w01((i-1)*n+1 : i*n,:) = h/g*[ddxr(1,cur_state+i); ddxr(2,cur_state+i); 0];
end
for i = 1 : N
    w02((i-1)*n+1 : i*n,:) = [[D/2;L/2;0]*1.0 - h/g*[ddxr(1,cur_state+i); ddxr(2,cur_state+i); 0]]*1;
%     w02((i-1)*n+1 : i*n,:) = [[D/2;L/2;0] - h/g*[ddxr(1,cur_state+i); ddxr(2,cur_state+i); 0]]*1;
%     w02((i-1)*n+1 : i*n,:) = -h/g*[ddxr(1,cur_state+i); ddxr(2,cur_state+i); 0];
end
w0 = [w01;w02];

