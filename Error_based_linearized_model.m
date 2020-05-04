
% for i = 1 : length(t)
%     A(:,:,i) = eye(3) + [ 0          ur(3,i)*T  -ur(2,i)*T;
%                          -ur(3,i)*T  0           ur(1,i)*T;
%                           0          0           0        ];
%     B(:,:,i) = eye(3,3)*T;
% end
%% Optimization
options = optimoptions('quadprog',...
     'Algorithm','interior-point-convex');
 
T0_e(:,:,1) = [ cos(x0(3,1)) sin(x0(3,1)) 0;
               -sin(x0(3,1)) cos(x0(3,1)) 0;
                0            0            1];
          
x0_e =  T0_e(:,:,1)*(xr(:,1) - x0); 
x_e_current = x0_e; 
T_e=T0_e;
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
        x_e_m1 = x0_e;
    else
        x_e_m1 = x_e_m1_dummy;
    end        
    [G0,E0,w0] = QP_constraints_ZMP(cur_state, ddxr, x_e_m1, N, A, B, h2, D, L, T, switch_zmp, switch_input, thr_input, switch_state, thr_state);
    
    %% QP
    [u_e,fval] = quadprog(H,2*F'*x_e_current,G0,w0+E0*x_e_current,[],[],[],[],[],options);  
    
    %% Simulation
    u_e_set(:,k) = u_e;
    x_e_set(:,k) = x_e_current;  
    x_e_m1_dummy = x_e_current;     
       
    x_state(:,k) = xr(:,k) - inv(T_e(:,:,k)) * x_e_current;    
    
    u(:,k) = [ur(1,k)*cos(xr(3,k)-x_state(3,k));
              ur(1,k)*cos(xr(3,k)-x_state(3,k));
              ur(3,k)]                           -   u_e(1:m) ;

    if model_switch == 1 %nonlinear model
        if k < length(t)-N  
            dx_state(:,k+1) = [u(1,k)*cos(x_state(3,k)) - u(2,k)*sin(x_state(3,k)); u(1,k)*sin(x_state(3,k)) + u(2,k)*cos(x_state(3,k)); u(3,k)];  %backward derivative    
            x_state_next = x_state(:,k) + dx_state(:,k+1)*T;    %backward derivative
            x_e_current = x_state_next - xr(:,k+1);
        end
    
    elseif model_switch == 2 %linear model    
        if k < length(t)-N  
            x_e_current_next = A(:,:,k)*x_e_current + B(:,:,k)*u_e(1:m);   
            
            T_e(:,:,k+1) = [ cos(xr(3,k+1)-x_e_current_next(3)) sin(xr(3,k+1)-x_e_current_next(3)) 0;
                            -sin(xr(3,k+1)-x_e_current_next(3)) cos(xr(3,k+1)-x_e_current_next(3)) 0;
                             0                                  0                                  1];
            
            dx_state(:,k+1) = ((xr(:,k+1) - inv(T_e(:,:,k+1)) * x_e_current_next) - x_state(:,k))/T;  %backward derivative
            x_e_current = x_e_current_next;  
        end
    end
    
   
    
end
