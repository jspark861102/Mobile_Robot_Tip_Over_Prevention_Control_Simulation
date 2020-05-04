function [t, xr, ur, ddxr] = ref_trajectory_mecanum_xy(T)

%% reference trajectory 
tb = [0:T:10];
x_ref = 0 + 0.7*sin(2*pi*tb/10);
y_ref = 0 + 0.7*sin(4*pi*tb/10);

t = [0:T:12];
x_ref = [zeros(1,100) x_ref zeros(1,100)];
y_ref = [zeros(1,100) y_ref zeros(1,100)];
% theta_ref = [zeros(1,100) y_ref zeros(1,100)];

dx_ref(1) = 0; dy_ref(1) = 0;
dx_ref(2:length(x_ref)) = (x_ref(2:end)-x_ref(1:end-1))/T;  %backward derivative
dy_ref(2:length(y_ref)) = (y_ref(2:end)-y_ref(1:end-1))/T;  %backward derivative
dx_ref(1) = dx_ref(2); dy_ref(1) = dy_ref(2); %first value is not used, but for smooth acceleration, v_ref(1)=v_ref(2)

ddx_ref(1) = 0; ddy_ref(1) = 0;
ddx_ref(2:length(dx_ref)) = (dx_ref(2:end)-dx_ref(1:end-1))/T;  %backward derivative
ddy_ref(2:length(dy_ref)) = (dy_ref(2:end)-dy_ref(1:end-1))/T;  %backward derivative
% ddx_ref(1) = ddx_ref(2); ddy_ref(1) = ddy_ref(2); 

theta_ref = atan2(dy_ref,dx_ref);
dtheta_ref = (dx_ref.*ddy_ref - dy_ref.*ddx_ref) ./ (dx_ref.^2 +  dy_ref.^2);

ddtheta_ref(1) = 0;
ddtheta_ref(2:length(dtheta_ref)) = (dtheta_ref(2:end)-dtheta_ref(1:end-1))/T;  %backward derivative
% ddtheta_ref(1) = ddtheta_ref(2);

%make zeros
% t = tb;
% t = [0:T:12];
% x_ref = [zeros(1,100) x_ref zeros(1,100)];
% y_ref = [zeros(1,100) y_ref zeros(1,100)];
% % theta_ref = [zeros(1,100) y_ref zeros(1,100)];

%define state and input reference
xr = [x_ref; y_ref; theta_ref];
ur = [dx_ref; dy_ref; dtheta_ref];

ddxr = [ddx_ref; ddy_ref; ddtheta_ref];

figure;
plot(t,x_ref,t,dx_ref)%,t,ddx_ref)

figure;
plot(t,y_ref,t,dy_ref,t,ddy_ref)

figure;
plot(t,theta_ref,t,dtheta_ref,t,ddtheta_ref)