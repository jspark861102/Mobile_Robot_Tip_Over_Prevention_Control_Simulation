function [t, xr, ur, ddxr] = ref_trajectory_mecanum_xy(T)

%% reference trajectory 
tb = [0:T:10];
dx_ref = -1.0*2*pi/10*cos(2*pi*tb/10);
dy_ref = -1.0*4*pi/10*cos(4*pi*tb/10);

tz=3;
t = [0:T:10+tz];
num=19;
dx_ref = [zeros(1,(tz*100-2*(num+1)/2)/2) 1.0*2*pi/10*(cos(2*pi.*[0:1/num:0.5])/2-0.5) dx_ref 1.0*2*pi/10*(-cos(2*pi.*[0:1/num:0.5])/2-0.5) zeros(1,(tz*100-2*(num+1)/2)/2)];
dy_ref = [zeros(1,(tz*100-2*(num+1)/2)/2) 1.0*4*pi/10*(cos(2*pi.*[0:1/num:0.5])/2-0.5) dy_ref 1.0*4*pi/10*(-cos(2*pi.*[0:1/num:0.5])/2-0.5) zeros(1,(tz*100-2*(num+1)/2)/2)];

theta_ref = [atan2(dy_ref(142),dx_ref(142))*ones(1,141)  atan2(dy_ref(142:end-141),dx_ref(142:end-141))  atan2(dy_ref(1161),dx_ref(1161))*ones(1,141)];
% theta_ref = zeros(1,length(dx_ref));
% theta_ref = atan2(dy_ref,dx_ref);
theta_ref = unwrap(theta_ref);


x_ref(1) = +0.0;
y_ref(1) = +0.0;
for i = 1 : length(t)-1
    x_ref(i+1) = x_ref(i) + dx_ref(i+1) * T; %backward integration
    y_ref(i+1) = y_ref(i) + dy_ref(i+1) * T; %backward integration
end

ddx_ref(1) = 0; ddy_ref(1) = 0;
ddx_ref(2:length(dx_ref)) = (dx_ref(2:end)-dx_ref(1:end-1))/T;  %backward derivative
ddy_ref(2:length(dy_ref)) = (dy_ref(2:end)-dy_ref(1:end-1))/T;  %backward derivative

dtheta_ref(1) = 0;
dtheta_ref(2:length(theta_ref)) = (theta_ref(2:end)-theta_ref(1:end-1))/T;  %backward derivative
% dtheta_ref = (dx_ref.*ddy_ref - dy_ref.*ddx_ref) ./ (dx_ref.^2 +  dy_ref.^2);

ddtheta_ref(1) = 0;
ddtheta_ref(2:length(dtheta_ref)) = (dtheta_ref(2:end)-dtheta_ref(1:end-1))/T;  %backward derivative

%define state and input reference
xr = [x_ref; y_ref; theta_ref];
ur = [dx_ref; dy_ref; dtheta_ref];

ddxr = [ddx_ref; ddy_ref; ddtheta_ref];

% figure;
% plot(x_ref,y_ref)
% 
% figure;
% plot(t,x_ref,t,dx_ref,t,ddx_ref)
% 
% figure;
% plot(t,y_ref,t,dy_ref,t,ddy_ref)
% 
% figure;
% plot(t,theta_ref,t,dtheta_ref,t,ddtheta_ref)