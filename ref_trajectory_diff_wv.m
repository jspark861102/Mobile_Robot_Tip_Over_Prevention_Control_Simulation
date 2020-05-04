function [t, xr, ur, ddxr] = ref_trajectory_diff_wv(L, v_max, road_width,dt,T)

%% reference trajectory
circle_r = (road_width - L*cos(pi/4)) / (1-cos(pi/4));

t = 0 : dt : pi *circle_r/(2*abs(v_max)) +0.5 ;

%make input reference
num = 20;
v_ref = v_max*[-cos(pi .* [0:1/num:1])/2+0.5 1*ones(1,length(t) - 2*(num+1)) +cos(pi .* [0:1/num:1])/2+0.5];
w_ref = 1*(v_max/circle_r)*[-cos(pi .* [0:1/num:1])/2+0.5 1*ones(1,length(t) -2*(num+1)) +cos(pi .* [0:1/num:1])/2+0.5];

t = 0 : dt : pi *circle_r/(2*abs(v_max)) +0.5 +3; %to add zeros in front and back

v_ref = [zeros(1,150) v_ref zeros(1,150)];
w_ref = [zeros(1,150) w_ref zeros(1,150)];

%make state reference
theta_ref(1) = 0;
for i = 1 : length(t)-1   
    theta_ref(i+1) = theta_ref(i) + w_ref(i+1) * dt; %backward integration    
end

dx_ref = v_ref.*cos(theta_ref);
dy_ref = v_ref.*sin(theta_ref);
x_ref(1) = 0;
y_ref(1) = 0;
for i = 1 : length(t)-1
    x_ref(i+1) = x_ref(i) + dx_ref(i+1) * dt; %backward integration
    y_ref(i+1) = y_ref(i) + dy_ref(i+1) * dt; %backward integration
end
ddx_ref(1) = 0;
ddy_ref(1) = 0;
ddtheta_ref(1) = 0;% %% input limit
ddx_ref(2:length(dx_ref)) = (dx_ref(2:end)-dx_ref(1:end-1))/T;  %backward derivative
ddy_ref(2:length(dx_ref)) = (dy_ref(2:end)-dy_ref(1:end-1))/T;  %backward derivative
ddtheta_ref(2:length(dx_ref)) = (w_ref(2:end)-w_ref(1:end-1))/T;%backward derivative

%define state and input reference
xr = [x_ref; y_ref; theta_ref];
ur = [v_ref; w_ref];

ddxr = [ddx_ref; ddy_ref; ddtheta_ref];