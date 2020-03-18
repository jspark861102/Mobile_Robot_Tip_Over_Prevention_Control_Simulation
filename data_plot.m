
%% plot
%trajectory
figure;
plot(xr(1,:),xr(2,:),'r','LineWidth',2)
hold on
plot(x_state(1,:),x_state(2,:),'--b')
title('trajectory')
legend('reference','real')
grid on

%state & input
figure;
subplot(2,2,1)
plot(t,v_ref,'b','LineWidth',2)
hold on
plot(t(1:end-N),u(1,:),'--k')
plot(t,w_ref,'r','LineWidth',2)
plot(t(1:end-N),u(2,:),'--k')
legend('v ref','v real','w ref','w real')
title('input')
grid on

subplot(2,2,2)
plot(t(1:end-N),dx_state(1,1:end),'b')
hold on
plot(t(1:end-N),dx_state(2,1:end),'r')
plot(t(1:end-N),dx_state(3,1:end),'k')
title('dstate')
legend('x','y','theta')
grid on

subplot(2,2,3)
plot(t(1:end-N),x_tilt_set(1,:),'b')
hold on
plot(t(1:end-N),x_tilt_set(2,:),'r')
plot(t(1:end-N),x_tilt_set(3,:),'k')
title('x tilt')
legend('x','y','theta')
grid on

subplot(2,2,4)
plot(t(1:end-N),u_tilt_set(1,:),'b')
hold on
plot(t(1:end-N),u_tilt_set(2,:),'r')
title('u tilt')
legend('v','w')
grid on

%zmp
figure;
subplot(2,2,1)
% figure;
plot(t(1:end-N),zmp(1,1:end),'b')
hold on
plot([t(1) t(end-N)], [D/2 D/2],'k')
plot([t(1) t(end-N)], [-D/2 -D/2],'k')
title('zmp x')
grid on

subplot(2,2,2)
plot(t(1:end-N),zmp(2,1:end),'b')
hold on
plot([t(1) t(end-N)], [L/2 L/2],'k')
plot([t(1) t(end-N)], [-L/2 -L/2],'k')
title('zmp y')
grid on

subplot(2,2,3)
plot(t(1:end-N),zmp_tilt(1,1:end),'b')
hold on
plot(t(1:end), (h2/g*ddxr(1,1:end)+D/2),'k')
plot(t(1:end), (h2/g*ddxr(1,1:end)-D/2),'k')
title('zmp tilt x')
grid on

subplot(2,2,4)
plot(t(1:end-N),zmp_tilt(2,1:end),'b')
hold on
plot(t(1:end), (h2/g*ddxr(2,1:end)+L/2),'k')
plot(t(1:end), (h2/g*ddxr(2,1:end)-L/2),'k')
title('zmp tilt y')
grid on