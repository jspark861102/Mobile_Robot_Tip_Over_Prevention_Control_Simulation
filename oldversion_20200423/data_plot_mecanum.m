%% plot
figure(1);
subplot(2,3,1)
plot(xr(1,:),xr(2,:),'r','LineWidth',2)
hold on
plot(x_state(1,:),x_state(2,:),'--b','LineWidth',2)
title('trajectory')
legend('reference','real','Location','Best')
grid on
xlabel('x')
ylabel('y')
if ref_traj_switch == 1
    axis([-0.5 1.0 0 2.5])
elseif ref_traj_switch == 2
    axis([-0.9 0.7 -0.9 0.7])
end
set(gca,'fontsize', 16);

subplot(2,3,4)
plot(t(1:end-N),x_tilt_set(1,:),'b','LineWidth',2)
hold on
plot(t(1:end-N),x_tilt_set(2,:),'r','LineWidth',2)
plot(t(1:end-N),x_tilt_set(3,:),'k','LineWidth',2)
if switch_state == 1
    plot([t(1) t(end-N)], [thr_state_plus(1) thr_state_plus(1)],'--b')
    plot([t(1) t(end-N)], [thr_state_plus(1) thr_state_plus(1)],'--b')
    plot([t(1) t(end-N)], [thr_state_plus(2) thr_state_plus(2)],'--r')
    plot([t(1) t(end-N)], [thr_state_plus(2) thr_state_plus(2)],'--r')
    plot([t(1) t(end-N)], [thr_state_plus(3) thr_state_plus(3)],'--k')
    plot([t(1) t(end-N)], [thr_state_plus(3) thr_state_plus(3)],'--k')
end
title('x tilt')
legend('x','y','theta','Location','Best')
grid on
xlabel('t(sec)')
ylabel('x tilt')
set(gca,'fontsize', 16);

subplot(2,3,2)
plot(t,ur(1,:),'k','LineWidth',2)
hold on
plot(t(1:end-N),u(1,:),'--b','LineWidth',2)
plot(t,ur(2,:),'g','LineWidth',2)
plot(t(1:end-N),u(2,:),'--c','LineWidth',2)
plot(t,ur(3,:),'m','LineWidth',2)
plot(t(1:end-N),u(3,:),'--r','LineWidth',2)
legend('v ref','v real','vy ref', 'vy real', 'w ref','w real','Location','Best')
title('input')
grid on
xlabel('t(sec)')
ylabel('u')
set(gca,'fontsize', 16);

subplot(2,3,5)
plot(t(1:end-N),u_tilt_set(1,:),'b','LineWidth',2)
hold on
plot(t(1:end-N),u_tilt_set(2,:),'k','LineWidth',2)
plot(t(1:end-N),u_tilt_set(3,:),'r','LineWidth',2)
if switch_input == 1
    plot([t(1) t(end-N)], [thr_input_plus(1) thr_input_plus(1)],'--b','LineWidth',2)
    plot([t(1) t(end-N)], [thr_input_minus(1) thr_input_minus(1)],'--b','LineWidth',2)
    plot([t(1) t(end-N)], [thr_input_plus(2) thr_input_plus(2)],'--r','LineWidth',2)
    plot([t(1) t(end-N)], [thr_input_minus(2) thr_input_minus(2)],'--r','LineWidth',2)
end
title('u tilt')
legend('v', 'vy', 'w','Location','Best')
grid on
xlabel('t(sec)')
ylabel('u tilt')
set(gca,'fontsize', 16);

% subplot(2,3,3)
% plot(t(1:end-N),zmp(1,1:end),'b','LineWidth',2)
% hold on
% plot(t(1:end-N),zmp(2,1:end),'r','LineWidth',2)
% plot([t(1) t(end-N)], [D/2 D/2],'--b','LineWidth',2)
% plot([t(1) t(end-N)], [-D/2 -D/2],'--b','LineWidth',2)
% plot([t(1) t(end-N)], [L/2 L/2],'--r','LineWidth',2)
% plot([t(1) t(end-N)], [-L/2 -L/2],'--r','LineWidth',2)
% title('zmp')
% legend('zmp x','zmp y','Location','Best')
% grid on
% xlabel('t(sec)')
% ylabel('zmp')
% set(gca,'fontsize', 16);

% subplot(2,3,6)
% plot(t(1:end-N),zmp_tilt(1,1:end),'b','LineWidth',2)
% hold on
% plot(t(1:end-N),zmp_tilt(2,1:end),'r','LineWidth',2)
% plot(t(1:end-N), (-zmpr(1,:)+D/2),'--b','LineWidth',2)
% plot(t(1:end-N), (-zmpr(1,:)-D/2),'--b','LineWidth',2)
% plot(t(1:end-N), (-zmpr(2,:)+L/2),'--r','LineWidth',2)
% plot(t(1:end-N), (-zmpr(2,:)-L/2),'--r','LineWidth',2)
% title('zmp tilt')
% legend('zmp tilt x','zmp tilt y','Location','Best')
% grid on
% xlabel('t(sec)')
% ylabel('zmp tilt')
% set(gca,'fontsize', 16);

subplot(2,3,3)
plot(t(1:end-N),zmp(1,1:end),'b','LineWidth',2)
hold on
plot(t(1:end-N),zmp(2,1:end),'r','LineWidth',2)
plot(t(1:end-N),abs(cos(x_state(3,1:end)))* D/2 + abs(sin(x_state(3,1:end)))* L/2,'--b','LineWidth',2)
plot(t(1:end-N),abs(cos(x_state(3,1:end)))*-D/2 + abs(sin(x_state(3,1:end)))*-L/2,'--b','LineWidth',2)
plot(t(1:end-N),abs(sin(x_state(3,1:end)))* D/2 + abs(cos(x_state(3,1:end)))* L/2,'--r','LineWidth',2)
plot(t(1:end-N),abs(sin(x_state(3,1:end)))*-D/2 + abs(cos(x_state(3,1:end)))*-L/2,'--r','LineWidth',2)
title('zmp')
legend('zmp x','zmp y','Location','Best')
grid on
xlabel('t(sec)')
ylabel('zmp')
set(gca,'fontsize', 16);

subplot(2,3,6)
plot(t(1:end-N),zmp_tilt(1,1:end),'b','LineWidth',2)
hold on
plot(t(1:end-N),zmp_tilt(2,1:end),'r','LineWidth',2)
plot(t(1:end-N), -zmpr(1,:)+abs(cos(x_state(3,1:end)))* D/2 + abs(sin(x_state(3,1:end)))* L/2,'--b','LineWidth',2)
plot(t(1:end-N), -zmpr(1,:)+abs(cos(x_state(3,1:end)))*-D/2 + abs(sin(x_state(3,1:end)))*-L/2,'--b','LineWidth',2)
plot(t(1:end-N), -zmpr(2,:)+abs(sin(x_state(3,1:end)))* D/2 + abs(cos(x_state(3,1:end)))* L/2,'--r','LineWidth',2)
plot(t(1:end-N), -zmpr(2,:)+abs(sin(x_state(3,1:end)))*-D/2 + abs(cos(x_state(3,1:end)))*-L/2,'--r','LineWidth',2)
title('zmp tilt')
legend('zmp tilt x','zmp tilt y','Location','Best')
grid on
xlabel('t(sec)')
ylabel('zmp tilt')
set(gca,'fontsize', 16);


figure;
% plot(t(1:end-N),-h2/g*ddv(1,:),'b','LineWidth',2)
% hold on
% plot(t(1:end-N),-h2/g*ddv(2,:),'r','LineWidth',2)
plot(t(1:end-N),-h2/g*du(1,:),'b','LineWidth',2)
hold on
plot(t(1:end-N),-h2/g*du(2,:),'r','LineWidth',2)
plot([t(1) t(end-N)], [D/2 D/2],'--b','LineWidth',2)
plot([t(1) t(end-N)], [-D/2 -D/2],'--b','LineWidth',2)
plot([t(1) t(end-N)], [L/2 L/2],'--r','LineWidth',2)
plot([t(1) t(end-N)], [-L/2 -L/2],'--r','LineWidth',2)
title('zmp')
legend('zmp x','zmp y','Location','Best')
grid on
xlabel('t(sec)')
ylabel('zmp')
set(gca,'fontsize', 16);