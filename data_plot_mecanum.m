%% plot
figure;
subplot(2,2,1)
plot(xr(1,:),xr(2,:),'r','LineWidth',2)
hold on
plot(x_state(1,:),x_state(2,:),'--b','LineWidth',2)
plot(xr(1,1),xr(2,1),'or','LineWidth',3)
plot(xr(1,end),xr(2,end),'^r','LineWidth',3)
plot(x_state(1,1),x_state(2,1),'o','LineWidth',3)
plot(x_state(1,end),x_state(2,end),'^b','LineWidth',3)
title('trajectory')
legend('reference trajectory','real trajecotry','Location','Best')
grid on
xlabel('x')
ylabel('y')
if ref_traj_switch == 1
    axis([-0.5 1.0 0 2.5])
elseif ref_traj_switch == 2
    axis([-1.2 1.2 -1.2 1.5])
end
set(gca,'fontsize', 16);

% figure;
subplot(2,2,3)
plot(t(1:end-N),x_tilt_set(1,:),'b','LineWidth',2)
hold on
plot(t(1:end-N),x_tilt_set(2,:),'r','LineWidth',2)
plot(t(1:end-N),x_tilt_set(3,:),'k','LineWidth',2)
if switch_state == 1
    plot([t(1) t(end-N)], [thr_state_plus(1) thr_state_plus(1)],'--b')
    plot([t(1) t(end-N)], [thr_state_plus(2) thr_state_plus(2)],'--r')
    plot([t(1) t(end-N)], [thr_state_plus(3) thr_state_plus(3)],'--k')
    plot([t(1) t(end-N)], [thr_state_plus(1) thr_state_plus(1)],'--b')
    plot([t(1) t(end-N)], [thr_state_plus(2) thr_state_plus(2)],'--r')    
    plot([t(1) t(end-N)], [thr_state_plus(3) thr_state_plus(3)],'--k')
end
title('x tilt')
legend('x','y','theta','Location','Best')
grid on
xlabel('t(sec)')
ylabel('x tilt')
set(gca,'fontsize', 16);

% figure;
subplot(2,2,2)
plot(t,ur(1,:),'b','LineWidth',2)
hold on
plot(t(1:end-N),u(1,:),'--b','LineWidth',2)
plot(t,ur(2,:),'r','LineWidth',2)
plot(t(1:end-N),u(2,:),'--r','LineWidth',2)
plot(t,ur(3,:),'k','LineWidth',2)
plot(t(1:end-N),u(3,:),'--k','LineWidth',2)
legend('vx ref','vx real','vy ref', 'vy real', 'w ref','w real','Location','Best')
title('input')
grid on
xlabel('t(sec)')
ylabel('u')
set(gca,'fontsize', 16);

% figure;
subplot(2,2,4)
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
legend('vx', 'vy', 'w','Location','Best')
grid on
xlabel('t(sec)')
ylabel('u tilt')
set(gca,'fontsize', 16);

figure;
plot(t(1:end-N),ddx_state(1,1:end),'b','LineWidth',2)
hold on
plot(t(1:end-N),ddx_state(2,1:end),'r','LineWidth',2)
plot(t(1:end-N),(abs(cos(x_state(3,1:end)))* D/2 + abs(sin(x_state(3,1:end)))* L/2)*g/h2,'--b','LineWidth',2)
plot(t(1:end-N),(abs(sin(x_state(3,1:end)))* D/2 + abs(cos(x_state(3,1:end)))* L/2)*g/h2,'--r','LineWidth',2)
plot(t(1:end-N),mu(1)*g *ones(1,length(t)-N), ':k','LineWidth',2)
plot(t(1:end-N),(abs(cos(x_state(3,1:end)))*-D/2 + abs(sin(x_state(3,1:end)))*-L/2)*g/h2,'--b','LineWidth',2)
plot(t(1:end-N),(abs(sin(x_state(3,1:end)))*-D/2 + abs(cos(x_state(3,1:end)))*-L/2)*g/h2,'--r','LineWidth',2)
plot(t(1:end-N),-mu(1)*g *ones(1,length(t)-N), ':k','LineWidth',2)
plot(t(1:end-N),mu(2)*g *ones(1,length(t)-N), ':k','LineWidth',2)
plot(t(1:end-N),-mu(2)*g *ones(1,length(t)-N), ':k','LineWidth',2)
title('constraint')
legend('ddx','ddy','zmpx constraint','zmpy constraint','slip constraint','Location','Best')
grid on
xlabel('t(sec)')
ylabel('zmp')
set(gca,'fontsize', 16);
