%linearization test
x(:,1) = xr(:,1);
xx(:,1) = xr(:,1);
u = [sin(2*pi*t);2*sin(2*pi*t)];
for i = 1 : length(t)
    x(:,i+1) = A(:,:,i)*x(:,i) + B(:,:,i)*u(:,i);
    dxx(:,i) = [u(1,i)*cos(xx(3,i)); u(1,i)*sin(xx(3,i)); u(2,i)];
    xx(:,i+1) = xx(:,i) + dxx(:,i)*T;
end
figure;plot(x(2,:),'b')
hold on
plot(xx(2,:),'--r')