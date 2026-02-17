clear; close all

%% q1

function dydt = arenstorf(t,y,z)
%arenstorf := sets up the ivp problem
% input:
%   t <- time
%   y <- states
%   z <- params
% output:
%   dydt <- states

x=y(1);
y2=y(2);
dx=y(3);
dy=y(4);

D1 = ((x+z(1))^2 + y2^2)^(3/2);
D2 = ((x-z(2))^2 + y2^2)^(3/2);

dx_2 = x + 2*dy - z(2)*(x+z(1))/D1 - z(1)*(x-z(2))/D2;
dy_2 = y2 - 2*dx - z(2)*y2/D1 - z(1)*y2/D2;

dydt = [dx; dy; dx_2; dy_2];

end

colors = 1/255*[0 0 255; 160 32 240; 255 0 0; 255 165 0; 0 255 255; 0 0 0];
Tmax = 17.0652165601579625588917206249;
tspan = [0 Tmax];
init = [0.994, 0, 0, -2.0015851063790852240537862224];
mu = 0.012277471;
params = [mu, 1-mu];

opts = odeset('RelTol',1e-12,'AbsTol',1e-12); %epsi = 1e-12

tic;
[~,y] = ode45(@(t,y) arenstorf(t,y,params), tspan, init);
time(1) = toc;
figure(1)
plot(y(:,1),y(:,2),'LineWidth',1.5,'Color',colors(1,:))
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
title('DOPRI5(4)')
legend(['Time= ', num2str(time(1)), ' sec'],'Location','northwest')
set(gca,'FontSize',14)
axis equal;

figure(2)
time = -1.0*ones(1,3);
Tmax=100; tspan=[0 Tmax];

%testing ode45
tic;
[~,y] = ode45(@(t,y) arenstorf(t,y,params), tspan, init);
time(1) = toc;
subplot(1,3,1)
plot(y(:,1),y(:,2),'LineWidth',1.5,'Color',colors(1,:))
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
title('ode45')
legend(['Time= ', num2str(time(1)), ' sec'],'Location','northwest')
set(gca,'FontSize',14)

%testing ode78
tic;
[~,y] = ode78(@(t,y) arenstorf(t,y,params), tspan, init, opts);
time(2) = toc;
subplot(1,3,2)
plot(y(:,1),y(:,2),'LineWidth',1.5,'Color',colors(2,:))
xlabel('$x$','Interpreter','latex')
% ylabel('$y$','Interpreter','latex')
title('ode78')
legend(['Time= ', num2str(time(2)), ' sec'],'Location','northwest')
set(gca,'FontSize',14)

%testing ode89
tic;
[~,y] = ode89(@(t,y) arenstorf(t,y,params), tspan, init, opts);
time(3) = toc;
subplot(1,3,3)
plot(y(:,1),y(:,2),'LineWidth',1.5,'Color',colors(3,:))
xlabel('$x$','Interpreter','latex')
% ylabel('$y$','Interpreter','latex')
title('ode89')
legend(['Time= ', num2str(time(3)), ' sec'],'Location','northwest')
set(gca,'FontSize',14)

%% 

clear; clc

%% q2

colors = 1/255*[0 0 255; 160 32 240; 255 0 0; 255 165 0; 0 255 255; 0 0 0];
real = linspace(-5,5,1001);
im = linspace(-5,5,1001);
[x,y] = meshgrid(real, im);

z = x + 1i*y;

euler = 1 + z;
midpoint = 1 + z + z.^2/2;
kutta3 = 1 + z + z.^2/2 + z.^3/6;
rk4 = 1 + z + z.^2/2 + z.^3/6 + z.^4/24;
dopri5 = 1 + z + z.^2/2 + z.^3/6 + z.^4/24 + z.^5/120 + z.^6/600;

%% q2 plotting

figure(3)
hold on;

contour(x,y,abs(euler), [1 1], 'LineColor', colors(1,:), 'LineWidth', 2.0);
contour(x,y,abs(midpoint), [1 1], 'LineColor', colors(2,:), 'LineWidth', 2.0);
contour(x,y,abs(kutta3), [1 1], 'LineColor', colors(3,:), 'LineWidth', 2.0);
contour(x,y,abs(rk4), [1 1], 'LineColor', colors(4,:), 'LineWidth', 2.0);
contour(x,y,abs(dopri5), [1 1], 'LineColor', colors(5,:), 'LineWidth', 2.0);

xlabel('Re($z$)','Interpreter','latex','FontSize',14);
ylabel('Im($z$)','Interpreter','latex','FontSize',14);
title('Regions of Absolute Stability','FontSize',16);
legend('Forward Euler', 'Midpoint with Euler Predictor','3-stage, 3rd-order RK', 'RK4', 'DOPRI5(4)','Location','east','FontSize',11);
grid on;
box on;
axis equal;
xlim([-4.5 4.5]);
ylim([-4.5 4.5]);