clear; close all;

n = 99; % # of int pts
h = 2/(n+1); % step size space
xi = (-1+h:h:1-h)';

%% ode function

function dydt = rhs(t, y, n, h, xi)
    u = y(1:n);
    xf = y(n+1);
    u_full = [0;u;0]; % sol vector with endpts

    % 1st derivative du/d(xi)
    du = zeros(n,1);
    for j = 1:n
        du(j) = (u_full(j+2)-u_full(j))/(2*h);
    end

    % 2nd der d^2u/d(xi)^2
    d2u = zeros(n,1);
    for j = 1:n
        d2u(j) = (u_full(j)-2*u_full(j+1)+u_full(j+2))/h^2;
    end

    du_left = (1/(2*h))*(4*u(1)-u(2)); % du/d(xi) at left boundary
    du_right = (1/(2*h))*(-4*u(n)+u(n-1)); % du/d(xi) at left boundary

    dxi_dt = ((1+xi)*du_right + (1-xi)*du_left); % d(xi)/dt
    dudt = (1/xf^2)*(-0.5*du.*dxi_dt + u.*d2u + du.^2);
    dxf_dt = -(1/(2*xf))*(du_right - du_left);

    dydt = [dudt; dxf_dt];
end

%% initial conds

%IC 1 (u = 1-x^2, |x|<1; 0, otherwise)
u0_1 = @(x) max(1-x.^2, 0);
u0_ic1 = u0_1(xi);
xf0_ic1 = 1.0;
y0_ic1 = [u0_ic1; xf0_ic1];

% IC 2 (u = 1 - 0.99cos(2pi x), |x|<1; 0, otherwise)
u0_func2 = @(x) max(1-0.99*cos(2*pi*x), 0);
u0_ic2 = u0_func2(xi);
xf0_ic2 = 1.0;
y0_ic2 = [u0_ic2; xf0_ic2];

%% numerical sol

tspan = [0,1.2];
t2= 0.1:0.1:1.2;
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);

sol1 = ode15s(@(t,y) rhs(t,y,n,h,xi), tspan, y0_ic1, opts);
sol2 = ode15s(@(t,y) rhs(t,y,n,h,xi), tspan, y0_ic2, opts);

Y1 = deval(sol1, t2);
Y2 = deval(sol2, t2);

xi_full = [-1; xi; 1];

%% plotting IC1 sols

figure(1);
for k = 1:length(t2)
    u_plot = [0; Y1(1:n, k); 0];
    plot(xi_full, u_plot, 'DisplayName', sprintf('t=%.1f', t2(k)),'LineWidth',1.25);
    hold on;
end
xlabel('\xi'); ylabel('u');
title('u(\xi,t) = 1-x^2');
legend show; grid on;

figure(2);
for k = 1:length(t2)
    u_plot = [0; Y1(1:n, k); 0];
    u_max_k = max(u_plot);
    plot(xi_full, u_plot/u_max_k, 'DisplayName', sprintf('t=%.1f', t2(k)),'LineWidth',1.25);
    hold on;
end
plot(xi_full, 1 - xi_full.^2, 'k--', 'LineWidth', 2, 'DisplayName', '1-\xi^2');
xlabel('\xi'); ylabel('u/u_{max}');
ylim([0,1.05])
title('u(\xi,t)/u_{max}(t)');
legend show; grid on;

%% plotting IC2 sols

figure(3);
for k = 1:length(t2)
    u_plot = [0; Y2(1:n, k); 0];
    plot(xi_full, u_plot, 'DisplayName', sprintf('t=%.1f', t2(k)),'LineWidth',1.25);
    hold on;
end
xlabel('\xi'); ylabel('u');
title('u(\xi,t) = 1 - 0.99cos(2\pix)');
legend show; grid on;

figure(4);
for k = 1:length(t2)
    u_plot = [0; Y2(1:n, k); 0];
    u_max_k = max(u_plot);
    plot(xi_full, u_plot/u_max_k, 'DisplayName', sprintf('t=%.1f', t2(k)),'LineWidth',1.25);
    hold on;
end
plot(xi_full, 1 - xi_full.^2,'k--','LineWidth',2, 'DisplayName', '1-\xi^2');
xlabel('\xi'); ylabel('u/u_{max}');
ylim([0,1.05])
title('u(\xi,t)/u_{max}(t)');
legend show; grid on;