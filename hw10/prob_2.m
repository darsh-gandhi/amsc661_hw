clear; close all;

%% baseline params/ICs

syms x
u0 = piecewise(x<0, 0, (0<x) & (x<1), 2, (1<x) & (x<2), 1, x>2, 0);
x_vals = linspace(-1,7.5,1001);
N=length(x_vals);
dx = x_vals(2)-x_vals(1);

% figure(1) % initial condition
% fplot(u0,[x_vals(1) x_vals(end)],'LineWidth',1.5)
% ylim([-0.5 2.5])

%% lax fried

u = zeros(1,N);
for i=1:N
    if x_vals(i)>0 && x_vals(i)<=1
        u(i) = 2;
    elseif x_vals(i)>1 && x_vals(i)<=2
        u(i) = 1;
    end
end

figure;
plot(x_vals,u,'LineWidth',1.5,'DisplayName','t=0')
ylim([-0.5 2.5])
xlim([-1 7.5])
legend('t=0', 'Location','northeast')
title('Lax-Friedrichs')
hold on
plot_times = [0.5, 1.5, 2.5, 3.5, 5];
curr_plot=1;

dt_fried = dx/2;
u_new = u;
t=0;
while t <= 5
    for j=2:N-1
        u_new(j) = 0.5*(u(j-1) + u(j+1)) - (dt_fried/(4*dx))*(u(j+1)^2 - u(j-1)^2);
    end

    u = u_new;
    t = t+dt_fried;

    if curr_plot <= length(plot_times) && t >= plot_times(curr_plot)
        plot(x_vals,u,'LineWidth',1.5,'DisplayName',['t = ', num2str(plot_times(curr_plot))])
        curr_plot = curr_plot+1;
    end
end

%% richtmyer

dt_richt = 0.4*dx;
u = zeros(1,N);
for i=1:N
    if x_vals(i)>0 && x_vals(i)<=1
        u(i) = 2;
    elseif x_vals(i)>1 && x_vals(i)<=2
        u(i) = 1;
    end
end

figure;
plot(x_vals,u,'LineWidth',1.5,'DisplayName','t=0')
ylim([-0.5 2.5])
xlim([-1 7.5])
legend('t=0', 'Location','northeast')
title('Richtmyer')
hold on
plot_times = [0.5, 1.5, 2.5, 3.5, 5];
curr_plot=1;

u_new = u;
u_half = zeros(1,N-1);
t=0;

while t<=5
    for j=1:N-1
        u_half(j) = 0.5*(u(j) + u(j+1)) - (dt_richt/(2*dx))*0.5*(u(j+1)^2 - u(j)^2);
    end
    
    for j=2:N-1
        u_new(j) = u(j) - (dt_richt/dx)*0.5*(u_half(j)^2 - u_half(j-1)^2);
    end

    % u_new(1)=0;
    % u_new(N)=0;
    u=u_new;
    t=t+dt_richt;

    if curr_plot <= length(plot_times) && t >= plot_times(curr_plot)
        plot(x_vals,u,'LineWidth',1.5,'DisplayName',['t = ', num2str(plot_times(curr_plot))])
        curr_plot = curr_plot+1;
    end
end

%% maccormack

dt_mac = 0.4*dx;
u = zeros(1,N);
for i=1:N
    if x_vals(i)>0 && x_vals(i)<=1
        u(i) = 2;
    elseif x_vals(i)>1 && x_vals(i)<=2
        u(i) = 1;
    end
end

figure;
plot(x_vals,u,'LineWidth',1.5,'DisplayName','t=0')
% ylim([-0.5 2.5])
xlim([-1 7.5])
legend('t=0', 'Location','northeast')
title('MacCormack')
hold on
plot_times = [0.5, 1.5, 2.5, 3.5, 5];
curr_plot=1;

u_new = u;
u_star = zeros(1,N-1);
t=0;

while t<=5
    for j=1:N-1
        u_star(j) = u(j) - (dt_mac/dx)*0.5*(u(j+1)^2 - u(j)^2);
    end

    for j=2:N-1
        u_new(j) = 0.5*(u(j) + u_star(j)) - (dt_mac/(2*dx))*0.5*(u_star(j)^2 - u_star(j-1)^2);
    end

    u=u_new;
    t=t+dt_mac;

    if curr_plot <= length(plot_times) && t >= plot_times(curr_plot)
        plot(x_vals,u,'LineWidth',1.5,'DisplayName',['t = ', num2str(plot_times(curr_plot))])
        curr_plot = curr_plot+1;
    end
end

%% exact

function u_exact = exact_solution(x_vals, t)
    u_exact = zeros(size(x_vals));
    
    if t == 0
        for i = 1:length(x_vals)
            x = x_vals(i);
            if x > 0 && x <= 1
                u_exact(i) = 2;
            elseif x > 1 && x <= 2
                u_exact(i) = 1;
            end
        end

        return
    end
    
    if t < 1
        shock1 = 1 + 1.5*t; % left shock
        shock2 = 2 + 0.5*t; % right shock
        for i = 1:length(x_vals)
            x = x_vals(i);
            if x <= 0
                u_exact(i) = 0;
            elseif x <= 2*t 
                u_exact(i) = x/t;
            elseif x < shock1
                u_exact(i) = 2;
            elseif x < shock2
                u_exact(i) = 1;
            else
                u_exact(i) = 0;
            end
        end
    else
        shock = 2.5*sqrt(t);
        for i = 1:length(x_vals)
            x = x_vals(i);
            if x <= 0
                u_exact(i) = 0;
            elseif x <= shock
                u_exact(i) = x/t;
            else
                u_exact(i) = 0;
            end
        end
    end
end

figure;
plot_times_exact = [0, 0.5, 1.5, 2.5, 3.5, 5];
for k = 1:length(plot_times_exact)
    u_ex = exact_solution(x_vals, plot_times_exact(k));
    plot(x_vals, u_ex, 'LineWidth', 1.5, 'DisplayName', ['t = ' num2str(plot_times_exact(k))])
    hold on
end
ylim([-0.5 2.5])
xlim([-1 7.5])
legend('Location','northeast')
title('Exact')