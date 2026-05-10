clear; close all

%%

N=1501;
x_vals = linspace(-1,5,N); % x space
dx = x_vals(2)-x_vals(1); % spatial step

u0 = zeros(1,N); %IC
for i=1:N
    if x_vals(i)>=0 && x_vals(i)<=1
        u0(i) = 1;
    end
end

colors = 1/255*[0 0 255 0.8*255; 160 32 240 0.8*255; 255 0 0 0.8*255; 255 165 0 0.8*255; 0 255 255 0.8*255; 0 0 0 0.8*255]; % blue purple red orange cyan black
colors_trans = 1/255*[0 0 255 0.4*255; 160 32 240 0.4*255; 255 0 0 0.4*255; 255 165 0 0.4*255; 0 255 255 0.4*255; 0 0 0 0.4*255]; % blue purple red orange cyan black

%% Exact

plot_times_exact = [1, 2, 3, 4, 5, 6];

function u_exact = exact_solution(x_vals, t)
    u_exact = zeros(size(x_vals));
    
    if t == 0
        for i = 1:length(x_vals)
            x = x_vals(i);
            if x >= 0 && x <= 1
                u_exact(i) = 1;
            end
        end

        return
    end
    
    if t <= 2
        shock = 1 + 0.5*t;
        for i = 1:length(x_vals)
            x = x_vals(i);
            if x < 0
                u_exact(i) = 0;
            elseif x <= t 
                u_exact(i) = x/t;
            elseif x < shock
                u_exact(i) = 1;
            else
                u_exact(i) = 0;
            end
        end
    else
        shock = sqrt(2*t);
        for i = 1:length(x_vals)
            x = x_vals(i);
            if x < 0
                u_exact(i) = 0;
            elseif x <= shock
                u_exact(i) = x/t;
            else
                u_exact(i) = 0;
            end
        end
    end
end

%% lax-fried

u=u0;

figure;
plot(x_vals,u,'LineWidth',1.5,'DisplayName','t=0')
ylim([-0.5 1.5])
xlim([-1 5])
legend('t=0', 'Location','northeast')
title('Lax-Friedrichs')
hold on
plot_times = [1, 2, 3, 4, 5, 6];
curr_plot=1;

dt_fried = dx/2;
u_new = u;
t=0;
while t <= 6
    for j=2:N-1
        u_new(j) = 0.5*(u(j-1) + u(j+1)) - (dt_fried/(4*dx))*(u(j+1)^2 - u(j-1)^2);
    end

    u = u_new;
    t = t+dt_fried;

    if curr_plot <= length(plot_times) && t >= plot_times(curr_plot)
        plot(x_vals,u,'LineWidth',1.5,'Color',colors(curr_plot,:))

        u_ex = exact_solution(x_vals, plot_times_exact(curr_plot));
        plot(x_vals, u_ex, 'LineWidth', 1.5,'Color',colors_trans(curr_plot,:))
        hold on

        curr_plot = curr_plot+1;
    end
end
legend('t=0','t=1','','t=2','','t=3','','t=4','','t=5','','t=6','')

%% Richtmyer

u=u0;

figure;
plot(x_vals,u,'LineWidth',1.5,'DisplayName','t=0')
ylim([-0.5 1.5])
xlim([-1 5])
legend('t=0', 'Location','northeast')
title('Richtmyer')
hold on
plot_times = [1, 2, 3, 4, 5, 6];
curr_plot=1;

dt_richt = 0.4*dx;
u_new = u;
u_half = zeros(1,N-1);
t=0;

while t<=6
    for j=1:N-1
        u_half(j) = 0.5*(u(j) + u(j+1)) - (dt_richt/(2*dx))*0.5*(u(j+1)^2 - u(j)^2);
    end
    
    for j=2:N-1
        u_new(j) = u(j) - (dt_richt/dx)*0.5*(u_half(j)^2 - u_half(j-1)^2);
    end

    u=u_new;
    t=t+dt_richt;

    if curr_plot <= length(plot_times) && t >= plot_times(curr_plot)
        plot(x_vals,u,'LineWidth',1.5,'Color',colors(curr_plot,:))

        u_ex = exact_solution(x_vals, plot_times_exact(curr_plot));
        plot(x_vals, u_ex, 'LineWidth', 1.5,'Color',colors_trans(curr_plot,:))
        hold on

        curr_plot = curr_plot+1;
    end
end
legend('t=0','t=1','','t=2','','t=3','','t=4','','t=5','','t=6','')

%% MacCormack

u=u0;

figure;
plot(x_vals,u,'LineWidth',1.5,'DisplayName','t=0')
% ylim([-0.5 1.5])
xlim([-1 5])
legend('t=0', 'Location','northeast')
title('MacCormack')
hold on
plot_times = [1, 2, 3, 4, 5, 6];
curr_plot=1;

dt_mac = 0.4*dx;
u_new = u;
u_star = zeros(1,N-1);
t=0;

while t<=6
    for j=1:N-1
        u_star(j) = u(j) - (dt_mac/dx)*0.5*(u(j+1)^2 - u(j)^2);
    end

    for j=2:N-1
        u_new(j) = 0.5*(u(j) + u_star(j)) - (dt_mac/(2*dx))*0.5*(u_star(j)^2 - u_star(j-1)^2);
    end

    u=u_new;
    t=t+dt_mac;

    if curr_plot <= length(plot_times) && t >= plot_times(curr_plot)
        plot(x_vals,u,'LineWidth',1.5,'Color',colors(curr_plot,:))

        u_ex = exact_solution(x_vals, plot_times_exact(curr_plot));
        plot(x_vals, u_ex, 'LineWidth', 1.5,'Color',colors_trans(curr_plot,:))
        hold on

        curr_plot = curr_plot+1;
    end
end
legend('t=0','t=1','','t=2','','t=3','','t=4','','t=5','','t=6','')

%% Godunov

dt_god = 0.4*dx;
u=u0;

figure;
plot(x_vals, u, 'LineWidth', 1.5, 'DisplayName', 't = 0')
ylim([-0.5 1.5])
legend('Location', 'northeast')
title('Godunov')
hold on
plot_times = [1, 2, 3, 4, 5, 6];
curr_plot = 1;

u_new = u;
t = 0;

while t <= 6
    F = zeros(1, N-1);
    for j = 1:N-1
        uL = u(j);
        uR = u(j+1);
        if uL <= uR
            if uL <= 0 && 0 <= uR
                F(j) = 0;
            elseif uR < 0
                F(j) = uR^2/2;
            else
                F(j) = uL^2/2;
            end
        else
            F(j) = max(uL^2, uR^2) / 2;
        end
    end

    for j = 2:N-1
        u_new(j) = u(j) - (dt_god/dx)*(F(j) - F(j-1));
    end

    u = u_new;
    t = t + dt_god;

    if curr_plot <= length(plot_times) && t >= plot_times(curr_plot)
        plot(x_vals,u,'LineWidth',1.5,'Color',colors(curr_plot,:))

        u_ex = exact_solution(x_vals, plot_times_exact(curr_plot));
        plot(x_vals, u_ex, 'LineWidth', 1.5,'Color',colors_trans(curr_plot,:))
        hold on

        curr_plot = curr_plot+1;
    end
end
legend('t=0','t=1','','t=2','','t=3','','t=4','','t=5','','t=6','')
