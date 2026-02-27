clear; close;

%% hw3 prob 2

L    = 1e4;
phi  = @(t) sin(t + pi/4);
dphi = @(t) cos(t + pi/4);

% Exact solution
exact = @(t, y0) exp(-L*t)*(y0 - phi(0)) + phi(t);

% h values: h = 10^(-p), p in {1, 1+d, 1+2d, ..., 6}, d = 5/24
d = 5/24;
p_vals = 1 : d : 6;
h_vals = 10.^(-p_vals);
Nh = length(h_vals);

Tmax = 10;

y0_a = sin(pi/4);

err_DIRK2_a  = zeros(1, Nh);
err_DIRKo3_a = zeros(1, Nh);

for ih = 1:Nh
    h = h_vals(ih);
    Nsteps = ceil(Tmax / h);
    t_arr  = (0:Nsteps) * h;

    % DIRK2
    y_DIRK2  = run_DIRK2(y0_a, h, Nsteps, L, phi, dphi);
    y_exact  = exact(t_arr, y0_a);
    err_DIRK2_a(ih) = max(abs(y_DIRK2 - y_exact));

    % DIRKo3
    y_DIRKo3 = run_DIRKo3(y0_a, h, Nsteps, L, phi, dphi);
    err_DIRKo3_a(ih) = max(abs(y_DIRKo3 - y_exact));
end

%%
% slope 1
C1_DIRK2 = err_DIRK2_a(1)/h_vals(1);
C1_DIRKo3 = err_DIRKo3_a(1)/h_vals(1);
C1_avg = (C1_DIRK2 + C1_DIRKo3)/2;

%slope 2
C2 =err_DIRK2_a(end)/h_vals(end)^2;

figure;
loglog(h_vals, err_DIRK2_a,'b-o', 'LineWidth', 2,'MarkerSize', 7, 'DisplayName','DIRK2');
hold on;
loglog(h_vals, err_DIRKo3_a,'r-s', 'LineWidth', 2, 'MarkerSize', 7, 'DisplayName','DIRKo3');
loglog(h_vals, C1_avg*h_vals, 'k--', 'LineWidth', 1.5,'DisplayName','Slope 1');
loglog(h_vals, C2*h_vals.^2, 'k:','LineWidth', 1.5,'DisplayName','Slope 2');
xlabel('h', 'FontSize', 14);
ylabel('e(h)', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 12);
set(gca, 'FontSize', 13);
grid on;

%%
y0_b = sin(pi/4) + 10;

err_DIRK2_b  = zeros(1, Nh);
err_DIRKo3_b = zeros(1, Nh);

for ih = 1:Nh
    h = h_vals(ih);
    Nsteps = ceil(Tmax / h);
    t_arr  = (0:Nsteps) * h;

    y_DIRK2  = run_DIRK2(y0_b, h, Nsteps, L, phi, dphi);
    y_exact  = exact(t_arr, y0_b);
    err_DIRK2_b(ih) = max(abs(y_DIRK2 - y_exact));

    y_DIRKo3 = run_DIRKo3(y0_b, h, Nsteps, L, phi, dphi);
    err_DIRKo3_b(ih) = max(abs(y_DIRKo3 - y_exact));
end

%%

C1b_DIRK2 = err_DIRK2_b(1)/h_vals(1);
C1b_DIRKo3 = err_DIRKo3_b(1)/h_vals(1);

figure;
loglog(h_vals, err_DIRK2_b,'b-o','LineWidth', 2, 'MarkerSize', 7, 'DisplayName', 'DIRK2');
hold on;
loglog(h_vals, err_DIRKo3_b, 'r-s','LineWidth', 2, 'MarkerSize', 7,'DisplayName', 'DIRKo3');
loglog(h_vals, C1b_DIRK2*h_vals,'b--', 'LineWidth', 1.5, 'DisplayName','Slope 1');
loglog(h_vals, C1b_DIRKo3*h_vals,'r--', 'LineWidth', 1.5, 'DisplayName','Slope 2');
xlabel('h', 'FontSize', 14);
ylabel('e(h)','FontSize', 14);
legend('Location', 'northwest', 'FontSize', 12);
set(gca, 'FontSize', 13);
grid on

%%
h_vec = [1e-1, 1e-2, 1e-3];
Tmax_plot = [10, 1];
% line_colors = {'b', 'r', 'g'};
colors = 1/255*[0 0 255; 160 32 240; 255 0 0; 255 165 0; 0 255 255; 0 0 0];

for iTmax = 1:2
    Tv = Tmax_plot(iTmax);

    figure;
    for im = 1:2
        subplot(1,2,im);
        hold on;
        for ih = 1:3
            h = h_vec(ih);
            Ns = ceil(Tv / h);
            t_arr = (0:Ns) * h;
            y_ex= exact(t_arr, y0_b);

            if im == 1
                y_num = run_DIRK2(y0_b, h, Ns, L, phi, dphi);
                method_str = 'DIRK2';
            else
                y_num = run_DIRKo3(y0_b, h, Ns, L, phi, dphi);
                method_str = 'DIRKo3';
            end

            e_t = abs(y_num - y_ex);

            semilogy(t_arr, e_t, 'Color',colors(ih,:), 'LineWidth', 1.5,'DisplayName', sprintf('h = 10^{%d}', round(log10(h))));
        end
        xlabel('t', 'FontSize', 14);
        ylabel('|e(t)|', 'FontSize', 14);
        title(sprintf('%s,  T_{max} = %d', method_str, Tv), 'FontSize', 14);
        legend('Location', 'best', 'FontSize', 12);
        set(gca, 'FontSize', 13);
        grid on;
    end
end

function rhs = func(t, y, L, phi, dphi)
    rhs = -L*(y - phi(t)) + dphi(t);
end

function y_all = run_DIRK2(y0, h, Nsteps, L, phi, dphi)
    gamma = 1 - 1/sqrt(2);
    tol = 1e-14;
    itermax = 50;

    y_all = zeros(1, Nsteps+1);
    y_all(1) = y0;
    y = y0;

    for j = 1:Nsteps
        t_n = (j-1)*h;
        y = DIRK2step_PR(y, h, t_n, gamma, L, phi, dphi, tol, itermax);
        y_all(j+1) = y;
    end
end

function ynew = DIRK2step_PR(y, h, tn, gamma, L, phi, dphi, tol, itermax)
    % stage1
    t1 = tn + gamma*h;
    k1 = func(t1, y, L, phi, dphi);
    for j = 1:itermax
        aux = y + h*gamma*k1;
        F = k1 - func(t1, aux, L, phi, dphi);

        DF   = 1 + h*gamma*L;
        k1   = k1 - F/DF;
        if abs(k1 - func(t1, y + h*gamma*k1, L, phi, dphi)) < tol
            break;
        end
    end

    % stage 2
    t2 = tn + h;
    y2 = y + h*(1-gamma)*k1;
    k2 = k1;
    for j = 1:itermax
        aux  = y2 + h*gamma*k2;
        F    = k2 - func(t2, aux, L, phi, dphi);
        DF   = 1 + h*gamma*L;
        k2   = k2 - F/DF;
        if abs(k2 - func(t2, y2 + h*gamma*k2, L, phi, dphi)) < tol
            break;
        end
    end

    ynew = y2 + h*gamma*k2;
end

function y_all = run_DIRKo3(y0, h, Nsteps, L, phi, dphi)
    gamma = 0.5 + sqrt(3)/6;
    tol = 1e-14;
    itermax = 50;

    y_all = zeros(1, Nsteps+1);
    y_all(1) = y0;
    y = y0;

    for j = 1:Nsteps
        t_n = (j-1)*h;
        y = DIRKo3step_PR(y, h, t_n, gamma, L, phi, dphi, tol, itermax);
        y_all(j+1) = y;
    end
end

function ynew = DIRKo3step_PR(y, h, tn, gamma, L, phi, dphi, tol, itermax)
    % stage 1
    t1 = tn + gamma*h;
    k1 = func(t1, y, L, phi, dphi);
    for j = 1:itermax
        aux = y + h*gamma*k1;
        F = k1 - func(t1, aux, L, phi, dphi);
        DF = 1 + h*gamma*L;
        k1 = k1 - F/DF;
        if abs(k1 - func(t1, y + h*gamma*k1, L, phi, dphi)) < tol
            break;
        end
    end

    % stage 2
    t2 = tn + (1-gamma)*h;
    y2 = y + h*(1-2*gamma)*k1;   % explicit part for stage 2
    k2 = k1;
    for j = 1:itermax
        aux = y2 + h*gamma*k2;
        F = k2 - func(t2, aux, L, phi, dphi);
        DF = 1 + h*gamma*L;
        k2 = k2 - F/DF;
        if abs(k2 - func(t2, y2 + h*gamma*k2, L, phi, dphi)) < tol
            break;
        end
    end

    % update
    ynew = y + h*(0.5*k1 + 0.5*k2);
end