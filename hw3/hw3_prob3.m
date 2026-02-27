%% problem 3

Tmax = 100;
tol  = 1e-14;
itermax = 50;
y0 = [1.0; 0.0; 0.0];

hvals = [1e-3, 1e-2, 1e-1];
method_names = ["DIRK2", "DIRKo3", "BDF2"];
colors = 1/255*[0 0 255; 160 32 240; 255 0 0; 255 165 0; 0 255 255; 0 0 0];
cpu_times = zeros(3, 3);

figure;
for i = 1:length(hvals)
    h = hvals(i);
    Nsteps = ceil(Tmax/h);
    t = (0:Nsteps)'*h;

    % dirk2
    tic
    dirk2_sol = zeros(Nsteps+1,3);
    dirk2_sol(1,:) = y0';
    for j = 1:Nsteps
        dirk2_sol(j+1,:) = DIRK2step(dirk2_sol(j,:)', h, tol, itermax)';
    end
    cpu_times(1,i) = toc;

    % dirk3
    tic
    dirk3_sol = zeros(Nsteps+1,3);
    dirk3_sol(1,:) = y0';
    for j = 1:Nsteps
        dirk3_sol(j+1,:) = DIRKo3step(dirk3_sol(j,:)', h, tol, itermax)';
    end
    cpu_times(2,i) = toc;

    % bdf2
    tic
    bdf2_sol = zeros(Nsteps+1,3);
    bdf2_sol(1,:) = y0';
    bdf2_sol(2,:) = DIRK2step(bdf2_sol(1,:)', h, tol, itermax)';
    for j = 2:Nsteps
        bdf2_sol(j+1,:) = BDF2step(bdf2_sol(j,:)', bdf2_sol(j-1,:)', h, tol, itermax)';
    end
    cpu_times(3,i) = toc;

    fsz = 14;
    var_names = {'x','y','z'};
    sols = {dirk2_sol, dirk3_sol, bdf2_sol};

    for k = 1:3
        subplot(3,3,(k-1)*3 + i);
        hold on;
        for l = 1:3
            plot(t, sols{l}(:,k),'Color', colors(l,:),'LineWidth', 1.5,'DisplayName', method_names(l));
        end
        xlabel('t','FontSize',fsz);
        ylabel(var_names{k},'FontSize',fsz);
        title(sprintf('h = 10^{%d}', round(log10(h))),'FontSize',fsz);
        if k==1 && i==1
            legend('Location','best','FontSize',10);
        end
        set(gca,'FontSize',fsz);
        box on;
    end
end

figure;
fsz = 16;
for l = 1:3
    loglog(hvals, cpu_times(l,:), '-o', 'Color', colors(l,:), 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', method_names(l));
    hold on;
end
xlabel('h','FontSize',fsz);
ylabel('CPU time (s)','FontSize',fsz);
title("CPU Time",'FontSize',fsz);
legend('Location','best','FontSize',fsz);
set(gca,'FontSize',fsz,'XScale','log','YScale','log');
grid on;
box on;

function dy = func(y)
    a = 0.04; b = 1.0e4; c = 3.0e7;
    byz = b*y(2)*y(3);
    cy2 = c*y(2)^2;
    ax  = a*y(1);
    dy = [-ax + byz; ax - byz - cy2; cy2];
end

function J = Jac(y)
    a = 0.04; b = 1.0e4; c = 3.0e7;
    J = zeros(3);
    J(1,1) = -a;
    J(1,2) =  b*y(3);
    J(1,3) =  b*y(2);
    J(2,1) =  a;
    J(2,2) = -b*y(3) - 2*c*y(2);
    J(2,3) = -b*y(2);
    J(3,2) =  2*c*y(2);
end

function ynew = DIRK2step(y, h, tol, itermax)
    gamma = 1.0 - 1.0/sqrt(2);

    % stage 1
    k1 = func(y);
    for j = 1:itermax
        k1 = NewtonDIRK2(y, h, k1, gamma);
        if norm(k1 - func(y + h*gamma*k1)) < tol
            break;
        end
    end

    % stage 2
    y2 = y + h*(1-gamma)*k1;
    k2 = k1;
    for j = 1:itermax
        k2 = NewtonDIRK2(y2, h, k2, gamma);
        if norm(k2 - func(y2 + h*gamma*k2)) < tol
            break;
        end
    end

    ynew = y2 + h*gamma*k2;
end

function knew = NewtonDIRK2(y, h, k, gamma)
    aux = y + h*gamma*k;
    F   = k - func(aux);
    DF  = eye(3) - h*gamma*Jac(aux);
    knew = k - DF\F;
end

function ynew = DIRKo3step(y, h, tol, itermax)
    gamma = 0.5 + sqrt(3)/6;

    % stage 1
    k1 = func(y);
    for j = 1:itermax
        k1 = NewtonDIRKo3(y, h, k1, gamma);
        if norm(k1 - func(y + h*gamma*k1)) < tol
            break;
        end
    end

    % stage 2
    y2 = y + h*(1 - 2*gamma)*k1;
    k2 = func(y2);
    for j = 1:itermax
        k2 = NewtonDIRKo3(y2, h, k2, gamma);
        if norm(k2 - func(y2 + h*gamma*k2)) < tol
            break;
        end
    end

    ynew = y + h*(0.5*k1 + 0.5*(y2/h + k2 - y/h));
    Y2 = y2 + h*gamma*k2;
    k2val = func(Y2);
    ynew = y + h*(0.5*k1 + 0.5*k2val);
end

function knew = NewtonDIRKo3(y, h, k, gamma)
    aux = y + h*gamma*k;
    F   = k - func(aux);
    DF  = eye(3) - h*gamma*Jac(aux);
    knew = k - DF\F;
end

function ynew = BDF2step(yn, ynm1, h, tol, itermax)
    rhs = 2*yn - 0.5*ynm1;
    ynew = yn;
    for j = 1:itermax
        F  = 1.5*ynew - rhs - h*func(ynew);
        DF = 1.5*eye(3) - h*Jac(ynew);
        delta = DF\F;
        ynew = ynew - delta;
        if norm(delta) < tol
            break;
        end
    end
end