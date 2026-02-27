clear; close all

%% problem 1 RAS plot

gamma = 1-sqrt(2)/2;
colors = 1/255*[0 0 255; 160 32 240; 255 0 0; 255 165 0; 0 255 255; 0 0 0];
real = linspace(-5,5,1001);
im = linspace(-5,5,1001);
[x,y] = meshgrid(real, im);
z = x + 1i*y;

num = 2 + (2-4*gamma)*z + (2*gamma^2 - 4*gamma + 1)*z.^2; % numerator
denom = 2*(1 - 2*gamma*z + gamma^2*z.^2); % denominator

R_z = num./denom; % R(z)
absR = abs(R_z); % |R(z)|

mask = double(absR <= 1);
mask(mask == 0) = NaN;

% Plot
figure(1)
pcolor(x, y, mask);
shading flat;
colormap([1 0.7 0.7]);
hold on

%boundary and axes
contour(x, y, absR, [1 1], 'k', 'LineWidth', 1.5);
xline(0, 'k--', 'LineWidth', 1);
yline(0, 'k--', 'LineWidth', 1);

xlabel('$\Re(z)$','Interpreter','latex'); ylabel('$\Im(z)$','Interpreter','latex');
title(sprintf('RAS of DIRK with \\gamma = 1 - \\surd2/2'));
legend('Stability region','Location', 'northwest');
grid on;
axis equal;
xlim([-4.5 4.5]);
ylim([-4.5 4.5]);

%% problem 4 - mu=100, Tmax=200, DIRK2 fixed

mu = 1e2;
Tmax = 200;
h = 1e-3;
y0 = [2; 0];

tic
[t1, sol1] = DIRK2_fixed(y0, h, Tmax, mu);
cpu1 = toc;

figure;
subplot(1,2,1);
plot(t1, sol1(:,1), 'b', 'LineWidth', 1.5);
hold on
plot(t1, sol1(:,2), 'r', 'LineWidth', 1.5);
xlabel('t','FontSize',14);
ylabel('x, y','FontSize',14);
title('Numerical Solution','FontSize',14);
legend('x(t)','y(t)','Location','northwest');
set(gca,'FontSize',13);
grid on; box on;

subplot(1,2,2);
plot(sol1(:,1), sol1(:,2), 'b', 'LineWidth', 1.5);
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
title('Phase Space','FontSize',14);
set(gca,'FontSize',13);
grid on; box on;

%% problem 4 - mu=1000, Tmax = 2000, DIRK2 fixed

mu2 = 1e3;
Tmax2 = 2e3;
h2 = 1e-3;

tic
[t2, sol2] = DIRK2_fixed(y0, h2, Tmax2, mu2);
cpu2 = toc;

figure;
subplot(1,2,1);
plot(t2, sol2(:,1), 'b', 'LineWidth', 1.5);
hold on
plot(t2, sol2(:,2), 'r', 'LineWidth', 1.5);
xlabel('t','FontSize',14);
ylabel('x, y','FontSize',14);
title('Numerical Solution','FontSize',14);
legend('x(t)','y(t)','Location','northwest');
set(gca,'FontSize',13);
grid on; box on;

subplot(1,2,2);
plot(sol2(:,1), sol2(:,2), 'b', 'LineWidth', 1.5);
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
title('Phase Space','FontSize',14);
set(gca,'FontSize',13);
grid on; box on;

%% problem 4 - mu=10^6, Tmax = 2*10^6, Adaptive DIRK2

mu3 = 1e6;
Tmax3 = 2e6;
atol = 1e-5;
rtol = 1e-5;
h_init = 1.0;

tic
[t3, sol3, nsteps3, nrej3] = DIRK2_adaptive(y0, h_init, Tmax3, mu3, atol, rtol);
cpu3 = toc;

figure;
subplot(1,2,1);
plot(t3, sol3(:,1), 'b', 'LineWidth', 1.5);
hold on
plot(t3, sol3(:,2), 'r', 'LineWidth', 1.5);
xlabel('t','FontSize',14);
ylabel('x, y','FontSize',14);
title('Numerical Solution','FontSize',14);
legend('x(t)','y(t)','Location','northwest');
set(gca,'FontSize',13);
grid on; box on;

subplot(1,2,2);
plot(sol3(:,1), sol3(:,2), 'b', 'LineWidth', 1.5);
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
title('Phase Space','FontSize',14);
set(gca,'FontSize',13);
grid on; box on;