clear; close all

%% funcs

function H = hamiltonian(u,v,x,y)
    H = 0.5*u^2 + 0.5*v^2 - 1/sqrt(x^2+y^2);
end

function [gx,gy] = gradU(x,y)
    denom = (x^2 + y^2)^1.5;
    gx = x/denom;
    gy = y/denom;
end

function z_n1 = stoermer_verlet(z_n,h)
%stoermer_verlet := one-step implementation of stoermer-verlet
%   input: z_n <- states of ode system at time n;
%       h   <- stepsize
%   output: z_n1 <- stats of ode system at time n+1

    % first half step
    [gx,gy] = gradU(z_n(3),z_n(4));
    u_n5 = z_n(1) - 0.5*h*gx;
    v_n5 = z_n(2) - 0.5*h*gy;
    x_n1 = z_n(3) + h*u_n5;
    y_n1 = z_n(4) + h*v_n5;
    
    % second half step
    [gx2,gy2] = gradU(x_n1,y_n1);
    u_n1 = u_n5 - 0.5*h*gx2;
    v_n1 = v_n5 - 0.5*h*gy2;
    
    z_n1 = [u_n1,v_n1,x_n1,y_n1];

end

%% ics and params
a = 4/3;
T = 2*pi*a^(3/2);   % period
h = 0.01*T;         %stepsize
t = 0:h:10*T;

%ICs
u = 0;
v = 0.5;
x = 2;
y = 0;

%% apply scheme

%initialize important vectors as -1s
z = -1.0*ones(length(t),4);
z(1,:) = [u,v,x,y];
hamil = -1.0*ones(length(t),1);
hamil(1) = hamiltonian(u,v,x,y);

for i=2:length(t)
    z(i,:) = stoermer_verlet(z(i-1,:),h);
    hamil(i) = hamiltonian(z(i,1),z(i,2),z(i,3),z(i,4));
end

%% plotting

colors = 1/255*[56 182 255; 255 222 89; 0 191 99]; % light blue, yellow, green

figure(1)
plot(z(:,3),z(:,4),'LineWidth',1.0,'Color',colors(1,:))
hold on
% plot(2,0,'r*','MarkerSize',18)
% hold on
% plot(-2/3,0,'r*','MarkerSize',18)
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
xlim([-1 2+1/3])
xline(2,'--','Color',colors(3,:),'LineWidth',2.0);
xline(-2/3,'--','Color',colors(2,:),'LineWidth',2.0);
legend('orbit','x_{max}','x_{min}')
grid on;
set(gca,'FontSize',14)

figure(2)
plot(t,hamil,'LineWidth',1.5, 'Color',colors(3,:))
xlabel('Time')
ylabel('$H(u,v,x,y)$','Interpreter','latex')
xlim([0,10*T+1])
set(gca,'FontSize',14)