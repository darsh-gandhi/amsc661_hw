clear; close all

%% prob 5 part c

function dydt = grav_ode(t,y)
%GRAV_ODE := defines the ode we care about
% Inputs:
%   t <- time
%   y <- states
%   z <- params
% Outputs:
%   y <- states

dydt = zeros(4,1);

denom = y(1)^2 + y(2)^2;

dydt(1) = y(3);
dydt(2) = y(4);
dydt(3) = -y(1)/denom;
dydt(4) = -y(2)/denom;

end

function u = num_scheme(f,h,tk,uk,uk1)
%NUM_SCHEME := numerical integrator developed in earlier part of problem
% Inputs:
%   h <- stepsize
%   tk <- time
%   uk <- u_k (approx sol at u(t_k))
%   uk1 <- u_{k-1} (approx sol at u(t_{k-1}))
% Outputs:
%   y <- states


%   fn <- function eval at t_n
%   fn1 <- function eval at t_{n-1}

u = -4*uk - 5*uk1 + h*(4*grav_ode(tk,uk) + 2*grav_ode(tk-h,uk1));

end

% initializations
N=[20 40 80]; % subintervals
tf_norm1 = -1.0*ones(1,3);
figure(1)
hold on

for j = 1:length(N)
    N_j = N(j);
    h=2*pi/N_j; % step size
    t0 = 0; tf = 4*pi; tspan = (t0:h:tf)';
    init = [1,0,0,1]'; % init cond
    
    u = -1.0*ones(4,tf/h);
    u(:,1)=init; %u(0)
    u(:,2)=init;
    uk1=init; %u(h)
    uk = init;
    for i=3:tspan(end)/h
        tk = i*h;
        uout = num_scheme(@(t,y)grav_ode(t,y),h,tk,uk,uk1);
        u(:,i) = uout;
        uk1=uk;
        uk = uout;
    end
    
    tf_norm1(j) = norm(u(:,end));
    
    plot(u(1,:),u(2,:),'LineWidth',2.0)
end
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
legend('N=20','N=40','N=80')
hold off

%%

function u = midpoint_fe_predict(f,h,tk,uk)
%MIDPOINT_FE_PREDICT := midpoint rule with forward euler predictor
% Inputs:
%   h <- stepsize
%   tk <- time
%   uk <- u_k
% Outputs:
%   u <- u_{k+1}

p=grav_ode(tk,uk);
u = uk + h*grav_ode(tk + h/2, uk + h/2*p);

end


N=[20 40 80]; % subintervals
tf_norm2 = -1.0*ones(1,3);
figure(2)
hold on

for j = 1:length(N)
    N_j = N(j);
    h=2*pi/N_j; % step size
    t0 = 0; tf = 4*pi; tspan = (t0:h:tf)';
    init = [1,0,0,1]'; % init cond
    
    u = -1.0*ones(4,tf/h);
    u(:,1)=init; %u(0)
    u(:,2)=init;
    uk1=init; %u(h)
    uk = init;
    for i=3:tspan(end)/h
        tk = i*h;
        uout = midpoint_fe_predict(@(t,y)grav_ode(t,y),h,tk,uk);
        u(:,i) = uout;
        uk1=uk;
        uk = uout;
    end
    
    tf_norm2(j) = norm(u(:,end));

    plot(u(1,:),u(2,:),'LineWidth',2.0)
end
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
legend('N=20','N=40','N=80')
hold off