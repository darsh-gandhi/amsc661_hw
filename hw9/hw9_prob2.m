clear; close all;

%% base params

x0 = -5:0.5:5; x0((end+1)/2) = [];
rho_0 = [log(0.1); log(0.9)];

t=linspace(0,10,101);

%% plotting

figure
plot(zeros(1,101),t,'k','LineWidth',2.0)
hold on
for i=1:length(x0)
    if x0(i) < 0
        characteristic = x0(i)*ones(1,length(t)) - (rho_0(1) + 1)*t;
        plot(characteristic,t, 'LineWidth',1.5)
    else
        characteristic = x0(i)*ones(1,length(t)) - (rho_0(2) + 1)*t;
        plot(characteristic,t, 'LineWidth',1.5)
    end
end
grid on
xlim([-5.1 5.1])
xlabel('x')
ylabel('t')


%% prob 2c

clear;

t_prime_eq = @(x) x*(1 + 1.8/pi * atan(x)) + 0.9/pi;
root = fzero(t_prime_eq,0)

t_shock = @(x) (0.5 + (1.8/pi)*atan(x))/((0.9)/(pi*(1+x^2))); t_shock(root)