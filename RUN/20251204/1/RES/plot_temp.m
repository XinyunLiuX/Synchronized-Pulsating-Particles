% clc; clear all
r0 = 2^(1/6);

plot(T, X(:,2)-X(:,1)); hold on
yline(r0)
legend(['\beta = ' num2str(p.beta)])
set(gca, 'FontSize', 15)