% clc; clear all
% synchronized_pulsating_particles_local
r0 = 2^(1/6);
yyaxis left
h1 = plot(T, X(:,2)-X(:,1)); 
yline(r0)
legend(['\beta = ' num2str(p.beta)],'')
set(gca, 'FontSize', 15)
yyaxis right
h2 = plot(T, 1 + p.epsilon*sin(p.omega*T));

legend([h1 h2], ...
       {[['\beta = ' num2str(p.beta)],''], ...
        'radius'}, ...
       'Location','best');


title(['\omega='  num2str(p.omega), '\epsilon =' num2str(p.epsilon)])