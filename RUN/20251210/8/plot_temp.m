% clc; clear all
% synchronized_pulsating_particles_local

% t = 1;
% scatter(X(t,:), Y(t,:), '.')
% axis([0,L,L/2-0.1,L/2+0.1])

% addpath('RES')
% filename = 'res_N=20_rho=2.00_beta=1.4_epsilon=0.07_omega=4.80';
% load([filename, '.mat'])
% %%
% subplot 511
% plot(T, x(:,1))
% subplot 512
% plot(T, x(:,2))
% subplot 513
% plot(T, x(:,3))
% subplot 514
% plot(T, x(:,4))
% subplot 515
% plot(T, x(:,5))
% 
% figure

clc; clear all;
folder = 'RES/';
files = dir([folder, '*.mat']);
numFile = length(files);
Omega = zeros(numFile,1);
Epsilon = zeros(numFile,1);
Mean = zeros(numFile,1);
Std = zeros(numFile,1);

for i = 1:numFile

    fileName = files(i).name;
    load([folder, fileName])

    Omega(i) = p.omega;
    Epsilon(i) = p.epsilon;

    [up, lo] = envelope(x(:,1), 10, 'peak');
    secondPart = up-lo;
    secondPart = secondPart(floor(length(secondPart)/2):end);
    Mean(i) = mean(secondPart);
    Std(i) = std(secondPart);

end
figure
scatter(Omega, Mean, 20, 'k', 'filled');
savefig([folder fileName(1:end-15), '_mean.fig'])

figure
scatter(Omega, Std, 20, 'k', 'filled');
savefig([folder fileName(1:end-15), '_std.fig'])


% %%
% VideoOn = 1;
% if VideoOn
%     Vid = VideoWriter([filename, '.mp4'], 'MPEG-4');
%     Vid.FrameRate = 20;
%     Vid.Quality = 95;
%     open(Vid)
% end
% 
% % Create figure and scatter plot ONCE
% fig = figure('Color', 'w', 'Visible','on');
% set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);  % full screen
% particlePlot = scatter(x(1,:), y(1,:), 20, 'filled');
% axis equal;
% axis([0 p.L 0 p.L]);
% box on;
% set(gca, 'XTick', [], 'YTick', []);
% title('2D Lennard-Jones Simulation');
% 
% % Add time stamp text
% timeText = text(0.02*p.L, 0.95*p.L, sprintf('t = %.2f', T(1)), ...
%     'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k', ...
%     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
% 
% drawnow;
% 
% 
% 
% % Animate by updating scatter data only
% for k = 1:30:length(T)
%     particlePlot.XData = x(k,:);
%     particlePlot.YData = y(k,:);
%     timeText.String = sprintf('t = %.2f', T(k)); % update time
%     drawnow;   
%     if VideoOn
%         frame = getframe(fig);
%         writeVideo(Vid, frame);
%     end
% end
% 
% if VideoOn
%     close(Vid) 
%     fprintf('Video "%s.mp4" successfully created.\n', filename);
% end