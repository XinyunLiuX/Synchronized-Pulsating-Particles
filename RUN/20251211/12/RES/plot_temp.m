clc; clear all
load('res_N=20_rho=2.00_beta=1.4_epsilon=0.04_omega=4.60.mat')

figure
subplot 511
plot(T, x(:,1))
subplot 512
plot(T, x(:,2))
subplot 513
plot(T, x(:,3))
subplot 514
plot(T, x(:,4))
subplot 515
plot(T, x(:,5))


%%

VideoOn = 1;
if VideoOn
    exportfineName = sprintf('res_N=%d_rho=%.2f_beta=%.2f_epsilon=%.2f_omega=%.2f', p.N, p.rho, p.beta, p.epsilon, p.omega);
    vid = VideoWriter([exportfineName '.mp4'],'MPEG-4');
    vid.FrameRate = 40;
    vid.Quality = 95;
    open(vid);
end


T_start = floor(19*p.M/20);
% Create figure and scatter plot ONCE
fig = figure('Color', 'w', 'Visible','on');
set(fig, 'Units', 'normalized', 'OuterPosition', [0 0.4 0.5 0.25]);
particlePlot = scatter(x(T_start,:), y(T_start,:), 30, 'filled', 'k');
axis equal;
axis([0 p.L p.L/2 - 1 p.L/2+1]);
box on;
set(gca, 'XTick', [], 'YTick', []);
title('2D Lennard-Jones Simulation','FontSize', 20);

% Add time stamp text
timeText = text(0.02*p.L, 0.95*p.L, sprintf('t = %.2f', T(1)), ...
    'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

drawnow;



% Animate by updating scatter data onlyclose(vid)
for k = T_start:1:length(T)
    k
    particlePlot.XData = x(k,:);
    particlePlot.YData = y(k,:);
    timeText.String = sprintf('t = %.2f', T(k)); % update time
    drawnow;   
    if VideoOn
        % ---- Save high-resolution PNG ----
        pngName = 'videomaker.jpg';
        exportgraphics(fig, pngName, 'Resolution', 150);  % TRUE 300 dpi
        img = imread(pngName);
        writeVideo(vid, img);
    end
end


close(vid)
