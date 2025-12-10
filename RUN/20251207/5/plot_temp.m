% clc; clear all
% synchronized_pulsating_particles_local

t = 1;
% scatter(X(t,:), Y(t,:), '.')
% axis([0,L,L/2-0.1,L/2+0.1])
subplot 511
plot(T, X(:,1))
subplot 512
plot(T, X(:,2))
subplot 513
plot(T, X(:,3))
subplot 514
plot(T, X(:,4))
subplot 515
plot(T, X(:,5))

%%

    VideoOn = 0;

    % Create figure and scatter plot ONCE
    fig = figure('Color', 'w', 'Visible','on');
    set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);  % full screen
    particlePlot = scatter(X(1,:), Y(1,:), 20, 'filled');
    axis equal;
    axis([0 p.L 0 p.L]);
    box on;
    set(gca, 'XTick', [], 'YTick', []);
    title('2D Lennard-Jones Simulation');

    % Add time stamp text
    timeText = text(0.02*p.L, 0.95*p.L, sprintf('t = %.2f', T(1)), ...
        'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

    drawnow;



    % Animate by updating scatter data only
    for k = 1:10:length(T)
        particlePlot.XData = X(k,:);
        particlePlot.YData = Y(k,:);
        timeText.String = sprintf('t = %.2f', T(k)); % update time
        drawnow;   
        if VideoOn
            frame = getframe(fig);
            writeVideo(v, frame);
        end
    end

    

