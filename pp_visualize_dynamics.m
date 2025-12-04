function pp_visualize_dynamics(folder, fileName)
    VideoOn = 1;
    load(fullfile(folder, fileName))
    if isfile(fullfile(folder, [fileName, '.mp4']))
        disp(['File exists. Skip ' fileName])
        return
    end
    fprintf(1, 'Now reading %s\n', fileName);
    numFrames = size(X, 1);
    
    % Create figure and scatter plot ONCE
    fig = figure('Color', 'w', 'Visible','off');
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
    
    if VideoOn
        videoFilename = fullfile(folder, [fileName, '.mp4']);
        v = VideoWriter(videoFilename, 'MPEG-4');
        v.FrameRate = 30;         % frames per second (adjust if needed)
        open(v);
    end
    
    % Animate by updating scatter data only
    for k = 1:numFrames
        particlePlot.XData = X(k,:);
        particlePlot.YData = Y(k,:);
        timeText.String = sprintf('t = %.2f', T(k)); % update time
        drawnow;   
        if VideoOn
            frame = getframe(fig);
            writeVideo(v, frame);
        end
    end
    
    if VideoOn
        close(v);
        disp(['✅ Video saved to: ', videoFilename]);
    end
    
    %%
    % figure
    % plot(T,U,T,K,T,E, 'LineWidth', 1.2)
    % legend('Potential Energy', 'Kinetic Energy', 'Total Energy')
end

