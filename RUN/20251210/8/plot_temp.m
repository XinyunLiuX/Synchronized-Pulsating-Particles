% clc; clear all

% addpath('RES')
% filename = 'res_N=20_rho=2.00_beta=1.4_epsilon=0.07_omega=4.60';
% load([filename, '.mat'])
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

%%
clc; clear all;
folder = 'RES/';
files = dir([folder, '*.mat']);
numFile = length(files);
Omega = zeros(numFile,1);
Epsilon = zeros(numFile,1);
Mean = zeros(numFile,1);
Std = zeros(numFile,1);
Angular_Freq = zeros(numFile,1);

for i = 1:numFile

    fileName = files(i).name;
    load([folder, fileName])

    Omega(i) = p.omega;
    Epsilon(i) = p.epsilon;

    xpost = x(floor(p.M/2):end,1);

    [up, lo] = envelope(xpost, 10, 'peak');
    env = up - lo;
    Mean(i) = mean(env);
    Std(i) = std(env);
    
    Angular_Freq(i) = 2*pi*dominantFrequence(xpost - mean(xpost), p.dt);

end

exportfineName = sprintf('res_N=%d_rho=%.2f_beta=%.2f_epsilon=%.2f', p.N, p.rho, p.beta, p.epsilon);

figure
scatter(Omega, Mean, 20, 'k', 'filled');
savefig([folder exportfineName, '_mean.fig'])

figure
scatter(Omega, Std, 20, 'k', 'filled');
hold on
plot(Omega, Omega/2, 'r--');
hold off
savefig([folder exportfineName, '_std.fig'])

figure
scatter(Omega, Angular_Freq, 20, 'k', 'filled');
savefig([folder exportfineName, '_freq.fig'])


function f_dominant = dominantFrequence(xpost, dt)

Fs = 1/dt;               % Sampling Frequency (Hz)
L = length(xpost);       % Length of the signal

% --- 3. Compute the FFT ---
% Use a power of 2 for NFFT for faster computation, though optional
NFFT = 2^nextpow2(L);
X = fft(xpost, NFFT);

% --- 4. Generate the Single-Sided Spectrum P1 ---
% Compute the two-sided spectrum P2
P2 = abs(X/L);

% Take the single-sided spectrum P1
% Only need the first half (up to the Nyquist frequency)
P1 = P2(1:NFFT/2+1);
P1(2:end-1) = 2*P1(2:end-1); % Multiply by 2 (except DC and Nyquist components)

% --- 5. Create the Frequency Vector f ---
f = Fs*(0:(NFFT/2))/NFFT;

% --- 6. Find the Dominant Frequency ---
% Find the index 'I' corresponding to the maximum magnitude in P1
[~, I] = max(P1);

% The frequency at that index is the dominant frequency
f_dominant = f(I);

end


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