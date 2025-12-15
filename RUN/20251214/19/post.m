
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
    load([folder, fileName]);

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
savefig([folder exportfineName, '_std.fig'])

figure
scatter(Omega, Angular_Freq, 20, 'k', 'filled');
hold on
plot(Omega, Omega/2, 'r--');
holf off
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


