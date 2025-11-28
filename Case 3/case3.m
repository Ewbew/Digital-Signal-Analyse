
clear; clc; close all;

[x, fs] = audioread("tale_tone_48000.wav");

% Plot the time domain signal
t = (0:length(x)-1)/fs; % Time vector
figure;
plot(t, x);
title('Time Domain Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Compute and plot the amplitude spectrum
n = length(x);
f = (0:n-1)*(fs/n); % Frequency vector
xfft = abs(fft(x));
figure;
subplot(2,1,1); % full scale amplitude spectrum
plot(f, xfft);
xlim([0 fs/2]); % Limit x-axis to Nyquist frequency
title('Amplitude Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(2,1,2); % limited amplitude spectrum
plot(f, xfft);
xlim([0 3000]); % Limit x-axis to 3000Hz
title('Amplitude Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Plot the spectrogram
figure;
subplot(2,1,1);
spectrogram(x, 256, 250, 256, fs, 'yaxis');
title('Spectrogram');
xlabel('Time (s)');
ylabel('Frequency (kHz)');
colorbar;


% Compute and plot the second spectrogram with the same window size but
% limited y-axis
subplot(2,1,2);
spectrogram(x, 256, 250, 256, fs, 'yaxis');
title('Spectrogram with zoom on y-axis');
xlabel('Time (s)');
ylabel('Frequency (kHz)');
ylim([0 5]); % Limit y-axis to 5 kHz
colorbar;

%%
soundsc(x, fs);

%%

fp = 785;  % Pole frequency - aflæst i amplitude spektrum
rp = 0.9; % Pole radius

fn = 785; % Zero frequency 
rn = 1.0; % Zero radius

wp = 2*pi*fp/fs;
wn = 2*pi*fn/fs;

% Coefficients numerator 
b0 = 1;
b1 = -2*rn*cos(wn);
b2 = rn^2;
b = [b0 b1 b2];

% Coefficients denominator
a0 = 1;
a1 = -2*rp*cos(wp);
a2 = rp^2;
a = [a0 a1 a2];

% Pole-zero plane
figure;
zplane(b, a);
title('Pole-Zero Diagram');

% Filter characteristics
figure;
[H, f] = freqz(b, a, fs); % Get frequency response
f = f * fs / (2 * pi); % Convert normalized frequency to Hz
plot(f, 20*log10(abs(H))); % Plot magnitude response in dB
title('Filter Frequency Response');
xlabel('Frequency (Hz)'); % Set x-axis label to Frequency (Hz)
ylabel('Magnitude (dB)'); % Set y-axis label to Magnitude (dB)

%%

% Apply the filter to the signal x
filtered_x = filter(b, a, x);

% Plot the original and filtered signals for comparison
figure;
subplot(2,1,1);
plot(t, x);
title('Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t, filtered_x);
title('Filtered Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Compute and plot the frequency response of the original signal
[H_orig, f_orig] = freqz(x, 1, fs); % Get frequency response of original signal
f_orig = f_orig * fs / (2 * pi); % Convert normalized frequency to Hz

% Plot the frequency response of original and filtered signals
figure;
subplot(2,1,1);
plot(f_orig, abs(H_orig)); % Plot magnitude response
title('Original Signal Frequency Response');
xlabel('Frequency (Hz)');
xlim([0 4000]);
ylabel('Magnitude');

subplot(2,1,2);
[H_filt, f_filt] = freqz(filtered_x, 1, fs); % Get frequency response of filtered signal
f_filt = f_filt * fs / (2 * pi); % Convert normalized frequency to Hz
plot(f_filt, abs(H_filt)); % Plot magnitude response
title('Filtered Signal Frequency Response');
xlabel('Frequency (Hz)');
xlim([0 5000]);
ylabel('Magnitude');

% Play the filtered signal
soundsc(filtered_x, fs);
% Wait until the sound is done
pause(length(filtered_x)/fs);
% Play the original signal

soundsc(x, fs);