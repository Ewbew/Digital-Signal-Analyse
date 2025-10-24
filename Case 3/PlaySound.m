%% Case 3 – Opgave 1: Analyse af inputsignal
% Viser tids-signal, amplitude-spektrum og spektrogram
% Estimerer også tonefrekvensen vha. peak i spektrummet.

clear; clc; close all;

%% 1) Indlæs lyd
fname = 'tale_tone_48000.wav';
[x, fs] = audioread(fname);

% Gør mono hvis stereo
if size(x,2) > 1
    x = mean(x,2);
end

% Fjern DC (hjælper peak-detektion)
x = x - mean(x);

N = numel(x);
t = (0:N-1)/fs;

%% 2) Plot tids-signal (viser første 0.2 s for overblik + hele signalet)
figure('Name','Tidsdomæne');
subplot(2,1,1);
plot(t, x, 'LineWidth', 0.8); grid on;
xlabel('Tid [s]'); ylabel('Amplitude');
title('Hele signalet');

subplot(2,1,2);
Tview = 0.20;                         % vis de første 0.2 s
Nv = min(N, round(Tview*fs));
plot(t(1:Nv), x(1:Nv), 'LineWidth', 0.8); grid on;
xlabel('Tid [s]'); ylabel('Amplitude');
title(sprintf('Udsnit (0–%.0f ms)', 1000*Tview));

%% 3) Amplitudespektrum (enkeltsidet)
% Brug et Hann-vindue og zero-padding for pænere spektrum og lidt højere frekvensopløsning
w = hann(N);
xw = x .* w;

% Zero-pad til mindst 2^nextpow2 for finere opløsning
Nfft = 2^nextpow2(N);
X = fft(xw, Nfft);
f = (0:Nfft-1)*(fs/Nfft);

% Enkeltsidet spektrum
half = 1:floor(Nfft/2);
f1 = f(half);
Mag = abs(X(half)) / (sum(w)/2);   % normaliser cirka for vindues-energi

figure('Name','Amplitudespektrum');
plot(f1, Mag, 'LineWidth', 0.8); grid on;
xlim([0 fs/2]);
xlabel('Frekvens [Hz]'); ylabel('|X(f)| (rel.)');
title('Amplitudespektrum (enkeltsidet)');

%% 4) Estimér tonefrekvens (ignorer DC og meget lave frekvenser)
f_lowcut = 20;                                   % undgå DC og rumlen
idx = f1 >= f_lowcut;
Mag_roi = Mag(idx);
f_roi = f1(idx);

% Find største peak
% (Hvis du har Signal Processing Toolbox: brug findpeaks; ellers max)
try
    [pks, locs] = findpeaks(Mag_roi, 'NPeaks', 3, 'SortStr', 'descend', ...
                            'MinPeakDistance', round(5/(fs/Nfft))); % ~5 Hz afstand
    f_tone = f_roi(locs(1));
catch
    [~, k] = max(Mag_roi);
    f_tone = f_roi(k);
end

fprintf('Estimeret tonefrekvens: %.2f Hz\n', f_tone);

% Markér peak i plottet
hold on;
yline = interp1(f1, Mag, f_tone, 'linear', 'extrap');
plot(f_tone, yline, 'o', 'MarkerSize', 6, 'LineWidth', 1.0);
legend('Spektrum', sprintf('Peak ~ %.2f Hz', f_tone), 'Location', 'northeast');

%% 5) Spektrogram (tids-frekvens oversigt)
% Lidt længere vindue giver bedre frekvensopløsning; overlap for glat billede
winLen = 2048;
ovlp   = round(0.5*winLen);
nfft_s = 4096;

figure('Name','Spektrogram');
spectrogram(x, hann(winLen), ovlp, nfft_s, fs, 'yaxis');
title('Spektrogram');
colorbar;

%% 6) Gem figurer (valgfrit)
% saveas(figure(1), 'case3_time.png');
% saveas(figure(2), 'case3_spectrum.png');
% saveas(figure(3), 'case3_spectrogram.png');
