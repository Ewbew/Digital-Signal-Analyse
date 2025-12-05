
clc; close all; clear all;
%% Opgave 4 
%{
I skal analysere de optagne signaler i Matlab - dvs. plotte signalerne 
(ser de ud som forventet?).

Check frekvens-spektre etc..


Dernæst skal I benytte krydskorrelation til at finde tidsforskydningen.
OBS !! I må ikke benytte den indbyggede Matlab-funktion, den skal I selv implementere.
Forklar hvad I ser i plottet af krydskorrelationen.
Hvor præcist kan I måle afstand ?
Påvirkes systemet af andre signaler ? (fx. baggrundssnak / støj, refleksioner 
fra andre objekter,
direkte signalvej fra højttaler til mikrofon, ..)
%}





% Check signalerne – ser de ud som forventet?

% Get the path of the current script
path = matlab.desktop.editor.getActiveFilename;
thisFolder = fileparts(path);

% Define the folder containing .dat files
audioFolder = fullfile(thisFolder, 'Sonar_test_data_lyddodtrum', 'Sonar_test_data');

% Get a list of all .dat files in the folder
datFiles = dir(fullfile(audioFolder, '*.dat'));

% Sampling frequency (adjust if known)
Fs = 48000; % Example: 48 kHz, change to your actual sampling rate

for k = 1:length(datFiles)
    fileName = datFiles(k).name;
    filePath = fullfile(audioFolder, fileName);
    
    % Load the data
    data = load(filePath);
    
    % Time-domain plot
    figure;
    subplot(2,1,1);
    plot(data);
    [~, nameOnly, ~] = fileparts(fileName);
    title(['Time Domain: ', nameOnly], 'Interpreter', 'none');
    xlabel('Sample Index');
    ylabel('Amplitude');
    
    % Frequency-domain plot
    N = length(data);
    Y = fft(data);
    f = (0:N-1)*(Fs/N); % Frequency axis
    magnitude_dB = 20*log10(abs(Y)/N + eps); % Convert to dB, avoid log(0)
    
    subplot(2,1,2);
    plot(f(1:floor(N/2)), magnitude_dB(1:floor(N/2))); % Only positive frequencies
    title('Frequency Spectrum (Magnitude in dB)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
end



%%
% Original probe signal (chirp):
chirp = [32718, 32547, 32222, 31706, 30964, 29957, 28649, 27005, ...
    24995, 22595, 19788, 16569, 12947, 8946, 4609, 0, ...
    -4795, -9666, -14481, -19087, -23312, -26975, -29886, -31863, ...
    -32738, -32373, -30672, -27598, -23180, -17531, -10850, -3425, ...
    4370, 12093, 19250, 25330, 29842, 32356, 32555, 30274, ...
    25541, 18604, 9934, 214, -9704, -18868, -26312, -31164, ...
    -32764, -30761, -25202, -16569, -5774, 5919, 17018, 25997, ...
    31505, 32594, 28905, 20788, 9319, -3798, -16465, -26510, ...
    -32071, -31972, -26005, -15067, -1086, 13279, 25116, 31863, ...
    31891, 24952, 12354, -3212, -18194, -28974, -32764, -28378, ...
    -16673, -483, 15999, 28270, 32767, 27967, 14984, -2571,  ...
    -19509, -30572, -32099, -23322, -6826, 12093, 27058, 32767, ...
    26952, 11392, -8468, -25330, -32710, -27511, -11455, 9307, ...
    26439, 32758, 25339, 7022, -14481, -29736, -31737, -19261, ...
    2237, 22827, 32629, 26635, 7454, -15588, -30848, -30274, ...
    -13828, 10075, 28688, 31706, 17133, -7180, -27533, -32052];


%%
chirp = chirp(:);         % ensure column vector (samples × channels)
chirp = double(chirp);    % optional: use double for soundsc

soundsc(chirp/32768, 48000);
%%

filePath = fullfile(audioFolder, 'Sonar2m.dat');

% Load the data
data = load(filePath);

data = data';

%%

% Implement cross-correlation to find time delay
[r, lags] = crosscorr(data, chirp);
% Plot the crosscorrelation using r and lags
figure;
plot(lags, r);
title('Cross-Correlation between Data and Chirp Signal');
xlabel('Lags (samples)');
ylabel('Cross-Correlation Amplitude');
grid on;



