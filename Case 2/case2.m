%% Opgave 1 – Signal generation / kodning 
clear;
close all;
clc;
% A. Generer et lydsignal-array med "FSKgenerator" funktionen.

fstart = 2000; % transmission band frequency start
fend = 20000; % transmission band frequency end
Tsymbol = 0.5; % symbol duration in seconds
fs = 48000; % sampling frequency

sentence='hello world';
x = FSKgenerator(sentence, fstart, fend, Tsymbol, fs);
soundsc(x, fs)
%%
% B. Analyser signalet for at finde ud af, hvilke karakterer, som svarer
% til hvilke frekvenser. I skal se på signalet i både tids- og frekvens-domænet.%

figure(1)
N=length(x);    % antal sample 
plot((0:N-1)/fs,x)  %% Plot x signal i tide-domænet
xlim([1 1.009]);
xlabel("Tid [s]")
ylabel("Amplitude")
title("plot lydsignal tids-domænet")


figure(2)
plot(0:fs/N:fs-fs/N,abs(fft(x))) %% Plot x signal i frekvens-domænet
xlabel("Frekvens [Hz]")
ylabel("|X(f)|")
title("plot lydsignal frekvens-domænet")

figure(3)
plot(0:fs/N:fs-fs/N,abs(fft(x))) %% Plot x signal i frekvens-domænet
xlabel("Frekvens [Hz]")
ylabel("|X(f)|")
title("plot lydsignal frekvens-domænet")
xlim([500 10000])

h = FSKgenerator('h', fstart, fend, Tsymbol, fs);
e = FSKgenerator('e', fstart, fend, Tsymbol, fs);
l = FSKgenerator('l', fstart, fend, Tsymbol, fs);
o = FSKgenerator('o', fstart, fend, Tsymbol, fs);
space = FSKgenerator(' ', fstart, fend, Tsymbol, fs);

w = FSKgenerator('w', fstart, fend, Tsymbol, fs);
r = FSKgenerator('r', fstart, fend, Tsymbol, fs);
d = FSKgenerator('d', fstart, fend, Tsymbol, fs);

figure(4)
N =length(h);
x_akse = (0:fs/N:fs-fs/N);
plot(x_akse,abs(fft(h)))
hold on

N =length(e);
x_akse = (0:fs/N:fs-fs/N);
plot(x_akse,abs(fft(e)))

N =length(l);
x_akse = (0:fs/N:fs-fs/N);
plot(x_akse,abs(fft(l)))

N =length(o);
x_akse = (0:fs/N:fs-fs/N);
plot(x_akse,abs(fft(o)))

N =length(space);
x_akse = (0:fs/N:fs-fs/N);
plot(x_akse ,abs(fft(space)))

N =length(w);
x_akse = (0:fs/N:fs-fs/N);
plot(x_akse ,abs(fft(w)))

N =length(r);
x_akse = (0:fs/N:fs-fs/N);
plot(x_akse,abs(fft(r)))

N =length(d);
x_akse = (0:fs/N:fs-fs/N);
plot(x_akse ,abs(fft(d)))

xlabel("Frekvens [Hz]")
ylabel("|X(f)|")
title("plot lydsignal frekvens-domænet")
legend(["h","e","l", "o"," ", "w", "r", "d"])
xlim([1000 15000])
hold off

%%
%C. Analyser signalet vha. Short-Time Fourier Transform (kan læses om i bogen) – dvs. med 
% spektrogram-plot. Forklar trade-off imellem opløsningen i tid og frekvens.  %

figure(5)
spectrogram(x, hanning(512), 500, 1024, fs, 'yaxis')
ylim([0.5 10]);

figure(6)
spectrogram(x, hanning(2048), 500, 1024, fs, 'yaxis')
ylim([0.5 10]);
%%
%D. Eksperimenter med "FSKgenerator” funktionen for at få en forståelse af input
% parametrene.  
%% lydsignal ved Mindre båndbredde
figure(7)  
fstart = 1200; % transmission band frequency start
fend = 6000; % transmission band frequency end
Tsymbol = 0.5; % symbol duration in seconds
fs = 48000; % sampling frequency

sentence='hello world';
x = FSKgenerator(sentence, fstart, fend, Tsymbol, fs);
N=length(x);
plot(0:fs/N:fs-fs/N,abs(fft(x))) %% Plot x signal i frekvens-domænet
xlim([0 8000]);
xlabel("Frekvens [Hz]")
ylabel("|X(f)|")
title("frekvens-domænet til lydsignal ved Mindre båndbredde")
sound(x,fs)


%% lydsignal ved støre Tsymbol
figure(8)  
fstart = 2000; % transmission band frequency start
fend = 20000; % transmission band frequency end
Tsymbol = 0.8; % symbol duration in seconds
fs = 48000; % sampling frequency

sentence='hello world';
x = FSKgenerator(sentence, fstart, fend, Tsymbol, fs);
N=length(x);
plot(0:fs/N:fs-fs/N,abs(fft(x))) %% Plot x signal i frekvens-domænet
xlim([0 12000]);
xlabel("Frekvens [Hz]")
ylabel("|X(f)|")
title("frekvens-domænet til lydsignal ved støre Tsymbol")
sound(x,fs)

%% lydsignal ved mindre sampling frequency
figure(9) 
fstart = 500; % transmission band frequency start
fend = 1500; % transmission band frequency end
Tsymbol = 0.5; % symbol duration in seconds
fs = 4000; % sampling frequency

sentence='hello world';
x = FSKgenerator(sentence, fstart, fend, Tsymbol, fs);
N=length(x);
plot(0:fs/N:fs-fs/N,abs(fft(x))) %% Plot x signal i frekvens-domænet
xlim([0 2000]);
xlabel("Frekvens [Hz]")
ylabel("|X(f)|")
title("frekvens-domænet til lydsignal ved mindre sampling frequency")

sound(x,fs)

%% Opgave 2 – Dekodning 

% A. Send besked til en anden gruppe
% OBS: Vi sender en besked til os selv; vi afspiller nedenstående indkodede
% besked og optager på en telefon

fs = 30000;
fstart = 100;
fend = 10000;

lydsignal = FSKgenerator('hello world', 100, 10000, 0.2, fs); 

soundsc(lydsignal,fs);

%%

% Load audio file from the current work space folder 
path = matlab.desktop.editor.getActiveFilename;
thisFolder = fileparts(path);
audioFile = fullfile(thisFolder, 'hello_world_0m.m4a');
[y, fs] = audioread(audioFile);

%%
% Play said audio file
soundsc(y,fs);

%%
% Tilklip filen, så de stille perioder før og efter lydsignalet i
% optagelsen er fjernet
threshold = 0.01; 

mask = abs(y) > threshold; 

firstSample = find(mask, 1, 'first');
lastSample = find(mask, 1, 'last');

yTrimmed = y(firstSample:lastSample,:);
soundsc(yTrimmed,fs);

%% Visually inspecting the signal via spectrogram
figure(43)
spectrogram(yTrimmed, blackman(1000), 0,1000, fs, 'yaxis')

%% Decoding the signal using self-made algorithm
% Adding path to the current work space, so Matlab can access the function
addpath(thisFolder);

% received_message = FSKdecoder(signal, fs)
% mysymbolseq, fstart, fend, Tsymbol, fs

received_message = FSKdecoder(yTrimmed,100, 10000, 0.2, fs)
% Current version is quite slow, because it's O(N^2) complexity – could be
% cut down to 256, because those are the frequencies that we are interested
% in
 


%%
N = length(y);
t = (0:N-1)/fs;

figure(42)
plot(t,y)




%%


% Decode the received signal using a simple thresholding method
threshold = 0.5; % Define a threshold for detection
decodedChars = ''; % Initialize the decoded characters string
for i = 1:length(x)
    if x(i) > threshold
        decodedChars = [decodedChars, '1']; % Detected signal
    else
        decodedChars = [decodedChars, '0']; % No signal
    end
end


%% Opgave 3 – SNR-analyse af FSK-signal ved forskellige afstande

clear; close all; clc;

files = ["hello_world_0m.m4a", "hello_world_1m.m4a", "hello_world_2m.m4a", "hello_world_3m.m4a", "hello_world_4m.m4a"];
dists = [0 1 2 3 4];    % afstande i meter

fstart = 100; 
fend   = 10000;
Tsymbol = 0.2;  % symboltid (samme som ved optagelse)

SNRdB = nan(size(files));

figure; clf; hold on; grid on;
for k = 1:numel(files)
    % --- indlæs og trim stille perioder (samme som i opg. 2) ---
    path = matlab.desktop.editor.getActiveFilename;
    thisFolder = fileparts(path);
    audioFile = fullfile(thisFolder, files(k));
    [y, fs] = audioread(audioFile);

    threshold = 0.01;
    mask = abs(y) > threshold;
    yTrimmed = y(find(mask,1,'first'):find(mask,1,'last'));
    yTrimmed = yTrimmed - mean(yTrimmed);   % fjern DC-komponent

    % --- tag ét symbol til analyse ---
    N   = round(Tsymbol * fs);
    seg = yTrimmed(1:N);

    % --- FFT og power-spektret ---
    Y = fft(seg);
    P = abs(Y).^2 / (N*N);           % lineær effekt
    f = (0:N-1)*(fs/N);              % frekvensakse [Hz]

    % --- vælg kun FSK-båndet ---
    band = (f >= fstart) & (f <= fend);
    fB = f(band);
    PB = P(band);

    % --- beregn SNR (peak vs. median noise) ---
    [Ppk, ipk] = max(PB);
    excl = false(size(PB));
    excl(max(ipk-1,1):min(ipk+1,numel(PB))) = true;
    noise_lin = median(PB(~excl));
    SNRdB(k)  = 10*log10(Ppk / max(noise_lin, eps));

    % --- konverter til dB ---
    PdB = 10*log10(PB + eps);

    % --- plot i dB ---
    plot(fB, PdB, 'LineWidth', 1.2, 'DisplayName', ...
         sprintf('%s (SNR≈%.1f dB)', files(k), SNRdB(k)));
end

xlabel('Frekvens [Hz]');
ylabel('Effekt [dB]');
title('Power Spectrum (dB) for FSK-signal ved forskellige afstande');
legend('Location','best');
xlim([fstart fend]);


%% Opgave 3C – SNR som funktion af afstand 


% Til sidst plotter vi de beregnede SNR-værdier for hver afstand.
% Vi forventer, at SNR falder, når afstanden øges, da signalstyrken
% aftager hurtigere end støjniveauet ændrer sig.


figure; clf;
plot(dists, SNRdB, 'o-','LineWidth',1.5);
grid on;
xlabel('Afstand [m]');
ylabel('SNR [dB]');
title('SNR (peak/median) vs. afstand');




%% test for at tjekke decoder virker
clear; close all; clc;

fs = 30000;
fstart = 100;
fend = 10000;
Tsymbol = 0.2;

testmsg = 'hello';
x = FSKgenerator(testmsg, fstart, fend, Tsymbol, fs);
yhat = FSKdecoder(x, fstart, fend, Tsymbol, fs);

disp("Expected: " + testmsg);
disp("Decoded:  " + yhat);



%% Opgave 4A – Bit rate vs. fejlfri dekodning
clear; close all; clc;

% A) Opsæt test: sweep Tsymbol og mål fejlrate + bitrate

% Vi undersøger, hvor kort Tsymbol kan være, før dekoderen begynder at lave fejl.
% Antagelse: 1 symbol = 1 byte (ASCII) => 8 bits per symbol.
% Bitrate = 8 / Tsymbol  [bits/s]
%
% Parametre holdes som i opg. 2 for konsistens.
fs      = 30000;       % samplingsfrekvens
fstart  = 100;         % start på FSK-bånd
fend    = 10000;       % slut på FSK-bånd
msg     = 'hello world hello world';   % lidt længere teststreng
Ts_list = [0.50 0.30 0.20 0.15 0.10 0.08 0.06 0.05 0.04 0.03 0.025 0.020 0.010 0.005 0.002 0.0019];



bitrate = 8 ./ Ts_list;
FER     = nan(size(Ts_list));   % Frame Error Rate (fejl i hele tekst-strengen)
BER     = nan(size(Ts_list));   % Bit Error Rate (groft estimeret som char-fejl*8 / (len*8))

% Hjælpere til BER-estimat (pr. tegn): vi tæller tegn der er forkerte
count_char_errors = @(a,b) sum(a ~= b);

for k = 1:numel(Ts_list)
    Tsymbol = Ts_list(k);

    % --- Generér FSK-signal for hele beskeden
    x = FSKgenerator(msg, fstart, fend, Tsymbol, fs);

    % --- Dekod med vores egen decoder (samme som i opgave 2, men kun med 256 frekvenser)
    yhat = FSKdecoder(x, fstart, fend, Tsymbol, fs);

    % --- Fejlmåling
    % Frame Error Rate: 1 hvis bare ét tegn er forkert, ellers 0
    FER(k) = any(yhat ~= msg);

    % Bit Error Rate (grov): antal forkerte tegn * 8 / total bits
    n_char_err = count_char_errors(yhat, msg);
    BER(k) = (n_char_err * 8) / (length(msg) * 8);
end


%% B) Plot: bitrate vs. fejlrate (FER & BER)

figure;
subplot(2,1,1);
stem(bitrate, FER, 'filled'); grid on;
xlabel('Bitrate [bits/s]'); ylabel('FER');
title('Frame Error Rate vs. Bitrate (1 symbol = 1 byte)');
ylim([-0.05 1.05]);

subplot(2,1,2);
plot(bitrate, BER, 'o-','LineWidth',1.2); grid on;
xlabel('Bitrate [bits/s]'); ylabel('BER');
title('Bit Error Rate vs. Bitrate');

% Find maksimal fejlfri bitrate:
ok_idx = find(FER==0);
if ~isempty(ok_idx)
    max_ok_bitrate = max(bitrate(ok_idx));
    fprintf('Maks. fejlfri bitrate ≈ %.1f bits/s (Tsymbol ≈ %.3f s)\n', ...
            max_ok_bitrate, 8/max_ok_bitrate);
else
    fprintf('Ingen fejlfri bitrate i dette sweep – forøg Tsymbol eller SNR.\n');
end


%% C) Kort diskussion

% Vi ser, at når Tsymbol bliver kort (høj symbolrate/bitrate), begynder
% dekoderen at lave fejl. Det hænger sammen med DFT/vindue-teori:
% Færre samples pr. symbol => bredere hovedlobe (= ringere frekvensopløsning)
% og mere lækage, hvilket gør peaks mindre adskilte ved samme frekvensafstand.
% Dermed stiger sandsynligheden for, at nabofrekvenser forveksles – især ved
% begrænset SNR. (Se slides om vinduer/leakage/oploesning). 


%% Opgave 4B – Betydning af vindueslængde N for amplitude-spektret
clear; close all; clc;

fs = 30000;              % samplingsfrekvens (samme som i casen)
f0 = 4000;               % en tone inde i jeres bånd (kun ét symbol)
Ts_list = [0.50 0.20 0.10 0.05 0.02];   % symboltider => forskellige N
colors = lines(numel(Ts_list));

figure; hold on;
for k = 1:numel(Ts_list)
    Ts = Ts_list(k);
    N  = round(Ts*fs);               % vindueslængde (antal samples)
    n  = 0:N-1;
    x  = cos(2*pi*f0*n/fs);          % ren sinus (rektangulært vindue)

    % FFT (zero-padding for pæn kurve; opløsningen bestemmes stadig af N)
    Nfft = 2^nextpow2(max(4096, N));
    X = fft(x, Nfft);
    f = (0:Nfft-1)*(fs/Nfft);

    % Amplitude (normaliser peak til 0 dB for sammenligning)
    A = abs(X);
    A = A / max(A);
    AdB = 20*log10(A + 1e-12);

    plot(f(1:Nfft/2), AdB(1:Nfft/2), 'LineWidth', 1.2, 'Color', colors(k,:));

    % Teoretisk hovedlobebredde (rektangulært vindue): ca. 2*fs/N (nul-til-nul)
    dF_mainlobe = 2*fs/N;
    text(f0+200 + 150*k, -6-2*k, sprintf('N=%d  \\Delta f_{HL}≈ %.1f Hz', N, dF_mainlobe), ...
         'Color', colors(k,:), 'FontSize', 9);
end
grid on; xlim([0 9000]); ylim([-80 3]);
xlabel('Frekvens [Hz]'); ylabel('Amplitude [dB]');
title('Amplitude-spektre for forskellige N (rektangulært vindue)');
legend(arrayfun(@(Ts)sprintf('Ts=%.02f s (N=%d)', Ts, round(Ts*fs)), Ts_list, 'uni',0), 'Location','southwest');

%% Vindue-tradeoff: sammenlign rektangel vs. Hann ved fast N
figure; hold on;
Ts = 0.10; N = round(Ts*fs);
n  = 0:N-1;
x_rect = cos(2*pi*f0*n/fs);
x_hann = x_rect .* hann(N).';

Nfft = 2^nextpow2(max(4096, N));
f = (0:Nfft-1)*(fs/Nfft);

Xr = fft(x_rect, Nfft);  Xh = fft(x_hann, Nfft);
Ar = abs(Xr)/max(abs(Xr)); Ah = abs(Xh)/max(abs(Xh));
plot(f(1:Nfft/2), 20*log10(Ar(1:Nfft/2)+1e-12), 'LineWidth',1.2);
plot(f(1:Nfft/2), 20*log10(Ah(1:Nfft/2)+1e-12), 'LineWidth',1.2);
grid on; xlim([0 9000]); ylim([-100 3]);
xlabel('Frekvens [Hz]'); ylabel('Amplitude [dB]');
title(sprintf('Vindue-effekt ved N=%d (Ts=%.02f s), f_0=%g Hz', N, Ts, f0));
legend('Rektangel (smal HL, høje sidelobes)','Hann (bredere HL, lavere sidelobes)', 'Location','southwest');

%% Lille tabel til rapporten 
fprintf('\nTeoretisk hovedlobebredde (rektangulært vindue, nul-til-nul): Δf ≈ 2*fs/N\n');
for Ts = Ts_list
    N = round(Ts*fs);
    fprintf('Ts=%.3f s, N=%5d  ->  Δf≈ %.1f Hz\n', Ts, N, 2*fs/N);
end


%% Opgave 4C – SNR'ens betydning for systemets ydeevne
clear; close all; clc;

fs     = 30000;
fstart = 9000;
fend   = 10000;                % Bruger smallere bånd
Ts     = 0.01;                 % og kortere symboltid for at vise SNR betydning
msg    = 'hello world hello world';

% Liste af SNR-værdier i dB
SNR_dB = 0:5:60;               % fra 0 dB (meget støj) til 60 dB (næsten perfekt)
FER = zeros(size(SNR_dB));
BER = zeros(size(SNR_dB));

% Original FSK-signal
x = FSKgenerator(msg, fstart, fend, Ts, fs);
Px = mean(x.^2);               % signalets gennemsnitlige effekt

for k = 1:length(SNR_dB)
    snr_val = SNR_dB(k);
    
    % Beregn nødvendig støj-effekt for den ønskede SNR
    Pn = Px / (10^(snr_val/10));
    noise = sqrt(Pn) * randn(size(x));     % hvidt støjsignal
    
    % Tilføj støj
    y = x + noise;

    % Decode
    yhat = FSKdecoder(y, fstart, fend, Ts, fs);
    
    % Frame- og bitfejl
    FER(k) = any(yhat ~= msg);             % 1 hvis mindst én fejl
    BER(k) = sum(yhat ~= msg) / length(msg);
    
    fprintf('SNR = %2d dB | FER = %.2f | BER = %.3f\n', snr_val, FER(k), BER(k));
end

% Plot
figure;
subplot(2,1,1);
semilogy(SNR_dB, BER, 'o-','LineWidth',1.2);
xlabel('SNR [dB]'); ylabel('BER');
title('Bit Error Rate som funktion af SNR');
grid on; ylim([1e-4 1]);

subplot(2,1,2);
stem(SNR_dB, FER,'filled');
xlabel('SNR [dB]'); ylabel('FER');
title('Frame Error Rate som funktion af SNR');
grid on; ylim([-0.05 1.05]);
