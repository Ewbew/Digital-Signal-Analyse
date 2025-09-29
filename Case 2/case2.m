%% Opgave 1 – Signal generation / kodning 
clear;
close all;
clc;
% a) Generer et lydsignal-array med "FSKgenerator" funktionen.

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

