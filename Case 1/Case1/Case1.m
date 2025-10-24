%%CASE 1 - Indlæsning og opsætning af variable

clear; clc; close all;

S = load('vejecelle_data.mat');
x = S.vejecelle_data;

fs = 300; %Hz
N = length(x);
n = (0:N-1);
t = n/fs;

% tester


%% Opdeling af signal hhv belastet og ubelastet
plot(0:N-1,x);
%Ovenstående plot af rå data viser, at det transiente forløb stopper efter
%ca. 150 samples og overgangen fra belastet sker ved ca. 1000-1100 samples.
%Derfor tager vi forløbet med belastet vejecelle fra 150-1000 samples og
%ubelastet fra 1100-N
x_belastet = x(150:1000);
x_ubelastet = x(1100:N-1);

%% Opg 1 - Data analyse

mean_belastet = mean(x_belastet);
mean_ubelastet = mean(x_ubelastet);

var_belastet = var(x_belastet,1);
var_ubelastet = var(x_ubelastet,1);

std_belastet = std(x_belastet,1);
std_ubelastet = std(x_ubelastet,1);

fprintf('Belastet (1 kg):   mean = %.6g, var = %.6g, std = %.6g\n', mean_belastet, var_belastet, std_belastet);
fprintf('Ubelastet (0 kg):  mean = %.6g, var = %.6g, std = %.6g\n\n', mean_ubelastet, var_ubelastet, std_ubelastet);



%% Opg 1 – Histogrammer 
% Histogrammer med ±1σ, ±2σ, ±3σ linjer
figure;

% ---------- Ubelastet ----------
subplot(1,2,1)
histogram(x_ubelastet, 'Normalization', 'pdf',NumBins=50); hold on; grid on
title('Ubelastet');
xlabel('Ukendt enhed'); ylabel('Sandsynlighedstæthed');

% Lodrette linjer ved ±1σ, ±2σ, ±3σ
mu = mean_ubelastet;
sigma = std_ubelastet;
xline(mu, 'k-', 'LineWidth', 1.2, 'DisplayName','mean');
xline(mu + sigma, 'r--', 'LineWidth', 1.2, 'DisplayName','±1σ');
xline(mu - sigma, 'r--', 'LineWidth', 1.2, HandleVisibility='off');
xline(mu + 2*sigma, 'g--', 'LineWidth', 1.2, 'DisplayName','±2σ');
xline(mu - 2*sigma, 'g--', 'LineWidth', 1.2, HandleVisibility='off');
xline(mu + 3*sigma, 'b--', 'LineWidth', 1.2, 'DisplayName','±3σ');
xline(mu - 3*sigma, 'b--', 'LineWidth', 1.2, HandleVisibility='off');
legend('show');

% ---------- Belastet ----------
subplot(1,2,2)
histogram(x_belastet, 'Normalization', 'pdf',NumBins=50); hold on; grid on
title('Belastet (1 kg)');
xlabel('Ukendt enhed'); ylabel('Sandsynlighedstæthed');

mu = mean_belastet;
sigma = std_belastet;
xline(mu, 'k-', 'LineWidth', 1.2, 'DisplayName','mean');
xline(mu + sigma, 'r--', 'LineWidth', 1.2, 'DisplayName','±1σ');
xline(mu - sigma, 'r--', 'LineWidth', 1.2, HandleVisibility='off');
xline(mu + 2*sigma, 'g--', 'LineWidth', 1.2, 'DisplayName','±2σ');
xline(mu - 2*sigma, 'g--', 'LineWidth', 1.2, HandleVisibility='off');
xline(mu + 3*sigma, 'b--', 'LineWidth', 1.2, 'DisplayName','±3σ');
xline(mu - 3*sigma, 'b--', 'LineWidth', 1.2, HandleVisibility='off');
legend('show');

%Histogrammerne for både den belastede og ubelastede tilstand viser en tydelig klokkeform, hvilket tyder på, at målesignalet (og støjen) er tilnærmelsesvist normalfordelt.

%Den målte spredning (standardafvigelse) stemmer også godt overens med histogrammets bredde:
%Flertallet af samples ligger indenfor ca. ±1 til ±2 standardafvigelser fra middelværdien, som forventet for en normalfordeling.

%% Opg 1 - PLot af frekvensspektret

% FFT af ubelastet signal – brug detrend til at fjerne offset/trend
Nf = length(x_ubelastet);             % antal samples

% Fjern både DC-offset og evt. langsom drift (lineær trend)
x_fft_u = detrend(x_ubelastet);
x_fft_b = detrend(x_belastet);

% Beregn FFT
X_u = fft(x_fft_u);
X_b = fft(x_fft_b);

% Frekvensakse (0 ... fs)
f = (0:Nf-1)*(fs/Nf);

% Plot kun den positive halvdel af spektret
figure;
subplot(2,1,1);
plot(f(1:floor(Nf/2)), abs(X_u(1:floor(Nf/2)))/Nf, 'LineWidth', 1.2);
xlabel('Frekvens [Hz]');
ylabel('Amplitude');
title('Frekvensspektrum af ubelastet signal (med detrend)');
grid on;

% Plot only the positive half of the spectrum for the loaded signal
subplot(2,1,2);
plot(f(1:floor(Nf/2)), abs(X_b(1:floor(Nf/2)))/Nf, 'LineWidth', 1.2);
xlabel('Frekvens [Hz]');
ylabel('Amplitude');
title('Frekvensspektrum af belastet signal (med detrend)');
grid on;

%Spektret af det detrendede ubelastede signal fremstår bredbåndet og uden dominerende frekvenskomponenter, hvilket er karakteristisk for hvid støj.
%Selvom amplituden i spektret ikke er nul, er det netop forventeligt, da støjens energi er fordelt tilfældigt og jævnt over hele frekvensbåndet.
%Dette bekræfter, at den pålagte støj i signalet tilnærmelsesvist opfører sig som hvid støj.

%% Opg 1: Forskel i bit-niveau i gram


delta_counts = mean_ubelastet - mean_belastet;   
LSB_g = 1000 / delta_counts;                     
LSB_kg = 1 / delta_counts;                       
fprintf('LSB = %.4f g/count (%.6f kg/count)\n\n', LSB_g, LSB_kg);



%% OPGAVE 2: MA-filter

% Variansanalyse af MA-filtrering
M_vals = [10, 50, 100];

% Variansanalyse af MA-filtrering
for M = M_vals
    y_b = moving_average_filter(x_belastet, M);
    y_u = moving_average_filter(x_ubelastet, M);

    y_b = y_b(M:end); y_u = y_u(M:end);
    var_b = var(y_b); var_u = var(y_u);

    fprintf('MA M=%d | Var belastet: %.4f | Var ubelastet: %.4f\n', M, var_b, var_u);
end

% Histogrammer for filtrerede signaler
figure;
for i = 1:3
    M = M_vals(i);
    y_u = moving_average_filter(x_ubelastet, M);
    subplot(2,3,i); histogram(y_u,400); title(sprintf('MA M=%d',M));
    axis([1300, 1500,0,600])
    y_b = moving_average_filter(x_belastet, M);
    subplot(2,3,i+3); histogram(y_b,400); title(sprintf('MA M=%d',M));
    axis([1000, 1200,0,300])
end

%Add titles for the rows
annotation('textbox', [0.048, 0.68, 0.8, 0.1], 'String', 'Ubelastet', ...
           'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
           'HorizontalAlignment', 'left');

annotation('textbox', [0.05, 0.18, 0.8, 0.1], 'String', 'Belastet', ...
           'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
           'HorizontalAlignment', 'left');



%% Moving-average (midlingsfilter) med filterlængde 100
M = 100;                          % filterlængde
N = numel(x);

xmem = zeros(M-1,1);             
y    = zeros(N,1);

for nn = 1:N
    % Gennemsnit af (M-1) tidligere + den aktuelle prøve
    y(nn) = (sum(xmem) + x(nn)) / M;

    xmem = [x(nn); xmem(1:end-1)];
end

% Drop opstartsdel (de første M-1 samples er delvist udfyldt)
y_100M = y(M:end);

figure;
plot(y_100M, 'LineWidth', 2);
title('data med midlingsfilter (M = 100)');
ylabel('kg');
xlabel('Samples');
grid on;



%% Beregn den maksimale af FIR-midlingsfilteret:

% Definer den ønskede maksimale indsvingningstid i sekunder 
max_indsvingningstid_sec = 0.1; % 100 millisekunder 

% Beregn det tilsvarende antal samples baseret på samplingsfrekvensen 
max_indsvingningstid_samples = max_indsvingningstid_sec * fs; 

% Beregn den maksimale længde af FIR-midlingsfilteret 
max_fir_length = floor(max_indsvingningstid_samples);

disp(['Den maksimale længde af FIR-midlingsfilteret: ', num2str(max_fir_length), ...
    ' Samples']); 




%% Implementér et eksponentielt midlingsfilter med α = 0.1
alpha = 0.1;
yold  = 0;

for n = 1:length(x)
    y(n) = alpha .* x(n) + (1 - alpha) .* yold;
    yold = y(n);
end

subplot(3,1,1);  % Create a subplot for the first filter
plot(y);  % plot filtreret med eksponentielt midlingsfilter (alpha = 0.1)
title('data med eksponentielt midlingsfilter alpha=0.1');
ylabel('kg');
xlabel('Samples');

% Implementér et eksponentielt midlingsfilter med α = 0.5
alpha = 0.5;
yold  = 0;

for n = 1:length(x)
    y(n) = alpha .* x(n) + (1 - alpha) .* yold;
    yold = y(n);
end

subplot(3,1,2);  % Create a subplot for the second filter
plot(y);  % plot filtreret med eksponentielt midlingsfilter (alpha = 0.5)
title('data med eksponentielt midlingsfilter alpha=0.5');
ylabel('kg');
xlabel('Samples');

% Implementér et eksponentielt midlingsfilter med α = 0.9
alpha = 0.9;
yold  = 0;

for n = 1:length(x)
    y(n) = alpha .* x(n) + (1 - alpha) .* yold;
    yold = y(n);
end

subplot(3,1,3);  % Create a subplot for the third filter
plot(y);  % plot filtreret med eksponentielt midlingsfilter (alpha = 0.9)
title('data med eksponentielt midlingsfilter alpha=0.9');
ylabel('kg');
xlabel('Samples');


%% α-værdien, således at får samme støj-reduktion, som for 100. ordens FIR midlingsfilter
alpha = 0.0198;
yold  = 0;

for n = 1:length(x)
    y(n) = alpha.* x(n) + (1 - alpha) .* yold;
    yold = y(n);
end

figure; clf
plot(y,    'LineWidth', 1.6); hold on; grid on
plot(y_100M,'LineWidth', 1.6);
title('Eksponentielt midlingsfilter (\alpha=0.0198) vs. Midlingsfilter (M=100)')
ylabel('kg'); xlabel('Samples');
legend('Eksponentielt midlingsfilter (\alpha=0.0198)','Midlingsfilter M=100','Location','best');
hold off
%%
