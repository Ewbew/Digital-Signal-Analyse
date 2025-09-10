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

var_belastet = var(x_belastet);
var_ubelastet = var(x_ubelastet);

std_belastet = std(x_belastet);
std_ubelastet = std(x_ubelastet);

fprintf('Belastet (1 kg):   mean = %.6g, var = %.6g, std = %.6g\n', mean_belastet, var_belastet, std_belastet);
fprintf('Ubelastet (0 kg):  mean = %.6g, var = %.6g, std = %.6g\n\n', mean_ubelastet, var_ubelastet, std_ubelastet);



%% Opg 1 – Histogrammer 
% Histogrammer med ±1σ, ±2σ, ±3σ linjer
figure;

% ---------- Ubelastet ----------
subplot(1,2,1)
histogram(x_ubelastet, 'Normalization', 'pdf'); hold on; grid on
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
histogram(x_belastet, 'Normalization', 'pdf'); hold on; grid on
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
%x_fft = detrend(x_ubelastet);
x_fft = x_ubelastet - mean_ubelastet;

% Beregn FFT
X = fft(x_fft);


% Frekvensakse (0 ... fs)
f = (0:Nf-1)*(fs/Nf);

% Plot kun den positive halvdel af spektret
figure;
plot(f(1:floor(Nf/2)), abs(X(1:floor(Nf/2)))/Nf, 'LineWidth', 1.2);
xlabel('Frekvens [Hz]');
ylabel('Amplitude');
title('Frekvensspektrum af ubelastet signal (med detrend)');
grid on;

%Spektret af det detrendede ubelastede signal fremstår bredbåndet og uden dominerende frekvenskomponenter, hvilket er karakteristisk for hvid støj.
%Selvom amplituden i spektret ikke er nul, er det netop forventeligt, da støjens energi er fordelt tilfældigt og jævnt over hele frekvensbåndet.
%Dette bekræfter, at den pålagte støj i signalet tilnærmelsesvist opfører sig som hvid støj.

%% Opg 1: Forskel i bit-niveau i gram


delta_counts = mean_ubelastet - mean_belastet;   % = 298.44
LSB_g = 1000 / delta_counts;                     % = 3.3514 g/count
LSB_kg = 1 / delta_counts;                       % = 0.003351 kg/count

fprintf('LSB = %.4f g/count (%.6f kg/count)\n', LSB_g, LSB_kg);


%% Opgave 2 - MA-filter



% Create a new figure
figure('Name', 'Moving Average Filter Results', 'NumberTitle', 'off');

% Call the function for the first row (Belastet)
M = 10;
y_belastet = moving_average_filter(x_belastet, M);

subplot(2,3,1);
plot(x_belastet,'*')
hold on
plot(y_belastet,'*')
legend('Before filter', 'after filter')
title('M = 10');
xlim([0, 600]); % Set xlim for Belastet
ylim([800, 1600]); % Set ylim for all subplots

M = 50;
y_belastet = moving_average_filter(x_belastet, M);

subplot(2,3,2);
plot(x_belastet,'*')
hold on
plot(y_belastet,'*')
legend('Before filter', 'after filter')
title('M = 50');
xlim([0, 600]); % Set xlim for Belastet
ylim([800, 1600]); % Set ylim for all subplots

M = 100;
y_belastet = moving_average_filter(x_belastet, M);

subplot(2,3,3);
plot(x_belastet,'*')
hold on
plot(y_belastet,'*')
legend('Before filter', 'after filter')
title('M = 100');
xlim([0, 600]); % Set xlim for Belastet
ylim([800, 1600]); % Set ylim for all subplots

% Call the function for the second row (Ubelastet)
M = 10;
y_ubelastet = moving_average_filter(x_ubelastet, M);

subplot(2,3,4);
plot(x_ubelastet,'*')
hold on
plot(y_ubelastet,'*')
legend('Before filter', 'after filter')
title('M = 10');
xlim([0, 800]); % Set xlim for Ubelastet
ylim([800, 1600]); % Set ylim for all subplots

M = 50;
y_ubelastet = moving_average_filter(x_ubelastet, M);

subplot(2,3,5);
plot(x_ubelastet,'*')
hold on
plot(y_ubelastet,'*')
legend('Before filter', 'after filter')
title('M = 50');
xlim([0, 800]); % Set xlim for Ubelastet
ylim([800, 1600]); % Set ylim for all subplots

M = 100;
y_ubelastet = moving_average_filter(x_ubelastet, M);

subplot(2,3,6);
plot(x_ubelastet,'*')
hold on
plot(y_ubelastet,'*')
legend('Before filter', 'after filter')
title('M = 100');
xlim([0, 800]); % Set xlim for Ubelastet
ylim([800, 1600]); % Set ylim for all subplots

% Add titles for the rows
annotation('textbox', [0.03, 0.68, 0.8, 0.1], 'String', 'Belastet', ...
           'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
           'HorizontalAlignment', 'left');

annotation('textbox', [0.03, 0.18, 0.8, 0.1], 'String', 'Ubelastet', ...
           'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
           'HorizontalAlignment', 'left');


% % Create a new figure
% figure('Name', 'Moving Average Filter Results', 'NumberTitle', 'off');
% 
% % Define M values and corresponding x data
% M_values = [10, 50, 100];
% x_data = {x_belastet, x_ubelastet};
% titles = {'Belastet', 'Ubelastet'};
% 
% % Loop through each dataset (Belastet and Ubelastet)
% for j = 1:2
%     for i = 1:length(M_values)
%         M = M_values(i);
%         y_filtered = moving_average_filter(x_data{j}, M);
% 
%         subplot(2, 3, (j-1)*3 + i);
%         plot(x_data{j}, '*');
%         hold on;
%         plot(y_filtered, '*');
%         legend('Before filter', 'After filter');
%         title(['M = ', num2str(M)]);
% 
%         % Set axis limits
%         if j == 1
%             xlim([0, 600]); % For Belastet
%         else
%             xlim([0, 800]); % For Ubelastet
%         end
%         ylim([800, 1600]); % For all subplots
%     end
% end
% 
% % Add titles for the rows
% annotation('textbox', [0.1, 0.85, 0.8, 0.1], 'String', 'Belastet', ...
%            'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
%            'HorizontalAlignment', 'center');
% 
% annotation('textbox', [0.1, 0.45, 0.8, 0.1], 'String', 'Ubelastet', ...
%            'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
%            'HorizontalAlignment', 'center');


%% Plot af histogrammer og varians/støj- effekten
y_belastet_M10 = moving_average_filter(x_belastet,50);
y_belastet_M50 = moving_average_filter(x_belastet, 50);
y_belastet_M100 = moving_average_filter(x_belastet, 100);
y_ubelastet_M10 = moving_average_filter(x_ubelastet, 10);
y_ubelastet_M50 = moving_average_filter(x_ubelastet, 50);
y_ubelastet_M100 = moving_average_filter(x_ubelastet, 100);

% Example filtered signals (replace these with your actual signals)
filtered_signals = {y_belastet_M10, y_belastet_M50, y_belastet_M100, ...
                    y_ubelastet_M10, y_ubelastet_M50, y_ubelastet_M100};

% Create a new figure for histograms
figure('Name', 'Histograms of MA Filtered Signals', 'NumberTitle', 'off');

% Loop through each filtered signal to plot histograms and calculate variance
for i = 1:length(filtered_signals)
    subplot(2, 3, i);
    histogram(filtered_signals{i}, 'FaceColor', 'b', 'EdgeColor', 'k');
    title(['Signal ', num2str(i)]);
    xlabel('Value');
    ylabel('Frequency');

    % Calculate variance
    variance_value = var(filtered_signals{i});

    % Print variance to command window
    fprintf('Variance of Signal %d: %.4f\n', i, variance_value);
end

% Adjust layout
sgtitle('Histograms and Variance of MA Filtered Signals');
