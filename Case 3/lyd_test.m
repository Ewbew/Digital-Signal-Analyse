clc; close all; clear all;
%%

% Generate a composite audio signal, that is 5 seconds long with three
% different sine waves; f1 = 400 Hz, f2 = 785 Hz, f3 = 100 Hz.
% Then play it using soundsc:
fs = 44100; % Sampling frequency
t = 0:1/fs:10; % Time vector
f1 = 400; f2 = 785; f3 = 100; % Frequencies
compositeSignal = sin(2*pi*f1*t) + sin(2*pi*f2*t) + sin(2*pi*f3*t);
soundsc(compositeSignal,fs ); % Play the audi0