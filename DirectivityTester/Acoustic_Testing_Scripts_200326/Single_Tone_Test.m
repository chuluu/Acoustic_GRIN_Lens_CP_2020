clc
clear

%% Constants
load Sensitivity.mat;
Std_P = 2*10^-5;
V = [2.18,4.36,5.17,5.45]; % Voltage RMS Amplitude [0.5, 1, 1.5, 2]
Power = (V.^2)./(8); % Watts
fs = 51200;

%% Inputs
Fundamental_Freq = 1000;
T = 6;
Zeroes_Delay = 17000;
A = 0.5;

%% Sine Wave Generation
t_delay = 0:1/fs:(Zeroes_Delay/fs)+T; % gets timestep
t = 0:1/fs:T; % gets timestep

x = A.*sin(2*pi*Fundamental_Freq*t);
x_prime = [zeros(1,Zeroes_Delay),x]; % delay with zeros
sound(x_prime,fs);

%% Acquire Data
devices = daq.getDevices;

devices(3);

s = daq.createSession('ni');
addAnalogInputChannel(s,'cDAQ1Mod1', 0, 'Microphone');

s.Channels.Sensitivity = Sensitivity;
%s.Channels.MaxSoundPressureLevel = 137;
s.Channels

s.DurationInSeconds = T+0.5;
s.Rate = fs;

[y_recorded,t_recorded] = s.startForeground();

%% Process Data
win = hamming(length(y_recorded)); % Create Window
Windowed_signal = y_recorded.*win; % Window input signal
plot(t_recorded,y_recorded); hold on;
plot(t_recorded,Windowed_signal);

%% Take max amplitude
Max_Pressure = max(Windowed_signal);
dB_SPL_Sine = 20*log10(abs(Max_Pressure./(Std_P)))

%% Plot Full FFT
[Y,f] = MyFFT(Windowed_signal,fs);
figure
semilogx(f,20*log10(fftshift(abs(Y/fs))./Std_P)); hold on
title(['Single Tone Frequency Response Fundamental = ', num2str(Fundamental_Freq)]); 
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
xlim([100 10000])