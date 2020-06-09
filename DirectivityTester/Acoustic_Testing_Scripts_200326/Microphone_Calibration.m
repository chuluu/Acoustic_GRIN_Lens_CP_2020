%% Intro
clc
clear

load Sensitivity.mat;

disp(['To calibrate the microphone start with these steps: ']);
disp(['Check current Sensitivity: ', num2str(Sensitivity)]);
disp(['Check the dB setting, there should be a switch to adjust it, it should be 94 or 114']);
disp(['Press the On button, if it is quite or you must hold it, it probably needs a new battery.']);

disp(['Once the tone is constant, place Calibrate on microphone,']);
disp(['the microphone should insert into the hole']);

input(['Now hit enter on this program']);

disp(['Check the dB displayed in the command window: ']);

%%
fs = 51200;
T = 2;
devices = daq.getDevices;

devices(3);

s = daq.createSession('ni');
addAnalogInputChannel(s,'cDAQ1Mod1', 0, 'Microphone');

s.Channels.Sensitivity = Sensitivity;
% s.Channels.MaxSoundPressureLevel = 137;
s.Channels;

s.DurationInSeconds = T;
s.Rate = fs;

[y_recorded,t_recorded] = s.startForeground();

%% Process Data
win = hamming(length(y_recorded));
Windowed_signal = y_recorded.*win;
% plot(t_recorded,Windowed_signal);

%%
Max_Pressure = max(Windowed_signal);
dB_SPL_Sine = 20*log10(abs(Max_Pressure./(2*10^-5)))
disp(['If the value is not either 94, or 114, change the sensitivity']);
disp(['Typically this value is anywhere from 0.5 - 1']);
disp(['Increasing it will decrease the dB value typically']);
Sensitivity = input('New Sensitivity: ');
save('Sensitivity.mat','Sensitivity');