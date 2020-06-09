%% Setup
%{ 
This program plots a polar pattern utilizing Farina Techniques and impulse windowing over
a polar sweep of 0 - 180 in 5 degree steps (37 data points)

Hardware:
cDAQ - 9174: Used to interface NI modules together
NI 9250: Audio Recording module
NI 9401: 100ns, TTL 3.3V Digital Input/Output
NI 9263: Not used currently

QSC: GX3 300-Watt Power Amplifier: used to amplify signal coming from UR12
Steinberg: UR12 USB Audio Interface: Processes signal, better sound card

CNC Sheild + 4X DRV8825
Arduino: UNO R3 Board: Board that triggers things

External Software:
Arduino Code that handles moving motor control and pulsing for data
collection

%}

clear 
clc 
close all

%% Variables to Input
load Sensitivity.mat;
Std_P = 2*10^-5;
Number_of_Recordings = 4; % Number of recordings
Start_index = 2; % Start index for personal code
fs = 51200; % Sampling Rate
f1 = 100; % Start Frequency
f2 = 10000; % Stop Frequency
T = 1; % Duration of Signal

%% Create Stimulus
Zeroes_Delay = 17000;
w1 = 2*pi*f1;
w2 = 2*pi*f2;
t_delay = 0:1/fs:(Zeroes_Delay/fs)+T;     % gets timestep
t = 0:1/fs:T;     % gets timestep
A = 1; 
x = A.*sin((w1*T)/log(w2/w1).*(exp((t./T)*log(w2/w1))-1));
x_prime = [zeros(1,Zeroes_Delay),x]; % delay with zeros

%% Inverse Filtering exp manipulation
x_inv = flip(x);
k = exp((t*log(f2/f1))./T);
z_t = x_inv./k;

%% Initiate Devices
devices = daq.getDevices;
devices(3);
Digital_Read = daq.createSession('ni');
addDigitalChannel(Digital_Read,'cDAQ1Mod2', 'Port0/Line0:1', 'InputOnly');

Microphone_Recording_Session = daq.createSession('ni');
addAnalogInputChannel(Microphone_Recording_Session,'cDAQ1Mod1', 0, 'Microphone');
Microphone_Recording_Session.Channels.Sensitivity = Sensitivity; % Make sure this is calibrated
%Microphone_Recording_Session.Channels.MaxSoundPressureLevel = 137;
Microphone_Recording_Session.Channels % Checks Channels
Microphone_Recording_Session.DurationInSeconds = T+0.5;
Microphone_Recording_Session.Rate = fs;


%% The Test
y_recorded = zeros(Number_of_Recordings,1);
t_recorded = zeros(Number_of_Recordings,1);
A = [0, 0.5, 1, 1.5]; %0.25*ones(1,Number_of_Recordings);
for a = 1:Number_of_Recordings
    while true
    data = Digital_Read.inputSingleScan();
    data(1) = 1;
    if data(1)==1;
        disp(['Index: ',num2str(a)]);    % Delay until Pulse Found
        break;
    end
    end
    x = A(a).*sin((w1*T)/log(w2/w1).*(exp((t./T)*log(w2/w1))-1));
    x_prime = [zeros(1,Zeroes_Delay),x]; % delay with zeros
    sound(x_prime,fs);
    [y_recorded_val,t_recorded_val] = Microphone_Recording_Session.startForeground();
    for b = 1:length(y_recorded_val)
        y_recorded(b,a) = y_recorded_val(b);
    end
end

%% Build LP Filter
[b_stop,a_stop] = butter(4, 0.5,'low');
num_bins = length(y_recorded(:,1));
H_resp = freqz(b_stop,a_stop, floor(num_bins/2));
% plot([0:1/(num_bins/2-1):1], abs(H_resp));

%% Post Processing - Windowing and LP Filtering

win = hamming(length(y_recorded(:,1)));
dB_SPL_Sine = zeros(1,Number_of_Recordings);
Windowed_signal = zeros(length(y_recorded(:,1)),Number_of_Recordings);
y_t_Filtered_Windowed = zeros(length(y_recorded(:,1)),Number_of_Recordings);

% Window and Filter the signal
for a = Start_index:Number_of_Recordings
    Windowed_signal(:,a) = y_recorded(:,a).*win;
    Max_Pressure = max(Windowed_signal(:,a));
    dB_SPL_Sine(a) = 20*log10(abs(Max_Pressure./(Std_P)));
    [y_t_Filtered_Windowed(:,a)] = filter(b_stop,a_stop,Windowed_signal(:,a));
end

%% Post Processing - Inverse Filter Convolution
h_t_array = zeros(Number_of_Recordings,1);
for a = Start_index:Number_of_Recordings
    [h_t_val,time_h_t] = MyFFTConv(y_t_Filtered_Windowed(:,a)',z_t,fs);
    for b = 1:length(h_t_val)
        h_t_array(b,a) = h_t_val(b);
    end
end

%% Plotting Section
% Plot Spectrogram Farinas
fig1 = figure(1);
% Input
subplot(2,1,1); spectrogram(x,hamming(1600),1600/2,2048,fs,'yaxis'); 
set(gca,'YScale','log')
title('Farina Chirp Input');
% Recorded
subplot(2,1,2); spectrogram(Windowed_signal(:,2),hamming(1600),1600/2,2048,fs,'yaxis'); 
set(gca,'YScale','log')
title('Farina Chirp Output');


% Plot Recordings Filtered and Windowed
fig2 = figure(2);
for a = Start_index:Number_of_Recordings
    subplot(3,1,1); plot(t_recorded_val,y_recorded(:,a)); hold on;
end
xlabel('Time'); ylabel('Magnitude'); title('Recording');
for a = Start_index:Number_of_Recordings
   subplot(3,1,2);  plot(t_recorded_val,Windowed_signal(:,a)); hold on;
end
xlabel('Time'); ylabel('Magnitude'); title('Recording Windowed');

for a = Start_index:Number_of_Recordings
    subplot(3,1,3); plot(t_recorded_val,y_t_Filtered_Windowed(:,a)); hold on;
end
xlabel('Time'); ylabel('Magnitude'); title('Recording Windowed and Filterd');

% Plot Impulse Response
fig7 = figure(3);
% Spectrogram Plot Impulse
subplot(2,1,1); spectrogram(h_t_array(:,2),hamming(800),800/2,2048,fs,'yaxis'); 
set(gca,'YScale','log')
title('Farina Chirp Impulse Response');
% Plot Impulse
for a = Start_index:Number_of_Recordings
    subplot(2,1,2); plot(time_h_t,h_t_array(:,a)); hold on;
end
title('h(t) impulse response'); xlabel('Time'); ylabel('Magnitude');
ylim([-0.01 0.01]);

% Frequency Response with distortion of recording
fig9 = figure(4);
for a = Start_index:Number_of_Recordings
    [H, freq_H] = MyFFT(h_t_array(:,a),fs);
    semilogx(freq_H,fftshift(20*log10(abs(H)./(Std_P)))); hold on;
    xlim([100,10000]);
    grid on
end
title('Frequency Response With no Windowing'); xlabel('Frequency (Hz)'); ylabel('SPL (dB)');

%% Extra Plotting Things
% % Inverse Filter
% fig6 = figure(123);
% plot(t,z_t); title('Inverse Filter?'); xlabel('Time'); ylabel('Magnitude');
% Nice Plotting?
% figs = [fig1, fig2]; %fig3, fig4, fig5, fig6, fig7, fig8, fig9];   %as many as needed
% nfig = length(figs);
% frac = 1/nfig;
% for K = 1 : nfig
%   old_pos = get(figs(K), 'Position');
%   set(figs(K), 'Position', [(K-1)*frac, old_pos(2), frac, old_pos(4)]);
% end