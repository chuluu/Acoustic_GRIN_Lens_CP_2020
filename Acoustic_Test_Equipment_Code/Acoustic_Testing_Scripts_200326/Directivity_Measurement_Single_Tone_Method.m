%% Setup
%{ 
This program plots a polar pattern utilizing a single frequency tone over
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
%% Extra things and comments for me
%{ 

Pass: 266379 Baylor Comp
Settling time for motor platform is like 7s therefore need to run the test at 7s and take
data for 3s, so 15s total

%}

%% Variables to Input
Std_P = 2*10^-5;
Number_of_Recordings = 37;
load Sensitivity.mat;
Start_index = 1;
fs = 51200;
f = 1000;
T = 1;
Name_dB_SPL = ['1000Hz_Max_dB_Data_Reflections', '_200519.csv']; % Out Data name
Name_Raw_Data = ['1000Hz_Raw_Data_Reflections', '_200519.csv']; % Out Data name

%% Create Stimulus
Zeroes_Delay = 17000;
t_delay = 0:1/fs:(Zeroes_Delay/fs)+T; % gets timestep
t = 0:1/fs:T; % gets timestep

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
A = [1*ones(1,Number_of_Recordings)];
for a = 1:Number_of_Recordings
    % While loop to wait for single scan
    while true
    data = Digital_Read.inputSingleScan();
    if data(1)==1;
        disp(['Index: ',num2str(a)]);   % Delay until Pulse Found
        break;
    end
    end
    
    % play sine wave tone and record the signal (x = test tone, y = recorded
    % tone)
    x = A(a).*sin(2*pi*f*t);             % Create tone
    x_prime = [zeros(1,Zeroes_Delay),x]; % delay with zeros
    sound(x_prime,fs);
    [y_recorded_val,t_recorded_val] = Microphone_Recording_Session.startForeground();
    for b = 1:length(y_recorded_val)
        y_recorded(b,a) = y_recorded_val(b);
    end
end

%% Post Processing - Windowing, Max SPL, and FFT
win = hamming(length(y_recorded(:,1)));
Windowed_signal = zeros(length(y_recorded(:,1)),Number_of_Recordings);
dB_SPL_Sine = zeros(1,Number_of_Recordings);
Y_Array = zeros(length(y_recorded(:,1)),Number_of_Recordings);

% Window the signal
for a = Start_index:Number_of_Recordings
    % Windowing
    Windowed_signal_val = y_recorded(:,a).*win;
    Windowed_signal(:,a) = Windowed_signal_val;
    
    % Max Value Section
    Max_Pressure = max(Windowed_signal_val);
    dB_SPL_Sine(a) = 20*log10(abs(Max_Pressure./(Std_P)));
    
    % FFT Section
    [Y_val, freq_Y_val] = MyFFT(Windowed_signal_val,fs);
    Y_Array(:,a) = Y_val;
end

%% Plot Recordings
figure(1);
for a = Start_index:Number_of_Recordings
    plot(t_recorded_val,y_recorded(:,a)); hold on;
end
xlabel('Time'); ylabel('Magnitude'); title('Recording');


% Plot Windowed Signal
figure(2)
for a = 2:Number_of_Recordings
    plot(t_recorded_val,Windowed_signal(:,a)); hold on;
end
xlabel('Time'); ylabel('Magnitude'); title('Recording Windowed');

% Spectrogram Recording
figure(3)
spectrogram(Windowed_signal(:,2),hamming(1600),1600/2,2048,fs,'yaxis'); 
set(gca,'YScale','log')
title('Farina Chirp Output');

% Plot FFT Response of Signal
figure(4)
for a = 2:Number_of_Recordings
    plot(freq_Y_val,fftshift(abs(Y_Array(:,a)))); hold on;
end
xlabel('Frequency (Hz)'); ylabel('Magnitude'); title('FFT Of Signal');
csvwrite(Name_dB_SPL,dB_SPL_Sine') % Write for polar plotting and saving data
csvwrite(Name_Raw_Data,y_recorded') % Write for polar plotting and saving data


%% Polar Plotting Section - Data Setup
Data_3 = dB_SPL_Sine(1:37);
Data_3_pressure = (10.^(Data_3/20)) * (Std_P);
Omindirectional = ones(1,37);
Pressure_Magnitude = Data_3_pressure;
Angles = ((pi*5/180).*(0:1:36/2))';
Angles_pos = Angles';
Angles_2 = ((-pi*5/180).*(0:1:36/2))';
Angles_neg = fliplr(Angles_2');
Angles_total = [Angles_neg(1:18),Angles_pos];


Angles = ((pi*5/180).*(0:1:36))';
Angles_deg = ((5).*(0:1:36))';
Angles_rad = deg2rad(Angles_deg);
[dummy,Meas_num] = size(Pressure_Magnitude);
delta_angle = 5;

%% Polar Plotting Section - Calculate Polar Plot array and DI
data_set = 1;
Pressure_Magnitude_Dataset = Pressure_Magnitude(:);
max_pressure = max(Pressure_Magnitude_Dataset);
if (max_pressure ~= Pressure_Magnitude(19))
    disp('Check');
end

% DI Utilized Beranek Acoustic book:
Constant = (4*pi*(max_pressure^2)*57.3)/(2*pi);
Last_value = 0;
for a = 1:length(Angles_deg)
    sum = Last_value + (Pressure_Magnitude_Dataset(a)^2)*sin(Angles_rad(a))*delta_angle;
    Last_value = sum;
end
sum = abs(sum);
Q_f = Constant/sum;
DI = 10*log10(Q_f)

%% Polar Plotting Section -  Set to dB
Normalized_dB = 20.*log(Pressure_Magnitude_Dataset./max_pressure);
max_val = max(Normalized_dB);
min_val = min(Normalized_dB);
wiggle = max_val - min_val;
radius_limits = [min_val-.1*wiggle, max_val+.1*wiggle];

%% Polar Plotting Stuff
figure(5);
polarplot(Angles,Normalized_dB);  hold on;
title('My Polar Plot');
rlim(radius_limits);
legend('2500 Hz');
title('Directivity Baseline'); 

%% Omindirectional case check
polarplot(Angles,Omindirectional);
Constant = (4*pi*(1^2)*57.3)/(2*pi);
Last_value = 0;
for a = 1:length(Angles_deg)
    sum = Last_value + (1^2)*sin(Angles_rad(a))*delta_angle;
    Last_value = sum;
end
sum = abs(sum);
Q_f = Constant/sum;
DI = 10*log10(Q_f)
