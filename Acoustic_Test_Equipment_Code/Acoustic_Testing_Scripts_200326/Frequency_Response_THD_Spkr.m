clc
clear

%% Constants
Std_P = 2*10^-5; % Standard Reference Pressure in Pa
fs = 51200;      % Sampling Rate of the DAQ: NI 9250

%% Inputs
load Sensitivity.mat;
disp(['This Program will excite a speaker to extract the Total Harmonic Distortion']);
disp(['and Frequency Response of said speaker. plots are shown of results']);
disp([' ']);
disp(['Main Parameters are fstart, fend, Time of Sweep, and Amplitude']);
disp(['Deafualt: fstart = 100, fend = 10000, Time of Sweep = 2, Amplitude = 1']);

prompt = ['Do you want default Parameters? yes (1), no (2): '];
choice = input(prompt);

if choice == 2
    f1 = input('fstart (Hz): ');              % Freq Start
    f2 = input('fstart (Hz): ');              % Freq End
    T =  input('Time of Sweep (s): ');        % Total Time Sweep Takes
    A =  input('Amplitude: ');                % Amplitude for Sine Wave
    Zeroes_Delay = 17000; % Delays for code delay to work
    LW = 1.6;             % Line Width For Graphs
    Window_Sample_Length = 10000;
else
    Zeroes_Delay = 17000; % Delays for code delay to work
    f1 = 100;             % Freq Start
    f2 = 10000;           % Freq End
    LW = 1.6;             % Line Width For Graphs
    T = 2;                % Total Time Sweep Takes
    A = 1;                % Amplitude for Sine Wave
    LW = 1.6;             % Line Width For Graphs
    Window_Sample_Length = 10000;
end

disp(['Note there are special paramters, but change within program itself']);
dist_from_mic = input(' Distance of Speaker from mic (cm): ');

t_recording_name = ['t_recorded_A_',num2str(A),'_dist_',num2str(dist_from_mic),'_',date,'.csv'];
y_recording_name = ['y_recorded_A_',num2str(A),'_dist_',num2str(dist_from_mic),'_',date,'.csv'];

%% Farina Chirp Generation
w1 = 2*pi*f1;
w2 = 2*pi*f2;
t_delay = 0:1/fs:(Zeroes_Delay/fs)+T;   % gets timestep
t = 0:1/fs:T;                           % gets timestep
x = A.*sin((w1*T)/log(w2/w1).*(exp((t./T)*log(w2/w1))-1));
x_prime = [zeros(1,Zeroes_Delay),x];    % delay with zeros

figure(1); 
set(0,'DefaultFigureWindowStyle','docked')
spectrogram(x,hamming(1600),1600/2,2048,fs,'yaxis');  % Analyze input
set(gca,'YScale','log')
title('Farina Chirp Input');
%% Obtain and set up Devices
devices = daq.getDevices;

devices(3);

s = daq.createSession('ni');
addAnalogInputChannel(s,'cDAQ1Mod1', 0, 'Microphone');

s.Channels.Sensitivity = Sensitivity;
s.Channels

s.DurationInSeconds = T+0.5;                   % Length of Recording
s.Rate = fs;
sound(x_prime,fs);
[y_recorded,t_recorded] = s.startForeground(); % Record

Recording = audioplayer(y_recorded, s.Rate);
y_recorded = y_recorded';
t_recorded = t_recorded';

csvwrite(t_recording_name,t_recorded) % Write for polar plotting and saving data
csvwrite(y_recording_name,y_recorded) % Write for polar plotting and saving data

[loop_size,asdasd] = size(y_recorded);
if loop_size > 40
    loop_size = asdasd;
end

%% Inverse Filtering exp manipulation
x_inv = flip(x);
k = exp((t*log(f2/f1))./T);
z_t = x_inv./k;
% figure(41234);
% plot(t,z_t); title('Inverse Filter?'); xlabel('Time'); ylabel('Magnitude');

%% Time Domain Chirp Recording
% C = linspecer(4);
% axes('NextPlot','replacechildren', 'ColorOrder',C);
title_recording = ['Time Domain Farina Chirp Recording @ ', num2str(dist_from_mic), ' cm and A = ', num2str(A)];
figure(2);
subplot(2,1,1);
y_recorded_2 = zeros(1,length(y_recorded(1000:end-21000)));
t_recorded_2 = zeros(1,length(y_recorded(1000:end-21000)));

y_recorded_interm = y_recorded(1000:end-21000);
t_recorded_interm = t_recorded(1000:end-21000);
y_recorded_2 = y_recorded_interm;
t_recorded_2 = t_recorded_interm;
plot(t_recorded_interm,y_recorded_interm,'Linewidth',LW); hold on; 
title(title_recording,'Fontsize',14);
xlabel('Time (s)','Fontsize',14); ylabel('Magnitude','Fontsize',14);

%% Spectrogram Recording
title_FC_Spec_recording = ['Farina Chirp Recording Spectrogram @ ', num2str(dist_from_mic), ' cm and A = ', num2str(A)];
subplot(2,1,2);
spectrogram(y_recorded(1,:),hamming(1600),1600/2,2048,fs,'yaxis'); 
set(gca,'YScale','log')
title(title_FC_Spec_recording,'Fontsize',14);

%% Impulse Response Plotting
title_IR = ['h(t) impulse response @ ', num2str(dist_from_mic), ' cm and A = ', num2str(A)];
figure(3);
subplot(2,1,1);
[h_t,time_h_t] = MyFFTConv(y_recorded_2,z_t,fs);
h_t_array = h_t;
plot(time_h_t,h_t,'Linewidth',LW); hold on;
title(title_IR,'Fontsize',14); xlabel('Time (s)','Fontsize',14); ylabel('Magnitude','Fontsize',14);

xlim([0 5.6]);

%% Spectrogram of transfer function
title_IR_Spec = ['h(t) impulse response Spectrogram @ ', num2str(dist_from_mic), ' cm and A = ', num2str(A)];
subplot(2,1,2);
[s,f_h,t_h] = spectrogram(h_t,hamming(2000),2000/2,2048,fs,'yaxis');
spectrogram(h_t_array,hamming(2000),2000/2,2048,fs,'yaxis');
set(gca,'YScale','log')
title(title_IR_Spec,'Fontsize',14);
time_step = 2000/fs;
fundamental_time = 2.012/time_step;

%% Frequency Response Total (With Distortion)
FR_w_Distortion_title = ['SPL w/ Distortion FR @ ', num2str(dist_from_mic), ' cm and A = ', num2str(A)];
figure(4);
subplot(2,1,1);
[H_s, freq_H] = MyFFT(h_t_array,fs);
semilogx(freq_H,20*log10(fftshift(abs(H_s))./Std_P),'Linewidth',LW); hold on % SPL Plot 
title(FR_w_Distortion_title,'Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14); ylabel('Magnitude (dB)','Fontsize',14);
xlim([f1 f2]);
grid on;

%% Window the Fundamental Frequency Response and plot
FR_Fundy = ['SPL Fundamental FR @ ', num2str(dist_from_mic), ' cm and A = ', num2str(A)];
subplot(2,1,2);

% Find the Impulse Response Window
[~, Idx_h_t] = max(h_t);
t_val_fundy = time_h_t(Idx_h_t) - 0.01;
Time_window_start = find(time_h_t < t_val_fundy + 0.0001 & t_val_fundy - 0.0001 < time_h_t);
Time_window_start = Time_window_start(1);
% Time_window_start = find(time_h_t == 2.041113281250000); % Hard Coded
Time_window_end = Time_window_start + Window_Sample_Length;%find(time_h_t == 2.250761718750000);

h_t_Windowed = h_t_array(Time_window_start:Time_window_end);
time_h_t_Windowed = time_h_t(Time_window_start:Time_window_end);

% FFT the fundamental
[H_s_windowed_Fundy, freq_H_windowed_Fundy] = MyFFT(h_t_Windowed,fs);
semilogx(freq_H_windowed_Fundy,20*log10(fftshift(abs(H_s_windowed_Fundy))./Std_P),'Linewidth',LW); hold on % SPL Plot 
grid on;
xlim([f1 f2]);
title(FR_Fundy,'Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14); ylabel('Magnitude (dB)','Fontsize',14);

H_s_windowed_Fundy_array = H_s_windowed_Fundy;

%% Window the 3rd Harmonic
FR_3rd = ['SPL 3rd Harmonic FR @ ', num2str(dist_from_mic), ' cm and A = ', num2str(A)];
figure(5);
subplot(2,1,1);

% Find the Impulse Response Window
%Time_window_start = find(time_h_t == 1.552050781250000);
Distortion_Components = h_t(1:Time_window_start);
Distortion_Components_time = time_h_t(1:Time_window_start);
plot(Distortion_Components_time,Distortion_Components)

[~, Idx_h_t] = max(Distortion_Components);
t_val_3rd = time_h_t(Idx_h_t) - 0.01;

Time_window_start_3rd = find(Distortion_Components_time < t_val_3rd + 0.00001 & t_val_3rd - 0.00001 < Distortion_Components_time);
Time_window_start_3rd = Time_window_start_3rd(1);
% Time_window_start_3rd = find(time_h_t == 1.552050781250000); % Hard Coded
Time_window_end_3rd = Time_window_start_3rd + Window_Sample_Length;%find(time_h_t == 1.646933593750000);

% FFT the 3rd Harmonic
h_t_Windowed_3rd = Distortion_Components(Time_window_start_3rd:Time_window_end_3rd);
time_h_t_Windowed_3rd = Distortion_Components_time(Time_window_start_3rd:Time_window_end_3rd);


[H_s_windowed_3rd, freq_H_windowed_3rd] = MyFFT(h_t_Windowed_3rd,fs);
semilogx(freq_H_windowed_3rd,20*log10(fftshift(abs(H_s_windowed_3rd))./Std_P),'Linewidth',LW); hold on % SPL Plot 
grid on;
xlim([f1 f2]);

H_s_windowed_3rd_array = H_s_windowed_3rd;

title(FR_3rd,'Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14); ylabel('Magnitude (dB)','Fontsize',14);


%% THD Calculations!
title_THD = ['THD @ Only 3rd Harmonic ', num2str(dist_from_mic), ' cm and A = ', num2str(A)];
subplot(2,1,2);
THD = abs(H_s_windowed_3rd_array./H_s_windowed_Fundy_array);
semilogx(freq_H_windowed_3rd,fftshift(THD*100),'Linewidth',LW); hold on % SPL Plot 
ylim([0 30]);
xlim([f1 f2]);
grid on;

plot(freq_H_windowed_3rd,5*ones(1,length(freq_H_windowed_3rd)),'r','Linewidth',LW);
xlabel('Frequency (Hz)','Fontsize',14); ylabel('Percentage (%)','Fontsize',14);
title(title_THD,'Fontsize',14);
legend('data','5 % cutoff');

%% Compare Fundamental Windowing vs. non-fundy windowing and also 3rd harmonic and fundy
title_Fundy_Distortion = ['FR: Fundamental vs. w/ Distortion ', num2str(dist_from_mic), ' cm and A = ', num2str(A)];
figure(6);
subplot(2,1,1);
semilogx(freq_H,20*log10(fftshift(abs(H_s))./Std_P),'Linewidth',LW); hold on % SPL Plot 
semilogx(freq_H_windowed_Fundy,20*log10(fftshift(abs(H_s_windowed_Fundy))./Std_P),'Linewidth',LW); hold on % SPL Plot 
xlim([f1 f2]);
grid on;
title(title_Fundy_Distortion,'Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14); ylabel('Magnitude (dB)','Fontsize',14);
legend('FR with Distortion Content','FR just the Fundamental');

title_3rd_Fundy = ['FR: Fundamental vs. 3rd Harmonic ', num2str(dist_from_mic), ' cm and A = ', num2str(A)];
subplot(2,1,2);
semilogx(freq_H_windowed_Fundy,20*log10(fftshift(abs(H_s_windowed_Fundy_array(1,:)))./Std_P),'Linewidth',LW); hold on % SPL Plot 
semilogx(freq_H_windowed_3rd,20*log10(fftshift(abs(H_s_windowed_3rd_array(1,:)))./Std_P),'Linewidth',LW); hold on % SPL Plot 
xlim([f1 f2]);
grid on;
title(title_3rd_Fundy,'Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14); ylabel('Magnitude (dB)','Fontsize',14);
legend('Fundamental FR','3rd Harmonic FR');

xlim([f1 f2]);
