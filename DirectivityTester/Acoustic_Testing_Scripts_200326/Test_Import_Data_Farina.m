clc
clear

%% Farina Chirp Generation

LW = 1.6; % Linewidth
Std_P = 2*10^-5;
Zeroes_Delay = 17000;
fs = 51200;
f1 = 100;
f2 = 10000;

w1 = 2*pi*f1;
w2 = 2*pi*f2;
T = 2;
t_delay = 0:1/fs:(Zeroes_Delay/fs)+T;   % gets timestep
t = 0:1/fs:T;                           % gets timestep
A = 2; 
x = A.*sin((w1*T)/log(w2/w1).*(exp((t./T)*log(w2/w1))-1));
x_prime = [zeros(1,Zeroes_Delay),x];    % delay with zeros

figure(1); 
set(0,'DefaultFigureWindowStyle','docked')
spectrogram(x,hamming(1600),1600/2,2048,fs,'yaxis');  % Analyze input
set(gca,'YScale','log')
title('Farina Chirp Input','Fontsize',14);

%% Import Recorded Data
y_recorded_05 = importdata('y_recorded_A_05.csv');
t_recorded_05 = importdata('t_recorded_A_05.csv');
y_recorded_10 = importdata('y_recorded_A_1.csv');
t_recorded_10 = importdata('t_recorded_A_1.csv');
y_recorded_15 = importdata('y_recorded_A_15.csv');
t_recorded_15 = importdata('t_recorded_A_15.csv');
y_recorded_20 = importdata('y_recorded_A_2.csv');
t_recorded_20 = importdata('t_recorded_A_2.csv');

%% Get some recorded Stuff
y_recorded = [y_recorded_20;y_recorded_15;y_recorded_10;y_recorded_05];
t_recorded = [t_recorded_20;t_recorded_15;t_recorded_10;t_recorded_05];
[loop_size,asdasd] = size(y_recorded);

%% Check time responses
for a = 1:loop_size
    y_test(a,:) = [zeros(1,14200),y_recorded(a,:)];
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
figure(2);
set(0,'DefaultFigureWindowStyle','docked')
subplot(2,1,1);
for a = 1:loop_size
    y_recorded_interm = y_recorded(a,1000:end-21000);
    t_recorded_interm = t_recorded(a,1000:end-21000);
    y_recorded_2(a,:) = y_recorded_interm;
    t_recorded_2(a,:) = t_recorded_interm;
    plot(t_recorded_interm,y_recorded_interm,'Linewidth',LW); hold on; 
    title('Time Domain Farina Chirp Recording','Fontsize',14);
    xlabel('Time (s)','Fontsize',14); ylabel('Magnitude','Fontsize',14);
end
legend('A = 2','A = 1.5','A = 1.0','A = 0.5');

%% Spectrogram Recording
subplot(2,1,2);
spectrogram(y_recorded(1,:),hamming(1600),1600/2,2048,fs,'yaxis'); 
set(gca,'YScale','log')
title('Farina Chirp Recording Max Amplitude (A = 2)','Fontsize',14);


%% Impulse Response Plotting
figure(3);
subplot(2,1,1);
for a = 1:loop_size
    [h_t,time_h_t] = MyFFTConv(y_recorded_2(a,:),z_t,fs);
    h_t_array(a,:) = h_t;
    plot(time_h_t,h_t,'Linewidth',LW); hold on;
    title('h(t) impulse response','Fontsize',14); xlabel('Time (s)','Fontsize',14); ylabel('Magnitude','Fontsize',14);
end
legend('A = 2','A = 1.5','A = 1.0','A = 0.5');
xlim([0 5.6]);

%% Spectrogram of transfer function
subplot(2,1,2);
[s,f_h,t_h] = spectrogram(h_t,hamming(2000),2000/2,2048,fs,'yaxis');
spectrogram(h_t_array(1,:),hamming(2000),2000/2,2048,fs,'yaxis');
set(gca,'YScale','log')
title('Speaker Transfer Function Max Amplitude (A = 2)');
time_step = 2000/fs;
fundamental_time = 2.012/time_step;

%% Frequency Response Total (With Distortion)
figure(4);
subplot(2,1,1);
for a = 1:loop_size
    [H_s, freq_H] = MyFFT(h_t_array(a,:),fs);
    semilogx(freq_H,20*log10(fftshift(abs(H_s))./Std_P),'Linewidth',LW); hold on % SPL Plot 
    title('SPL w/ Distortion FR @ 0.4m','Fontsize',14);
    xlabel('Frequency (Hz)','Fontsize',14); ylabel('Magnitude (dB)','Fontsize',14);
    xlim([f1 f2]);
    grid on;
end
legend('A = 2','A = 1.5','A = 1.0','A = 0.5');

%legend('H(s) without windowing','H(s) with windowing','Location','best');


%% Window the Fundamental Frequency Response
subplot(2,1,2);
for a = 1:loop_size
    Time_window_start = find(time_h_t == 2.041113281250000);
    Time_window_end = Time_window_start + 10000;%find(time_h_t == 2.250761718750000);
    
    Time_window_start_3rd = find(time_h_t == 1.552050781250000);
    Time_window_end_3rd = Time_window_start + 10000;%find(time_h_t == 1.646933593750000);
    
    
    h_t_Windowed = h_t_array(a,Time_window_start:Time_window_end);
    time_h_t_Windowed = time_h_t(Time_window_start:Time_window_end);

    [H_s_windowed_Fundy, freq_H_windowed_Fundy] = MyFFT(h_t_Windowed,fs);
    semilogx(freq_H_windowed_Fundy,20*log10(fftshift(abs(H_s_windowed_Fundy))./Std_P),'Linewidth',LW); hold on % SPL Plot 
    grid on;
    xlim([f1 f2]);
    title('SPL Fundamental FR @ 0.4m','Fontsize',14);
    xlabel('Frequency (Hz)','Fontsize',14); ylabel('Magnitude (dB)','Fontsize',14);

    H_s_windowed_Fundy_array(a,:) = H_s_windowed_Fundy;
end
legend('A = 2','A = 1.5','A = 1.0','A = 0.5');

%% Window the 3rd Harmonic
figure(5);
subplot(2,1,1);
for a = 1:loop_size
    Time_window_start = find(time_h_t == 1.552050781250000);
    Time_window_end = Time_window_start + 10000;%find(time_h_t == 1.646933593750000);

    h_t_Windowed_3rd = h_t_array(a,Time_window_start:Time_window_end);
    time_h_t_Windowed_3rd = time_h_t(Time_window_start:Time_window_end);

    [H_s_windowed_3rd, freq_H_windowed_3rd] = MyFFT(h_t_Windowed_3rd,fs);
    semilogx(freq_H_windowed_3rd,20*log10(fftshift(abs(H_s_windowed_3rd))./Std_P),'Linewidth',LW); hold on % SPL Plot 
    grid on;
    xlim([f1 f2]);

    H_s_windowed_3rd_array(a,:) = H_s_windowed_3rd;
end
title('SPL 3rd Harmonic FR','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14); ylabel('Magnitude (dB)','Fontsize',14);
legend('A = 2','A = 1.5','A = 1.0','A = 0.5');


%% THD Calculations!

subplot(2,1,2);
for a = 1:loop_size
    THD = abs(H_s_windowed_3rd_array(a,:)./H_s_windowed_Fundy_array(a,:));
    semilogx(freq_H_windowed_3rd,fftshift(THD*100),'Linewidth',LW); hold on % SPL Plot 
    ylim([0 30]);
    xlim([f1 f2]);
    grid on;
end
plot(freq_H_windowed_3rd,5*ones(1,length(freq_H_windowed_3rd)),'r','Linewidth',LW);
xlabel('Frequency (Hz)','Fontsize',14); ylabel('Percentage (%)','Fontsize',14);
title('THD Only 3rd Harmonic','Fontsize',14);
legend('A = 2','A = 1.5','A = 1.0','A = 0.5','5% Cutoff Line');

%% Compare Fundamental Windowing vs. non-fundy windowing and also 3rd harmonic and fundy
figure(6);
subplot(2,1,1);
semilogx(freq_H,20*log10(fftshift(abs(H_s))./Std_P),'Linewidth',LW); hold on % SPL Plot 
semilogx(freq_H_windowed_Fundy,20*log10(fftshift(abs(H_s_windowed_Fundy))./Std_P),'Linewidth',LW); hold on % SPL Plot 
xlim([f1 f2]);
grid on;
title('Fundamental vs. Whole Data: Lowest Amplitude','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14); ylabel('Magnitude (dB)','Fontsize',14);
legend('FR with Distortion Content','FR just the Fundamental');

subplot(2,1,2);
semilogx(freq_H_windowed_Fundy,20*log10(fftshift(abs(H_s_windowed_Fundy_array(4,:)))./Std_P),'Linewidth',LW); hold on % SPL Plot 
semilogx(freq_H_windowed_3rd,20*log10(fftshift(abs(H_s_windowed_3rd_array(4,:)))./Std_P),'Linewidth',LW); hold on % SPL Plot 
xlim([f1 f2]);
grid on;
title('Fundamental vs. Whole Data: Lowest Amplitude','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14); ylabel('Magnitude (dB)','Fontsize',14);
legend('Fundamental FR','3rd Harmonic FR');

xlim([100 10000]);