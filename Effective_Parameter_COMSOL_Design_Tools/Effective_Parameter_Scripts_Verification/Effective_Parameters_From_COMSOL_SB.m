%% Sample B Plotting
clear variables
close all
spy 
title('He is fetching data','Fontsize', 16);
xlabel('Bark! Bark!','Fontsize', 16);

%% Data from Simulation COMSOL import for Sample B
d = 0.00953; % 9.53cm

Data = csvread('Paper_Replication_T_R_Test_Sample_B_200127.csv',5,0);
freq   = Data(:,1);
Mag    = Data(:,2:3); 
Angles = Data(:,4:5); %[-8.8247, 81.187]; % T , R
Angles_rad = deg2rad(Angles);

T = Mag(:,1).*cos(Angles_rad(:,1)) + 1i*Mag(:,1).*sin(Angles_rad(:,1));
R = Mag(:,2).*cos(Angles_rad(:,2)) + 1i*Mag(:,2).*sin(Angles_rad(:,2));
[Z,n,B_eff,p_eff] =effective_material_derivation(T,R,d,freq);

%% Import T & R data for Sample B
data_R_p = csvread('Paper_T_R_Data/Reflection_Phase_y_dir_paper.csv',0,0);
data_T_p = csvread('Paper_T_R_Data/Transmission_Phase_y_dir_paper.csv',0,0);
data_R_m = csvread('Paper_T_R_Data/Reflection_Magnitude_y_dir_paper.csv',0,0);
data_T_m = csvread('Paper_T_R_Data/Transmission_Magnitude_y_dir_paper.csv',0,0);

phase_interp_R = webplot_digitizer_interpolater(data_R_p, freq);
phase_interp_T = webplot_digitizer_interpolater(data_T_p, freq);

mag_interp_R = webplot_digitizer_interpolater(data_R_m, freq);
mag_interp_T = webplot_digitizer_interpolater(data_T_m, freq);

%% Import Paper Bulk and density Data for Sample B
Pap_Den = csvread('Paper_density_y_dir.csv', 1, 0);
Pap_Bulk = csvread('Paper_Bulk_Moudlus_y_dir.csv', 1, 0);

B_eff_pap = Pap_Bulk(:,2) + Pap_Bulk(:,3).*1i;
Freq_B_eff_pap = Pap_Bulk(:,1);
B_eff_pap = interp1(Freq_B_eff_pap, B_eff_pap, freq);

p_eff_pap = Pap_Den(:,2) + Pap_Den(:,3).*1i;
Freq_p_eff_pap = Pap_Den(:,1);
p_eff_pap = interp1(Freq_p_eff_pap, p_eff_pap, freq);

%% Percent Error
P_E_T_mag = percent_error(mag_interp_T,Data(:,2));
P_E_T_mag_avg = sum(P_E_T_mag)/length(P_E_T_mag)

P_E_R_mag = percent_error(mag_interp_R,Data(:,3));
P_E_R_mag_avg = sum(P_E_R_mag)/length(P_E_R_mag)

P_E_T_ang = percent_error(phase_interp_T,Data(:,4));
P_E_T_ang_avg = sum(P_E_T_ang)/length(P_E_T_ang)

P_E_R_ang = percent_error(phase_interp_R,Data(:,5));
P_E_R_ang_avg = sum(P_E_R_ang)/length(P_E_R_ang)

P_E_T_B_eff = percent_error(B_eff,B_eff_pap);
P_E_T_B_eff_avg = sum(P_E_T_B_eff)/length(P_E_T_B_eff)

P_E_T_p_eff = percent_error(p_eff,p_eff_pap);
P_E_T_p_eff_avg = sum(P_E_T_p_eff)/length(P_E_T_p_eff)

%% Plotting
figure(1);
C = linspecer(6);
axes('NextPlot','replacechildren', 'ColorOrder',C);
subplot(2,1,1); 
plot(freq,real(B_eff),'Linewidth',2.5); hold on
plot(freq,imag(B_eff),'Linewidth',2.5);
plot(freq,real(B_eff_pap),'Linewidth',2.5); hold on;
plot(freq,imag(B_eff_pap),'Linewidth',2.5);
title('Magnitude vs. Freq: Effective Relative Bulk Modulus Sample B','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Magnitude','Fontsize',14);
legend('K eff real Simulation','K eff imaginary Simulation','K eff real paper','K eff imaginary paper','Location','NortheastOutside');
grid on;

C = linspecer(6);
axes('NextPlot','replacechildren', 'ColorOrder',C);
subplot(2,1,2); 
plot(freq,real(p_eff),'Linewidth',2.5); hold on
plot(freq,imag(p_eff),'Linewidth',2.5); 
plot(freq,real(p_eff_pap),'Linewidth',2.5);
plot(freq,imag(p_eff_pap),'Linewidth',2.5);
title('Magnitude vs. Freq: Effective Relative Mass Density Sample B','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Magnitude','Fontsize',14);
legend('p eff real Simulation','p eff imaginary Simulation','p eff real paper','p eff imaginary paper','Location','NortheastOutside');
grid on;

figure(3);
C = linspecer(6);
axes('NextPlot','replacechildren', 'ColorOrder',C);
subplot(2,1,1); 
plot(freq,Data(:,2),'Linewidth',2.5); hold on
plot(freq,Data(:,3),'Linewidth',2.5); 
plot(freq,mag_interp_T,'Linewidth',2.5);
plot(freq,mag_interp_R,'Linewidth',2.5);
title('Magnitude Comparing Paper and Simulation T & R Sample B','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Magnitude','Fontsize',14);
legend('Simulation T Magnitude','Simulation R Magnitude','Paper T Magnitude','Paper R Magnitude','Location','NortheastOutside');
grid on;

C = linspecer(6);
axes('NextPlot','replacechildren', 'ColorOrder',C);
subplot(2,1,2); 
plot(freq,Data(:,4),'Linewidth',2.5); hold on
plot(freq,Data(:,5),'Linewidth',2.5); 
plot(freq,phase_interp_T,'Linewidth',2.5);
plot(freq,phase_interp_R,'Linewidth',2.5);
title('Phase Comparing Paper and Simulation T & R Sample B','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Phase (deg)','Fontsize',14);
legend('Simulation T Phase','Simulation R Phase','Paper T Phase','Paper R Phase','Location','NortheastOutside');
grid on;