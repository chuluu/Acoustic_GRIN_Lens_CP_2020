%% Intro

clear variables
close all
spy 
title('He is fetching data','Fontsize', 16);
xlabel('Bark! Bark!','Fontsize', 16);
%% Sample A
Data = csvread('Paper_T_R_Data/Paper_T_R_Values.csv',1,0);
Data = Data([1:end-1],:);

Angles_SA = Data(:,[4:5]);
Angles_SA_rad = deg2rad(Angles_SA);
Mag_SA    = Data(:,[2:3]); 
freq   = Data(:,1);
d = 0.00953; % 0.794cm

T_SA = Mag_SA(:,1).*cos(Angles_SA_rad(:,1)) + 1i*Mag_SA(:,1).*sin(Angles_SA_rad(:,1));
R_SA = Mag_SA(:,2).*cos(Angles_SA_rad(:,2)) + 1i*Mag_SA(:,2).*sin(Angles_SA_rad(:,2));
[Z_SA,n_SA,B_eff_SA,p_eff_SA] = effective_material_derivation(T_SA,R_SA,d,freq);

%% Import T & R Data for Sample B to derive effective density and bulk modulus
Pap_Den_SB = csvread('Paper_density_y_dir.csv', 1, 0);
Pap_Bulk_SB = csvread('Paper_Bulk_Moudlus_y_dir.csv', 1, 0);

B_eff_pap_SB = Pap_Bulk_SB(:,2) + Pap_Bulk_SB(:,3).*1i;
Freq_B_eff_pap_SB = Pap_Bulk_SB(:,1);
B_eff_pap_SB = interp1(Freq_B_eff_pap_SB, B_eff_pap_SB, freq);


p_eff_pap_SB = Pap_Den_SB(:,2) + Pap_Den_SB(:,3).*1i;
Freq_p_eff_pap_SB = Pap_Den_SB(:,1);
p_eff_pap_SB = interp1(Freq_p_eff_pap_SB, p_eff_pap_SB, freq);

data_R_p_SB = csvread('Paper_T_R_Data/Reflection_Phase_y_dir_paper.csv',0,0);
data_T_p_SB = csvread('Paper_T_R_Data/Transmission_Phase_y_dir_paper.csv',0,0);
data_R_m_SB = csvread('Paper_T_R_Data/Reflection_Magnitude_y_dir_paper.csv',0,0);
data_T_m_SB = csvread('Paper_T_R_Data/Transmission_Magnitude_y_dir_paper.csv',0,0);

phase_interp_R_SB = webplot_digitizer_interpolater(data_R_p_SB, freq);
phase_interp_T_SB = webplot_digitizer_interpolater(data_T_p_SB, freq);

phase_interp_R_SB_rad = deg2rad(phase_interp_R_SB);
phase_interp_T_SB_rad = deg2rad(phase_interp_T_SB);

mag_interp_R_SB = webplot_digitizer_interpolater(data_R_m_SB, freq);
mag_interp_T_SB = webplot_digitizer_interpolater(data_T_m_SB, freq);

d = 0.00953; % 9.53cm

T_pap_SB = mag_interp_T_SB.*cos(phase_interp_T_SB_rad) + 1i*mag_interp_T_SB.*sin(phase_interp_T_SB_rad);
R_pap_SB = mag_interp_R_SB.*cos(phase_interp_R_SB_rad) + 1i*mag_interp_R_SB.*sin(phase_interp_R_SB_rad);
[Z_pap_SB,n_pap_SB,B_eff_pap_deriv_SB,p_eff_pap_deriv_SB] = effective_material_derivation(T_pap_SB,R_pap_SB,d,freq);

%% percent error
P_E_B = percent_error(B_eff_pap_deriv_SB,B_eff_pap_SB);
P_E_B_avg = sum(P_E_B)/length(P_E_B)

P_E_p = percent_error(p_eff_pap_deriv_SB,p_eff_pap_SB);
P_E_p_avg = sum(P_E_p)/length(P_E_p)
%% effective Density and Bulk Modulus
figure(1);
C = linspecer(4);
axes('NextPlot','replacechildren', 'ColorOrder',C);
subplot(2,1,1); 
plot(freq,real(B_eff_SA),'Linewidth',2.5); hold on
plot(freq,imag(B_eff_SA),'Linewidth',2.5);
plot(freq,1.18*ones(1,length(freq)),'Linewidth',2.5); 
plot(freq,0*ones(1,length(freq)),'Linewidth',2.5); 
title('Effective Bulk Modulus Using Paper T & R Sample A','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Magnitude','Fontsize',14);
legend('B eff real the script','B eff imaginary the script','B eff real paper data','B eff imaginary paper data','Location','NortheastOutside');
grid on;

C = linspecer(4);
axes('NextPlot','replacechildren', 'ColorOrder',C);
subplot(2,1,2); 
plot(freq,real(p_eff_SA),'Linewidth',2.5); hold on
plot(freq,imag(p_eff_SA),'Linewidth',2.5); 
plot(freq,1.2*ones(1,length(freq)),'Linewidth',2.5); 
plot(freq,0*ones(1,length(freq)),'Linewidth',2.5); 
title('Effective Density Using Paper T & R  Sample A','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Magnitude (unitless ratio)','Fontsize',14);
legend('p eff real the script','p eff imaginary the script','p eff real paper data','p eff imaginary paper data','Location','NortheastOutside');
grid on;

figure(2);
C = linspecer(6);
axes('NextPlot','replacechildren', 'ColorOrder',C);
subplot(2,1,1); 
plot(freq,real(B_eff_pap_deriv_SB),'Linewidth',2.5); hold on
plot(freq,imag(B_eff_pap_deriv_SB),'Linewidth',2.5);
plot(freq,real(B_eff_pap_SB),'Linewidth',2.5); hold on;
plot(freq,imag(B_eff_pap_SB),'Linewidth',2.5);
title('Effective Bulk Modulus Using Paper T & R  Sample B','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Magnitude (unitless ratio)','Fontsize',14);
legend('B eff real the script','B eff imaginary the script','B eff real paper data','B eff imaginary paper data','Location','NortheastOutside');
grid on;

C = linspecer(6);
axes('NextPlot','replacechildren', 'ColorOrder',C);
subplot(2,1,2); 
plot(freq,real(p_eff_pap_deriv_SB),'Linewidth',2.5); hold on
plot(freq,imag(p_eff_pap_deriv_SB),'Linewidth',2.5);
plot(freq,real(p_eff_pap_SB),'Linewidth',2.5); hold on;
plot(freq,imag(p_eff_pap_SB),'Linewidth',2.5);
title('Effective Density Using Paper T & R  Sample B','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Magnitude','Fontsize',14);
legend('p eff real the script','p eff imaginary the script','p eff real paper data','p eff imaginary paper data','Location','NortheastOutside');
grid on;

%% T & R Plotting
figure(3);
C = linspecer(2);
axes('NextPlot','replacechildren', 'ColorOrder',C);
subplot(2,1,1); 
plot(freq,Mag_SA(:,1),'Linewidth',2.5); hold on
plot(freq,Mag_SA(:,2),'Linewidth',2.5);
title('Magnitude Plot for Sample A','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Magnitude','Fontsize',14);
legend('Transmission','Reflection','Location','NortheastOutside');
grid on;

axes('NextPlot','replacechildren', 'ColorOrder',C);
subplot(2,1,2); 
plot(freq,Angles_SA(:,1),'Linewidth',2.5); hold on;
plot(freq,Angles_SA(:,2),'Linewidth',2.5); 
title('Phase Plot for Sample A','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Phase (deg)','Fontsize',14);
legend('Transmission','Reflection','Location','NortheastOutside');
grid on;

figure(4);
C = linspecer(2);
axes('NextPlot','replacechildren', 'ColorOrder',C);
subplot(2,1,1); 
plot(freq,mag_interp_T_SB,'Linewidth',2.5); hold on
plot(freq,mag_interp_R_SB,'Linewidth',2.5);
title('Magnitude Plot for Sample B','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Magnitude','Fontsize',14);
legend('Transmission','Reflection','Location','NortheastOutside');
grid on;

C = linspecer(2);
axes('NextPlot','replacechildren', 'ColorOrder',C);
subplot(2,1,2); 
plot(freq,phase_interp_T_SB,'Linewidth',2.5); hold on;
plot(freq,phase_interp_R_SB,'Linewidth',2.5); 
title('Phase Plot for Sample B','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Phase (deg)','Fontsize',14);
legend('Transmission','Reflection','Location','NortheastOutside');
grid on;