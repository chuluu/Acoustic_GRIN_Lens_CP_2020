%% Sample B Plotting
clear variables
close all
spy 
title('He is fetching data','Fontsize', 16);
xlabel('Bark! Bark!','Fontsize', 16);

%% Data from Simulation COMSOL import for Sample B
d = 0.00953; % 9.53cm
Data = csvread('Paper_Replication_T_R_Test_Sample_A_200216.csv',5,0);



freq   = Data(:,1);
Mag    = Data(:,2:3); 
Angles = Data(:,4:5); %[-8.8247, 81.187]; % T , R
Angles_rad = deg2rad(Angles);

T = Mag(:,1).*cos(Angles_rad(:,1)) + 1i*Mag(:,1).*sin(Angles_rad(:,1));
R = Mag(:,2).*cos(Angles_rad(:,2)) + 1i*Mag(:,2).*sin(Angles_rad(:,2));
[Z,n,B_eff,p_eff] =effective_material_derivation(T,R,d,freq);

%% Sample A
Data_pap_SA = csvread('Paper_T_R_Data/Paper_T_R_Values.csv',1,0);
Data_pap_SA = Data_pap_SA([1:end-1],:);

Angles_SA = Data_pap_SA(:,[4:5]);
Angles_SA_rad = deg2rad(Angles_SA);
Mag_SA    = Data_pap_SA(:,[2:3]); 
frequency      = Data_pap_SA(:,1);

Angles_SA_interp = interp1(frequency,Angles_SA,freq);
Mag_SA_interp = interp1(frequency,Mag_SA,freq);


%% Percent Error
P_E_T_mag = percent_error(abs(T),Mag_SA_interp(:,1));
P_E_T_mag_avg = sum(P_E_T_mag)/length(P_E_T_mag)

P_E_R_mag = percent_error(abs(R),Mag_SA_interp(:,2));
P_E_R_mag_avg = sum(P_E_R_mag)/length(P_E_R_mag)

P_E_T_ang = percent_error(rad2deg(angle(T)),Angles_SA_interp(:,1));
P_E_T_ang_avg = sum(P_E_T_ang)/length(P_E_T_ang)

P_E_R_ang = percent_error(rad2deg(angle(R)),Angles_SA_interp(:,2));
P_E_R_ang_avg = sum(P_E_R_ang)/length(P_E_R_ang)

%% 

P_E_T_B_eff = percent_error(B_eff',1.18*ones(1,length(freq)));
P_E_T_B_eff_avg = sum(P_E_T_B_eff)/length(P_E_T_B_eff)

P_E_T_p_eff = percent_error(p_eff',1.2*ones(1,length(freq)));
P_E_T_p_eff_avg = sum(P_E_T_p_eff)/length(P_E_T_p_eff)

%% Plotting
figure(1);
C = linspecer(6);
axes('NextPlot','replacechildren', 'ColorOrder',C);
subplot(2,1,1); 
plot(freq,real(B_eff),'Linewidth',2.5); hold on
plot(freq,imag(B_eff),'Linewidth',2.5);
plot(freq,1.18*ones(1,length(freq)),'Linewidth',2.5); 
plot(freq,0*ones(1,length(freq)),'Linewidth',2.5); 
title('Magnitude vs. Freq: Effective Relative Bulk Modulus Sample A','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Magnitude (dimensionless)','Fontsize',14);
legend('K eff real Simulation','K eff imaginary Simulation',...
    'K eff real paper','K eff imaginary paper','Location','NortheastOutside');
grid on;

C = linspecer(6);
axes('NextPlot','replacechildren', 'ColorOrder',C);
subplot(2,1,2); 
plot(freq,real(p_eff),'Linewidth',2.5); hold on
plot(freq,imag(p_eff),'Linewidth',2.5); 
plot(freq,1.2*ones(1,length(freq)),'Linewidth',2.5); 
plot(freq,0*ones(1,length(freq)),'Linewidth',2.5); 
title('Magnitude vs. Freq: Effective Relative Mass Density Sample A','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Magnitude (dimensionless)','Fontsize',14);
legend('p eff real Simulation','p eff imaginary Simulation',...
    'p eff real paper','p eff imaginary paper','Location','NortheastOutside');
grid on;

figure(3);
C = linspecer(6);
axes('NextPlot','replacechildren', 'ColorOrder',C);
subplot(2,1,1); 
plot(freq,Data(:,2),'Linewidth',2.5); hold on
plot(freq,Data(:,3),'Linewidth',2.5); 
plot(freq,Mag_SA_interp(:,1),'Linewidth',2.5);
plot(freq,Mag_SA_interp(:,2),'Linewidth',2.5);
title('Magnitude Comparing Paper and Simulation T & R Sample A','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Magnitude (dimensionless)','Fontsize',14);
legend('Simulation T Magnitude','Simulation R Magnitude','Paper T Magnitude',...
    'Paper R Magnitude','Location','NortheastOutside');
grid on;

C = linspecer(6);
axes('NextPlot','replacechildren', 'ColorOrder',C);
subplot(2,1,2); 
plot(freq,Data(:,4),'Linewidth',2.5); hold on
plot(freq,Data(:,5),'Linewidth',2.5); 
plot(freq,Angles_SA_interp(:,1),'Linewidth',2.5);
plot(freq,Angles_SA_interp(:,2),'Linewidth',2.5);
title('Phase Comparing Paper and Simulation T & R Sample A','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Phase (deg)','Fontsize',14);
legend('Simulation T Phase','Simulation R Phase','Paper T Phase',...
    'Paper R Phase','Location','NortheastOutside');
grid on;