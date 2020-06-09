%% Intro
% Plots frequency vs index of refraction and acoustic impedance in order to
% for a single parameter 

clear variables
close all
spy 
title('He is fetching data','Fontsize', 16);
xlabel('Bark! Bark!','Fontsize', 16);

%% Data from Simulation COMSOL import for Sample B
d = 0.012; % 11.4cm

Data = csvread('Dissertation_Unit_Cell_Cross_n=2.csv',5,0);
freq   = Data(:,1);
Mag    = Data(:,2:3); 
Angles = Data(:,4:5); %[-8.8247, 81.187]; % T , R
Angles_rad = deg2rad(Angles);

T = Mag(:,1).*cos(Angles_rad(:,1)) + 1i*Mag(:,1).*sin(Angles_rad(:,1));
R = Mag(:,2).*cos(Angles_rad(:,2)) + 1i*Mag(:,2).*sin(Angles_rad(:,2));
[Z,n,B_eff,p_eff] =effective_material_derivation(T,R,d,freq);

%% Import Data from Paper
Dissertation_data_n = csvread('Dissertation_n=2_Graph.csv',0,0);
Dissertation_data_n_interp = webplot_digitizer_interpolater(Dissertation_data_n, freq);

Dissertation_data_Z = csvread('Dissertation_n=2_Z_Graph.csv',0,0);
Dissertation_data_Z_interp = webplot_digitizer_interpolater(Dissertation_data_Z, freq);

%% Percent Difference
P_E_n = percent_error(n,Dissertation_data_n_interp);
P_E_n_avg = sum(P_E_n)/length(P_E_n)

P_E_Z = percent_error(Z,Dissertation_data_Z_interp);
P_E_Z_avg = sum(P_E_Z)/length(P_E_Z)

%% Plotting
figure(1);
C = linspecer(6);
axes('NextPlot','replacechildren', 'ColorOrder',C);
subplot(2,1,1); 
plot(freq,real(n),'Linewidth',2.5); hold on
plot(freq,Dissertation_data_n_interp,'Linewidth',2.5);
title('Effective index of refraction','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Magnitude (dimensionless)','Fontsize',14);
legend('n eff Simulation','n eff Paper','Location','NortheastOutside');
grid on;

C = linspecer(6);
axes('NextPlot','replacechildren', 'ColorOrder',C);
subplot(2,1,2); 
plot(freq,real(Z),'Linewidth',2.5); hold on
plot(freq,Dissertation_data_Z_interp,'Linewidth',2.5);
title('Effective Acoustic Impedance','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Magnitude (dimensionless)','Fontsize',14);
legend('Z eff Simulation','Z eff Paper','Location','NortheastOutside');
grid on;