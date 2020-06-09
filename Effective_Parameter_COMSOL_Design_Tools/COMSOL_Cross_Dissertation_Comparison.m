%% Intro
% Compare Paper with my model

clear variables
close all
spy 
title('He is fetching data','Fontsize', 16);
xlabel('Bark! Bark!','Fontsize', 16);

%% Import Data
d = 0.012; % 11.4cm

Data = csvread('Dissertation_Unit_Cell_Cross_n=2.csv',5,0);
freq   = Data(:,1);
Mag    = Data(:,2:3); 
Angles = Data(:,4:5); %[-8.8247, 81.187]; % T , R
Angles_rad = deg2rad(Angles);

T = Mag(:,1).*cos(Angles_rad(:,1)) + 1i*Mag(:,1).*sin(Angles_rad(:,1));
R = Mag(:,2).*cos(Angles_rad(:,2)) + 1i*Mag(:,2).*sin(Angles_rad(:,2));
[Z,n,B_eff,p_eff] =effective_material_derivation(T,R,d,freq);

%% Plotting
figure(2);
C = linspecer(6);
axes('NextPlot','replacechildren', 'ColorOrder',C);
plot(freq,real(n),'Linewidth',2.5); hold on
plot(freq,real(Z),'Linewidth',2.5); 
title('Effective index of refraction','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Magnitude (unitless)','Fontsize',14);
legend('Index of Refraction (n)','Acoustic Impedance (Z)','Location','NortheastOutside');
grid on;