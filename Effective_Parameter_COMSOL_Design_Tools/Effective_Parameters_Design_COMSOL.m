%% Intro
% Goal: 
% Plots frequency vs index of refraction and acoustic impedance
% for a single a parameter pillar for a given unit cell perimeter

% Inputs:
% d - Thickness of unit cell
% T & R - Data for a single pillar/unit cell design and for a range of freq

% Ouputs:
% Z - Relative effective acoustic impedance
% n - Index of Refraction
% B_eff - Relative effetive acoustic bulk modulus
% p_eff - Relative effetive acoustic mass density

% Constraints: 
% Effective unit cell dimensions constrain the pillar dimensions

clear variables
close all
spy 
title('He is fetching data','Fontsize', 16);
xlabel('Bark! Bark!','Fontsize', 16);

%% Data from Simulation COMSOL import for Sample B
d = 0.013; % 11.4cm

Data = csvread('Design_MRX_UC_13_n=1.csv',5,0);
freq   = Data(:,1);
Mag    = Data(:,2:3); 
Angles = Data(:,4:5); %[-8.8247, 81.187]; % T , R
Angles_rad = deg2rad(Angles);

T = Mag(:,1).*cos(Angles_rad(:,1)) + 1i*Mag(:,1).*sin(Angles_rad(:,1));
R = Mag(:,2).*cos(Angles_rad(:,2)) + 1i*Mag(:,2).*sin(Angles_rad(:,2));
[Z,n,B_eff,p_eff] =effective_material_derivation(T,R,d,freq);

%% error bound
n_lower = n(19) - n(19)*0.05;
n_upper = n(19) + n(19)*0.05;

Z_upper = Z(19) + Z(19)*0.05;
Z_lower = Z(19) - Z(19)*0.05;

%% Plotting
figure(1);
C = linspecer(6);
axes('NextPlot','replacechildren', 'ColorOrder',C);
subplot(2,1,1); 
plot(freq,real(n),'Linewidth',2.5); hold on
plot(freq,n_upper*ones(1,length(freq)),'Linewidth',2.5); hold on
plot(freq,n_lower*ones(1,length(freq)),'Linewidth',2.5); hold on

%plot(freq,imag(n),'Linewidth',2.5);
title('Effective index of refraction vs. Frequency','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Magnitude (dimensionless)','Fontsize',14);
legend('n eff Simulation','n upper 5% bound',' n lower 5% bound','Location','NortheastOutside');
grid on;
ylim([0.5 1.5]);

C = linspecer(6);
axes('NextPlot','replacechildren', 'ColorOrder',C);
subplot(2,1,2); 
plot(freq,real(Z),'Linewidth',2.5); hold on
plot(freq,Z_upper*ones(1,length(freq)),'Linewidth',2.5); hold on
plot(freq,Z_lower*ones(1,length(freq)),'Linewidth',2.5); hold on
%plot(freq,imag(Z),'Linewidth',2.5); 
title('Effective Acoustic Impedance vs. Frequency','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Magnitude (dimensionless)','Fontsize',14);
legend('Z eff Simulation','Z upper 5% bound',' Z lower 5% bound','Location','NortheastOutside');
ylim([0.5 1.5]);
grid on;