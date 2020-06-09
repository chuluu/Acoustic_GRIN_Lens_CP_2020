%% Intro
% Goal: 
% Given the parameter sweep of pillar_a dimension, we can figure out the
% dimensions based on the secant plot we want. Compare desired index of
% refraction desired with calculated index of refraction. 

% Inputs:
% d - Thickness of unit cell
% T & R - Data for varying pillar design (a dimension) and for a range of freq

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
title('Fetching data!!!','Fontsize', 16);
xlabel('Bark! Bark!','Fontsize', 16);
pause(1)

%% Input 

d = 0.013;

% Pillar_a parameters broken up
Pillar_a1 = 1:0.1:3;
Pillar_a2 = 3.1:0.05:5.4;
Pillar_a3 = 5.5:0.01:6.2;
Pillar_a = [Pillar_a1,Pillar_a2,Pillar_a3];

% Data to input
Data = csvread('Design_MRX_UC_13_PA_1_6_2mm_Varying_Step_Size.csv',5,0);

%% Ideal index of refraction secant graphing
h = 0.1625;
n_o = 1.965;
n_h = 1.001;
y = linspace(-h,h,25);

alpha = (1/h).*acosh(n_o/n_h);
n_y = n_o.*sech(alpha*y);
stem(y,n_y)

%% Data from Simulation COMSOL import for Sample B
parameter = Data(:,1);
freq_orig   = Data(:,2);
Mag_orig    = Data(:,3:4);   
Angles_orig = Data(:,5:6); %[-8.8247, 81.187]; % T , R
b = 1;
c = 1;
e = 1;

% Parse data for script to work
for a = 2:length(parameter)
    if parameter(a) ~= parameter(a-1)
        b = b + 1;
        c = 1;
        e = e + 2;
    end
    freq(c,b) = freq_orig(a);
    Mag(c,e:e+1) = Mag_orig(a,1:2);
    Angles(c,e:e+1) = Angles_orig(a,1:2);
    c = c+1;
end

% Add first index to make work
Mag(2:26,1:2) = Mag(1:25,1:2);
Mag(1,1:2) = Mag_orig(1,1:2);
Angles(2:26,1:2) = Angles(1:25,1:2);
Angles(1,1:2) = Angles_orig(1,1:2);
Angles_rad = deg2rad(Angles);
freq(2:26,1) = freq(1:25,1);
freq(1,1) = freq_orig(1);

%% Calculate Effective Parameters for all parametric sweep things
[row,col] = size(freq);
b = 1;
for a = 1:col
    T(:,a) = Mag(:,b).*cos(Angles_rad(:,b)) + 1i*Mag(:,b).*sin(Angles_rad(:,b));
    R(:,a) = Mag(:,b+1).*cos(Angles_rad(:,b+1)) + 1i*Mag(:,b+1).*sin(Angles_rad(:,b+1));
    b = b + 2;
    [Z(:,a),n(:,a),B_eff(:,a),p_eff(:,a)] = effective_material_derivation(T(:,a),R(:,a),d,freq(:,1));
end

%% Sweep through secant pattern at 2.3kHz - compare index of refraction with ideal to find dimensions
n_2_3k = real(n(19,:));
half_n_y = n_y(1:(length(n_y)/2)+1);
Z_2_3k = real(Z(19,:));

b = 1;
del = 0.008;
for a = 1:length(n_2_3k)
    if half_n_y(b) - del < n_2_3k(a) && n_2_3k(a) < half_n_y(b) + del
        Pillar_a_ideal(b) = Pillar_a(a);
        n_2_3k_ideal(b) = n_2_3k(a);
        Z_big(b) = Z_2_3k(a);
        b = b+1;
    end
    if b == 7
        del = 0.01;
    end
end

Pillar_a_scale = Pillar_a_ideal./max(Pillar_a_ideal);
Pillar_a_scale = Pillar_a_scale';
Pillar_a_ideal_plot = Pillar_a_ideal';

%% Plotting
figure(1);
C = linspecer(6);
axes('NextPlot','replacechildren', 'ColorOrder',C);
subplot(2,1,1); 
plot(Pillar_a,real(n(19,:)),'Linewidth',2.5); hold on
title('Effective index of refraction','Fontsize',14);
xlabel('dimension a for pillar (mm)','Fontsize',14);
ylabel('Magnitude (dimensionless)','Fontsize',14);
xlim([min(Pillar_a),max(Pillar_a)]);
grid on;

C = linspecer(6);
axes('NextPlot','replacechildren', 'ColorOrder',C);
subplot(2,1,2); 
plot(Pillar_a,real(Z(19,:)),'Linewidth',2.5); hold on
title('Effective Relative Acoustic Impedance','Fontsize',14);
xlabel('dimension a for pillar (mm)','Fontsize',14);
ylabel('Magnitude (dimensionless)','Fontsize',14);
xlim([min(Pillar_a),max(Pillar_a)]);
grid on;

figure(3)
stem(y,n_y); hold on
plot(y,n_y,'Linewidth',1.6)
title('Secant Discretized Plot','Fontsize',14);
xlabel('y - axis (m)','Fontsize',14);
ylabel('Index of Refraction (dimensionless)','Fontsize',14);
ylim([n_h,n_o]);
grid on;

figure(4);
plot(Pillar_a,real(n(19,:)),'Linewidth',2.5); hold on;
stem(Pillar_a_ideal,half_n_y,'Linewidth',2.5);
title('Index of Refraction vs. Pillar a dimension','Fontsize',14);
xlabel('dimension a for pillar (mm)','Fontsize',14);
ylabel('Index of Refraction (dimensionless)','Fontsize',14);
grid on;