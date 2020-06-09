clc
clear all
%% Three boundary/ Transmission and Reflection coefficients
Z1 = 417.0731; %Air Impedance (Pa*s)/m
Z2 = (1.03*10^3)*2230; %ABS Impedance (Pa*s)/m
Z3 = 417.0731; %Air Impedance (Pa*s)/m

c1 = 343;
c2 = 2230; % m/s

freq = 2000; % Hz
k2 = (2*pi*freq)/c2; % in 1/m
l = 0.0159; % in meters

X1 = 1 + (Z1/Z3);
X2 = ((Z2/Z3) + (Z1/Z2));
T = (2)/((X1*cos(k2*l)) + (1i*X2*sin(k2*l)));
disp(['Transmission Coefficient: ', num2str(T)])


X1_N = 1 - (Z1/Z3);
X2_N = ((Z2/Z3) - (Z1/Z2));
R = ((X1_N*cos(k2*l)) + (1i*X2_N*sin(k2*l)))/((X1*cos(k2*l)) + (1i*X2*sin(k2*l)));
disp(['Reflection Coefficient: ', num2str(R)])

%% Eyeball parameters from Duke Paper
T_mag = 1; %mag
T_phase = -17.76; %deg

R_mag = 0.058; %mag
R_phase = 73.4; %deg

T_real = T_mag*cos(deg2rad(T_phase));
T_imag = T_mag*sin(deg2rad(T_phase));

R_real = R_mag*cos(deg2rad(R_phase));
R_imag = R_mag*sin(deg2rad(R_phase));

T = T_real + 1i*T_imag;
R = R_real + 1i*R_imag;

%% Calculate effective index of refraction and acoustic impedance
m = 0;
k0 = (2*pi*freq)/c1;

cos_expression = (1-(R^2)+(T^2))/(2*T);

n = acos(cos_expression)/(k0*l) + (2*pi*m)/(k0*l);
Z = sqrt(((1+R)^2 - (T^2))/((1-R)^2 - (T^2)));
disp(['Index of Refraction: ', num2str(n)])
disp(['Effective Acoustic Impedance: ', num2str(Z), ' (Pa*s)/m'])
disp(['Relative Real Effective Acoustic Impedance: ', num2str(real(Z/Z1))])



%% effective Density and Bluk Modulus
B_air = 0.101; %Pa
p_air = 1.225; %kg/m3

p_eff = n*Z;
B_eff = Z/n;
disp(['Relative Effective Density: ', num2str(p_eff)])

disp(['Relative Effective Density: ', num2str(abs(p_eff)/p_air)])
disp(['Relative Effective Bulk Modulus: ', num2str(abs(B_eff)/B_air)])