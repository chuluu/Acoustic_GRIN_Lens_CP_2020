%% Inputs 
T_mag = 0.97; %mag
T_phase = 92; %deg

R_mag = 0.04; %mag
R_phase = -9; %deg

l = 0.050;    %in m length of device
freq = 1000;  %Frequency

%% constants
c1 = 343;      %m/s (speed of air)
Z1 = 420;      %(Pa*s)/m (impedance of air)
B_air = 0.101; %Pa
p_air = 1.225; %kg/m3

%% Eyeball parameters from Duke Paper
T_real = T_mag*cos(deg2rad(T_phase));
T_imag = T_mag*sin(deg2rad(T_phase));

R_real = R_mag*cos(deg2rad(R_phase));
R_imag = R_mag*sin(deg2rad(R_phase));

T = T_real + 1i*T_imag;
R = R_real + 1i*R_imag;

%% Calculate effective index of refraction and acoustic impedance
m = 0;
k0 = (2*pi*freq)/c1;

n = (1/(k0*l))*acos((1/(2*T))*(1-(R^2)-(T^2)) + (2*pi*m));
Z = sqrt(((1+R)^2 - (T^2))/((1-R)^2 - (T^2)));
disp(['Index of Refraction: ', num2str(abs(n))])
disp(['Effective Acoustic Impedance: ', num2str(Z), ' (Pa*s)/m'])
disp(['Relative Effective Acoustic Impedance: ', num2str(Z/Z1), ' (Pa*s)/m'])



%% effective Density and Bluk Modulus
p_eff = n*Z;
B_eff = Z/n;
disp(['Effective Relative Density: ', num2str(abs(p_eff)/p_air)])
disp(['Effective Relative Bulk Modulus: ', num2str(abs(B_eff)/B_air)])