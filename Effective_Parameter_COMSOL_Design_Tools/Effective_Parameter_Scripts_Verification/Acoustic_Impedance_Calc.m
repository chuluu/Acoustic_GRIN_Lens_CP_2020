%% Calculate Acoustic Velocity of a Material
K = 1.825*10^9; %(in Pa)
p = 1180;  %(in kg/m^3)
V = sqrt(K/p);  %Solves for wave velocity
disp(['Acoustic Velocity: ', num2str(V),'m/s'])

%% Calculate Acoustic impedance
% Z = sqrt(p*B)
% Z = Z = p*V
V = 2730;
Z = p*V;  %(in (Pa*s)/m) or (in kg/(s*m^2)) same thing
disp(['Acoustic Impedance: ', num2str(Z),'(Pa*s)/m'])


%% Three boundary/ Transmission and Reflection coefficients
Z1 = 417.0731; %Air Impedance
Z2 = 3221400; %Acryllic Impedance
Z3 = 417.0731; 
wavelength = .1;
k = (2*pi)/wavelength;
l = .0055;

X1 = 1 + (Z1/Z3);
X2 = ((Z2/Z3) + (Z1/Z2));
T = (2)/((X1*cos(k*l)) + (j*X2*sin(k*l)));
T_mag = sqrt(real(T)^2 + imag(T)^2)

X1_N = 1 - (Z1/Z3);
X2_N = ((Z2/Z3) - (Z1/Z2));
R = ((X1_N*cos(k*l)) + (j*X2_N*sin(k*l)))/((X1*cos(k*l)) + (j*X2*sin(k*l)));
R_mag = sqrt(real(R)^2 + imag(R)^2)


%% Calculate effective index of refraction and acoustic impedance
m = 1;
n = (1/(k*l))*acos((1/(2*T))*(1-(R^2)+(T^2)) + (2*pi*m));
Z = sqrt(((1+R)^2 - (T^2))/((1-R)^2 - (T^2)));

n_mag = sqrt(real(n)^2 + imag(n)^2)
Z_mag = sqrt(real(Z)^2 + imag(Z)^2)

Z_norm = Z_mag/Z1


%% effective Density and Bluk Modulus
p_eff = n*Z;
p_eff_mag = sqrt(real(p_eff)^2 + imag(p_eff)^2)
B_eff = Z/n;
B_eff_mag = sqrt(real(B_eff)^2 + imag(B_eff)^2)


