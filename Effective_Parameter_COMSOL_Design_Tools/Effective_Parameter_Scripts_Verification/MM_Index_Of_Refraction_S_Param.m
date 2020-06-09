%% Parameters
freq = 100;
rho1 = 1.29; %kg/m3
rho2 = 1.03*10^3; %kg/m3

B1 = 141; %KPA
B2 =2;

c1 = 343;  %m/s
c2 = 2230; %m/s

d = 0.005;

%% Reflection and transmission ceoff derivation
Z = (rho2*c2)/(rho1*c1);
k = (2*pi*freq)/c1;
m = rho2/rho1;
n = c1/c2;

R = (tan(n*k*d)*((1/Z) - Z)*1i)/(2 - tan(k*n*d)*((1/Z)+Z)*1i);
T = 2/(cos(n*k*d)*(2 - tan(k*n*d)*((1/Z)+Z)*1i));

%% index of refraction derivation
r = sqrt(((R^2)-(T^2)-1)^2 - (4*(T^2)));
x = (1-R^2 + T^2 +r)/(2*T);
Z = r/(1-(2*R) +(R^2) - (T^2))
n = (-1i*log10(x)+ 2*pi*m)/(k*d);
