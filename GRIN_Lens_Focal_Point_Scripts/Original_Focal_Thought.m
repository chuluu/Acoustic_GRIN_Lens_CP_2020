%% Program to do calculations of metamaterial grin lens design
%n_o = center at y = 0
%n_h = outermost index of refraction
%a = gradient parameter
%h = how tall the device is 
clc 
clear all

%% Presets that I chose based on this design
freq = 2000;
c= 343;
lamda = c/freq;

h = 0.1143;
d = 0.035;
n_o = 2;
n_h = 1;

y_o = [0.05,0.01];
x = 0:0.01:d;
y = -0.15:0.01:0.15;    %plotting points of long width

%% Index of refraction Secant Calculations
a = (1/h).*acosh(n_o/n_h);
n_y = index_of_refraction_calculation(y,a,n_o);

%% Beam Trajectory
% Hyp Transforms
u_o = sinh(a*y_o);
ud_o = a*cosh(a*y_o);

%% Solving for focal point Setup
[H_a_d, H_f_d, Hd_a_d, Hd_f_d] = H_Parameters(d,a);
[H_a, H_f, Hd_a, Hd_f] = H_Parameters(x,a);

%% Solve for Focal Point
u_d_val = [H_f_d, H_a_d; Hd_f_d, Hd_a_d]*[u_o; ud_o];
u_d = u_d_val(1);
ud_d = u_d_val(2);

%%
y_d = (1/a)*asinh(u_d);
yd_d = (1/(a*cosh(a*y_d)))*ud_d;
n_y_d = index_of_refraction_calculation(y_d,a,n_o);

Num_x_f = 1-((yd_d^2)*((n_y_d^2) - 1));
Den_x_f = n_y_d^2;

x_f = -(y_d/yd_d)*sqrt(Num_x_f/Den_x_f)
%% Plotting
figure(1)
plot_function(y,n_y,'Hyperbolic Secant Index of Refraction Profile','y-axis (m)','index of refraction',14);
grid on;

%%  quick focal point calculation MIT
% f = 1/(n_o*a*d); % From an MIT picture on gradient optics

%%
% Num = (u_o.*Hd_f);
% Den = a.*cosh(asinh(u_o.*H_f));
% yx_d = Num./Den;
% 
% y_x = (1/a).*asinh((u_o.*H_f) + (ud_o.*H_a));
% den = a.*cosh(asinh((u_o.*H_f) + (ud_o.*H_a)));
% yd_x = ((u_o.*Hd_f) + (ud_o.*Hd_a))./(den);

% w = 2*pi*1000;
% real(1i*w*exp(1i*w*1))
% (1i*w*exp(1i*w*1)-1i*w*exp(-1i*w*1))/2
