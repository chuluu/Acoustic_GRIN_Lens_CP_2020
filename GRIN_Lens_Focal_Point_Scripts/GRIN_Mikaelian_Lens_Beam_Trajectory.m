%% Intro
% The Mikaelian Lens is the GRIN lens with a Secant index of refraction
% pattern. This program plots the lens beam trajectory

clear variables
close all
figure(1);
spy
title('Fetching data!','Fontsize', 16);
xlabel('Bark! Bark!','Fontsize', 16);
pause(1);

%% Input
freq = 2300;     % Frequency
c= 343;          % Speed of Sound
lamda = c/freq;  % wavelength
k = (2*pi)/lamda;
h = 0.1625; % 0.1143
d = lamda/4;  % 0.0.0354
n_o = 2;%1.965;         % Index of refraction profile high
n_h = 1;%1.007;         % Index of refraction profile low
yo = 0:0.01:h;  % Initial condition

a = (1/h).*asech(n_h/n_o);
x = 0:0.001:d*8;

%% Beam Trajectory Plotting
close(figure(1));
figure(1);
C = linspecer(length(yo));
axes('NextPlot','replacechildren', 'ColorOrder',C);
for b = 1:length(yo)
    y = (1./a).*asinh(sinh(a.*yo(b)).*cos(a.*x));
    plot(x,y,'Linewidth',1.6); hold on;
    plot(x,-y,'Linewidth',1.6);
end
rectangle('Position',[0 -h d 2*h])
title('y-axis vs. x-axis (Beam Trajectory of Lens)','Fontsize',14);
xlabel('x-axis (m)','Fontsize',14);
ylabel('y-axis (m)','Fontsize',14);

x_f = x(y<0.001 & y>-0.001);
focal_length = x_f - d;
disp(['focal Length: ', num2str(focal_length*100), ' cm']);
grid on
