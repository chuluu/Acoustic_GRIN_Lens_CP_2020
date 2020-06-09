%% Intro

clear variables
close all
spy 
title('He is fetching data','Fontsize', 16);
xlabel('Bark! Bark!','Fontsize', 16);

%% Presets that I chose based on this design
h = 0.15525; %Long width/2 in m
d = 0.054;  %shorter distance in m
n_o = 12.4;    % max index of refraction
n_h = 1;    % min index of refraction
y = linspace(-h,h,23);    % plotting points of long width

%function n = ideal_refraction(n_o,n_h,y)
a = (1/h).*acosh(n_o/n_h);
n_y = n_o.*sech(a*y); % index of refraction along the GRIN Lens
figure(1)
stem(y,n_y);

%% Presets that I chose based on this design
h = 0.15525; %Long width/2 in m
d = 0.054;  %shorter distance in m
n_o = 1.98;    % max index of refraction
n_h = 1.3;    % min index of refraction
y = linspace(-h,h,23);    % plotting points of long width

%function n = ideal_refraction(n_o,n_h,y)
a = (1/h).*acosh(n_o/n_h);
n_y = n_o.*sech(a*y); % index of refraction along the GRIN Lens
figure(2)
stem(y,n_y);
