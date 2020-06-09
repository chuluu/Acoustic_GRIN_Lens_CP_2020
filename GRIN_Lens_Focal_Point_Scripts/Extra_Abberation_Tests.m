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
%%
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


%% Abberation and intial condition focal point dependency testing
figure(2)
B2 = linspace(-0.000002,-0.25,4);
B1 = linspace(-0.0001,-1,4);%-0.0679;
for b = 1:length(B1)
    for e = 1:length(B2)
        for c = 1:length(yo)
            yd = (1./a).*asinh(sinh(a.*yo(c)).*cos(a.*d));
            yo_x_f = yo(c);

            n_y_d = index_of_refraction_calculation_yo_independent(yd,a,n_o,B1(b),B2(e));
            n_yo  = index_of_refraction_calculation_yo_independent(yo_x_f,a,n_o,B1(b),B2(e));

            x_f_d(c) = real(yd*sqrt((1/((n_y_d^2)-(n_yo^2))) - 1));
        end            
        diff_x_f = diff(x_f_d);
        diff_x_f = diff_x_f(2:end);
        avg_diff = sum(abs(diff_x_f))/length(diff_x_f);
        if (-0.05 < avg_diff & avg_diff < 0.05)
            plot(yo,x_f_d,'Linewidth',1.6); hold on
            diff_x_f;
        end
    end
end
title('focal length vs. y_o (aberration secant change)','Fontsize',14);
xlabel('y_o (m)','Fontsize',14);
ylabel('focal length (m)','Fontsize',14);
% %% Calculate focal length
% 
% for b = 1:length(yo)
%     y_d = (1./a).*asinh(sinh(a.*yo(b)).*cos(a.*d));
%     yd_d_num = -sinh(a*yo(b))*a*sin(a*d);
%     yd_d_den = a*cosh(asinh(sinh(a.*yo(b)).*cos(a.*d)));
%     yd_d = yd_d_num/yd_d_den;
%     n_y_d = index_of_refraction_calculation(y_d,a,n_o);
%     x_f_num = 1 - (yd_d^2)*((n_y_d^2) - 1);
%     x_f_den = n_y_d^2;
% 
%     x_f(b) = -(y_d/yd_d)*sqrt(x_f_num/x_f_den);
% end
% plot(yo,x_f)
% 
% %% Calc
% for b = 1:length(yo)
%     y_d = (1./a).*asinh(sinh(a.*yo(b)).*cos(a.*d));
%     B = sinh(a*yo(b));
%     yd_d_num = -a*sin(a*d);
%     yd_d_den = sqrt(((B^2)*(cos(a*d).^2)) + 1);
%     yd_d = yd_d_num/yd_d_den;
%     n_y_d = index_of_refraction_calculation(y_d,a,n_o);
%     x_f_num = 1 - (yd_d^2)*((n_y_d^2) - 1);
%     x_f_den = n_y_d^2;
% 
%     x_f(b) = -(y_d/yd_d)*sqrt(x_f_num/x_f_den);
% end
% plot(yo,real(x_f))
% 
% 
% 
% 
% %%
% figure(3);
% for b = 1:length(yo)
%     y = (1./a).*asinh(sinh(a.*yo(b)).*cos(a.*x) + cosh(a.*yo(b))*sin(a*x));
%     plot(x,y,'Linewidth',1.6); hold on;
%     plot(x,-y,'Linewidth',1.6);
% end
% rectangle('Position',[-d -h d 2*h])
% title('y-axis vs. x-axis (Beam Trajectory of Lens Problem)','Fontsize',14);
% xlabel('x-axis (m)','Fontsize',14);
% ylabel('y-axis (m)','Fontsize',14);
% grid on
