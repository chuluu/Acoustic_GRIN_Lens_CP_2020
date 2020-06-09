%% Program to do calculations of metamaterial grin lens design
% Program will show the convergence of discrete rows to continuous pattern
% for the focal point.
clc 
clear
x = figure(1);
spy 
title('He is fetching data','Fontsize', 16);

%% Inputs 
freq = 2300;     % Frequency
c= 343;          % Speed of Sound
lamda = c/freq;  % wavelength
h = 0.1625; % 0.1143
h_real = h*2;
d = 0.052;  % 0.0.0354
n_o = 2;%1.965;         % Index of refraction profile high
n_h = 1;%1.007;         % Index of refraction profile low
y_o = 0.21*h;  % Initial condition

%% Solve beam trajectory and focal point for row discretization
b_initial = 5;   % Initial starting point for b for loop
b_step = 1;      % Input step of b
max_row = 300;   % Max amount of steps

b = b_initial;   % Initial starting point for array to work
c = 1;           % A needed variable for array filling
while (b <= max_row)
    % Solve discretized n_y
    y = linspace(-h,h,b);    %plotting points of long width

    a = (1/h).*acosh(n_o/n_h);
    n_y = index_of_refraction_calculation(y,a,n_o);

    [H_a_d, H_f_d, Hd_a_d, Hd_f_d] = H_Parameters(d,a);

    % Hyp Transforms

    u_o = sinh(a.*y_o);
    ud_o = 0;
    alpha = a.*cosh(a.*y_o);

    % u coordinate transformation 
    u_d_val = [H_f_d, H_a_d; Hd_f_d, Hd_a_d]*[u_o; ud_o];
    u_d = u_d_val(1);
    ud_d = u_d_val(2);

    % beam trajectory
    y_d = (1/a)*asinh(u_d);
    yd_d = (1/(a*cosh(a*y_d)))*ud_d;
    
    % calculate needed n_y_d
    n_y_d = index_of_refraction_calculation(y_d,a,n_o);
    
    % Closest value of n_y_d with discretized secant pattern
    A = repmat(n_y_d,[1 length(n_y)]);
    [minValue,closestIndex] = min(abs(A-n_y));
    n_y_discretize = n_y(closestIndex);
    
    % Calculated focal point
    Num_x_f = 1-((yd_d^2)*((n_y_discretize^2) - 1)); 
    Den_x_f = n_y_discretize^2;

    x_f(c) = -real((y_d/yd_d)*sqrt(Num_x_f/Den_x_f));
    b = b + b_step;
    c = c + 1;
end

%% Determine Focal Point
Num_x_f = 1-((yd_d^2)*((n_y_d^2) - 1)); 
Den_x_f = n_y_d^2;
x_f_ideal = -real((y_d/yd_d)*sqrt(Num_x_f/Den_x_f));
number_of_rows = b_initial:b_step:max_row;

y = linspace(-h,h,25);            % plotting points of long width
y_ideal = linspace(-h,h,3000);    % plotting points of long width

a = (1/h).*acosh(n_o/n_h);
n_y = index_of_refraction_calculation(y,a,n_o);
n_y_ideal = index_of_refraction_calculation(y_ideal,a,n_o);

x_f_ideal*100

%% Percent difference
P_E = percent_error(x_f,x_f_ideal*ones(1,length(number_of_rows)));

P_E = ((x_f_ideal - x_f)/x_f_ideal)*100;
P_E(21)
find(P_E < 1.5);
number_of_rows(21)
%% Plotting
figure(1);
C = linspecer(2);
axes('NextPlot','replacechildren', 'ColorOrder',C);
plot(number_of_rows,x_f,'.-','Linewidth',2.5); hold on
plot(number_of_rows,x_f_ideal*ones(1,length(number_of_rows)),'Linewidth',2.5); 
title('xf vs. number of GRIN lens rows','Fontsize',14);
xlabel('Number of Rows [Discrete Points]','Fontsize',14);
ylabel('Focal Point (m)','Fontsize',14);
legend('Discretization','Continuous');
grid on;

figure(2);
C = linspecer(2);
axes('NextPlot','replacechildren', 'ColorOrder',C);
stem(y,n_y,'Linewidth',2.5); hold on
plot(y_ideal,n_y_ideal,'Linewidth',2.5); 
title('Discretized Secant Pattern','Fontsize',14);
xlabel('y-axis (transverse) (m)','Fontsize',14);
ylabel('index of refraction','Fontsize',14);
legend('Discretization','Continuous');
ylim([1 2]);
grid on;