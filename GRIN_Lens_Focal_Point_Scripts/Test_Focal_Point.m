%% Program to do calculations of metamaterial grin lens design
%n_o = center at y = 0
%n_h = outermost index of refraction
%a = gradient parameter
%h = how tall the device is 

% clc 
% clear all

%% Presets that I chose based on this design
freq = 3000;
c= 343;
lamda = c/freq;

h = (2*lamda)/2; % 0.1143
d = 0.31*lamda;  % 0.0.0354
n_o = 2;
n_h = 1;

x = 0:0.01:d;
y = linspace(-h,h,50);    %plotting points of long width

%% Index of refraction Secant Calculations
a = (1/h).*acosh(n_o/n_h);
n_y = index_of_refraction_calculation(y,a,n_o);

%% Solving for focal point Setup
[H_a_d, H_f_d, Hd_a_d, Hd_f_d] = H_Parameters(d,a);

%% Beam Trajectory
% Hyp Transforms
y_o = 0:0.01:0.12;
for lm = 1:length(y_o)
    u_o(lm) = sinh(a.*y_o(lm));
    ud_o(lm) = a.*cosh(a.*y_o(lm));

    % u coordinate transformation 
    u_d_val = [H_f_d, H_a_d; Hd_f_d, Hd_a_d]*[u_o(lm); ud_o(lm)];
    u_d(lm) = u_d_val(1);
    ud_d(lm) = u_d_val(2);

    % beam trajectory
    y_d(lm) = (1/a)*asinh(u_d(lm));
    yd_d(lm) = (1/(a*cosh(a*y_d(lm))))*ud_d(lm);
    n_y_d(lm) = index_of_refraction_calculation(y_d(lm),a,n_o);
    A = repmat(n_y_d(lm),[1 length(n_y)]);
    [minValue,closestIndex] = min(abs(A-n_y));
    
    if n_y(closestIndex) < n_y_d(lm)
        if closestIndex == length(y)-1
            closestIndex = closestIndex - 1;
        else 
            closestIndex = closestIndex + 1;
        end
        n_y_discretize = n_y(closestIndex);
    else
        n_y_discretize = n_y(closestIndex);
    end
    

    Num_x_f(lm) = 1-((yd_d(lm)^2)*((n_y_discretize^2) - 1)); 
    Den_x_f(lm) = n_y_discretize^2;

    x_f(lm) = d-real((y_d(lm)/yd_d(lm))*sqrt(Num_x_f(lm)/Den_x_f(lm)));
end

%% plotting
plot(y_o,-x_f); hold on;

%%
y = linspace(-h,h,20);    %plotting points of long width
y_ideal = linspace(-h,h,500);    %plotting points of long width

a = (1/h).*acosh(n_o/n_h);
n_y = index_of_refraction_calculation(y,a,n_o);
n_y_ideal = index_of_refraction_calculation(y_ideal,a,n_o);
stem(y,n_y); hold on
plot(y_ideal,n_y_ideal);
ylim([1 2]);