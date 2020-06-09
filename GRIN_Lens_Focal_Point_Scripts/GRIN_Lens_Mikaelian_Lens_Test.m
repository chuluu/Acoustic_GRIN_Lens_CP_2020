%% 
freq = 2000;     % Frequency
c= 343;          % Speed of Sound
lamda = c/freq;  % wavelength
k = (2*pi)/lamda;
h = 0.1625; % 0.1143
d = 0.25*lamda%0.052  % 0.0.0354
n_o = 3;%1.965;         % Index of refraction profile high
n_h = 1;%1.007;         % Index of refraction profile low
yo = [0,h*0.5,h];  % Initial condition
a = (1/h).*acosh(n_o/n_h);
x = 0:0.01:d*8;
%%
figure;
for b = 1:length(yo)
    y = (1./a).*asinh(sinh(a.*yo(b)).*cos(a.*x));
    plot(x,y); hold on;
    plot(x,-y);
end
plot(d*ones(1,15),linspace(-h, h, 15));
rectangle('Position',[0 -h d 2*h])

%%
figure;
for b = 1:length(yo)
    y = (1./a).*asinh(sinh(a.*yo(b)).*cos(a.*x) + cosh(a.*yo(b))*sin(a*x));
    plot(x,y); hold on;
    plot(x,-y);
end
plot(d*ones(1,15),linspace(-h, h, 15));
rectangle('Position',[0 -h d 2*h])
