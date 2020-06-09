%% 

r = 1;
Arc_Center = -0.2238;
theta = 0:5:180;
theta_rad = deg2rad(theta);
x = r.*cos(theta_rad) + Arc_Center;
y = r.*sin(theta_rad);
x = x';
y = y';
plot(x,y)