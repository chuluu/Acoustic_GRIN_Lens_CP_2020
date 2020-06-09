%% Intro
% Plot Directivity based on COMSOL Data This is for both no lens and with
% lens plotting.

clear variables
close all
spy 
title('Fetching data!!!','Fontsize', 16);
xlabel('Bark! Bark!','Fontsize', 16);
pause(1)

%% Set up Data no lens for directivity
Data = csvread('Arc_Dots_w_o_lens.csv',5,0);
Pressure_no_Lens = Data(2:end);

a = 1;
b = 1;
c = 1;
for a = 1:length(Pressure_no_Lens)
   if mod(a,2) == 1
       Pressure_no_Lens_even(b) = Pressure_no_Lens(a);
       b = b + 1;
   else
       Pressure_no_Lens_odd(c) = Pressure_no_Lens(a);
       c = c + 1;
   end
end
Pressure_no_Lens_odd = flip(Pressure_no_Lens_odd);
Pressure_no_Lens = [Pressure_no_Lens_even, Pressure_no_Lens_odd];
Pressure_no_Lens = Pressure_no_Lens';

%% Set up Data with lens for directivity
Data = csvread('Directivity_with_Lens_Spkr_Center_-0_4_Dist_Away.csv',5,0);
Data_1_pressure = Data(2:end);

a = 1;
b = 1;
c = 1;
for a = 1:length(Data_1_pressure)
   if mod(a,2) == 1
       Data_1_pressure_even(b) = Data_1_pressure(a);
       b = b + 1;
   else
       Data_1_pressure_odd(c) = Data_1_pressure(a);
       c = c + 1;
   end
end
Data_1_pressure_odd = flip(Data_1_pressure_odd);
Data_1_pressure = [Data_1_pressure_even, Data_1_pressure_odd];
Data_1_pressure = Data_1_pressure';

Omnidirectional = ones(1,37);
Pressure_Magnitude = horzcat(Pressure_no_Lens,Data_1_pressure);
delta_angle = 5;
%% Processing Setup Angles
Angles = ((pi*5/180).*(0:1:36/2))';
Angles_pos = Angles';
Angles_2 = ((-pi*5/180).*(0:1:36/2))';
Angles_neg = fliplr(Angles_2');
Angles_total = [Angles_neg(1:18),Angles_pos];


Angles = ((pi*5/180).*(0:1:36))';
Angles_deg = ((5).*(0:1:36))';
Angles_rad = deg2rad(Angles_deg);

% %% Omindirectional case check
% polarplot(Angles,log10(Omnidirectional)); hold on
% Constant = (4*pi*(1^2)*57.3)/(2*pi);
% Last_value = 0;
% for a = 1:length(Angles_deg)
%     sum = Last_value + (1^2)*sin(Angles_rad(a))*delta_angle;
%     Last_value = sum;
% end
% sum = abs(sum);
% Q_f = Constant/sum;
% DI = 10*log10(Q_f); %% importdata
% 
% disp(['DI for OmniDirectional: ', num2str(DI), ' dB']);
%% Loop Through plot setup
[dummy,Meas_num] = size(Pressure_Magnitude);
delta_angle = 5;

%% Data Frequency
for data_set = 1:Meas_num
    Pressure_Magnitude_Dataset = Pressure_Magnitude(:,data_set);
    max_pressure = max(Pressure_Magnitude_Dataset);
%     if (max_pressure ~= Pressure_Magnitude(19,data_set))
%         disp('Check');
%     end

    % DI Utilized Beranek Acoustic book:
    Constant = (4*pi*(max_pressure^2)*57.3)/(2*pi);
    Last_value = 0;
    for a = 1:length(Angles_deg)
        sum = Last_value + (Pressure_Magnitude_Dataset(a)^2)*sin(Angles_rad(a))*delta_angle;
        Last_value = sum;
    end
    sum = abs(sum);
    Q_f = Constant/sum;
    DI(data_set) = 10*log10(Q_f);
    
    % Set to dB
    Normalized_dB = 20.*log(Pressure_Magnitude_Dataset./max_pressure);
    max_val = max(Normalized_dB);
    min_val = min(Normalized_dB);
    wiggle = max_val - min_val;
    radius_limits = [min_val-.1*wiggle, max_val+.1*wiggle];
    Normalized_dB = Normalized_dB';
    
    % Plotting Stuff
    polarplot(Angles,Normalized_dB,'Linewidth',1.2); hold on
    title('My Polar Plot');
    rlim(radius_limits);
end
disp([' ']);
disp(['DI for W/o Lens: ', num2str(DI(1)), ' dB']);
disp([' ']);
disp(['DI for W/ Lens: ', num2str(DI(2)), ' dB']);


legend('No Lens', 'W/ Lens');
title('Directivity Test With and W/o Lens'); 
