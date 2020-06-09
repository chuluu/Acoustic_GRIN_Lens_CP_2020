%% importdata
Data_1 = importdata('1000Hz_Max_dB_Data_Reflections_200519.csv');
Data_1 = Data_1(1:37);
% Data_1 = Data_1';
Data_1_pressure = (10.^(Data_1/20)) * (2.*10.^-5);

% Data_2 = importdata('2500Hz_Max_dB_Data_Reflections_191221.csv');
% Data_2 = Data_2(:,1);
% Data_2_pressure = (10.^(Data_2/20)) * (2.*10.^-5);
% Data_3 = importdata('190925_4000Hz_Setup_Directivity.csv',',');
% Data_3 = Data_3(:,1);
% Data_3_pressure = (10.^(Data_3/20)) * (2.*10.^-5);
Omindirectional = ones(1,37);
%% Setup Arrays for data manipulation
Pressure_Magnitude = horzcat(Data_1_pressure);%,Data_2_pressure);%Data_3_pressure); % add the data using comma etc.
Angles = ((pi*5/180).*(0:1:36/2))';
Angles_pos = Angles';
Angles_2 = ((-pi*5/180).*(0:1:36/2))';
Angles_neg = fliplr(Angles_2');
Angles_total = [Angles_neg(1:18),Angles_pos];


Angles = ((pi*5/180).*(0:1:36))';
Angles_deg = ((5).*(0:1:36))';
Angles_rad = deg2rad(Angles_deg);
%% Loop Through plot many
[dummy,Meas_num] = size(Pressure_Magnitude);
HPBW_array = zeros(3,1);
delta_angle = 5;

%%
for data_set = 1:Meas_num
    Pressure_Magnitude_Dataset = Pressure_Magnitude(:,data_set);
    max_pressure = max(Pressure_Magnitude_Dataset);
    if (max_pressure ~= Pressure_Magnitude(19,data_set))
        disp('Check');
    end

    % DI Utilized Beranek Acoustic book:
    Constant = (4*pi*(max_pressure^2)*57.3)/(2*pi);
    Last_value = 0;
    for a = 1:length(Angles_deg)
        sum = Last_value + (Pressure_Magnitude_Dataset(a)^2)*sin(Angles_rad(a))*delta_angle;
        Last_value = sum;
    end
    sum = abs(sum);
    Q_f = Constant/sum;
    DI = 10*log10(Q_f)
    
    % Set to dB
    Normalized_dB = 20.*log(Pressure_Magnitude_Dataset./max_pressure);
    max_val = max(Normalized_dB);
    min_val = min(Normalized_dB);
    wiggle = max_val - min_val;
    radius_limits = [min_val-.1*wiggle, max_val+.1*wiggle];
    Normalized_dB = Normalized_dB';
    
    
    % HPBW:
%     dB_down = -6;
%     Index = find((dB_down+(dB_down*0.15)) < Normalized_dB & (dB_down-(dB_down*0.15)) > Normalized_dB);
%     Angles_dB_Down = Angles_deg(Index);
%     if (Angles_dB_Down(2) == Angles_dB_Down(1) + 5)
%         HPBW = Angles_dB_Down(3) - Angles_dB_Down(1);
%     else 
%         HPBW = Angles_dB_Down(2) - Angles_dB_Down(1);
%     end
%     HPBW_array(data_set) = HPBW;
    
    
    % Plotting Stuff
    polarplot(Angles,Normalized_dB); hold on
    title('My Polar Plot');
    rlim(radius_limits);
end

legend('1000 Hz');
title('Directivity Baseline'); 

% %% Omindirectional case check
% polarplot(Angles,Omindirectional);
% Constant = (4*pi*(1^2)*57.3)/(2*pi);
% Last_value = 0;
% for a = 1:length(Angles_deg)
%     sum = Last_value + (1^2)*sin(Angles_rad(a))*delta_angle;
%     Last_value = sum;
% end
% sum = abs(sum);
% Q_f = Constant/sum;
% DI = 10*log10(Q_f)