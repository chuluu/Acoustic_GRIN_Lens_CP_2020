%% Intro
clear variables
close all
% spy 
% title('Fetching data','Fontsize', 16);
% xlabel('Bark! Bark!','Fontsize', 16);

%% Choosing Which to run

disp(['What program do you want to run: ']);
prompt = ['Directivity_Measurement_Single_Tone_Method (1), Frequency_Response_THD_Spkr (2), Microphone_Calibration (3), else (4): '];
choice = input(prompt);

if (choice == 1)
    run('Directivity_Measurement_Single_Tone_Method.m');
elseif (choice == 2)
    run('Frequency_Response_THD_Spkr.m');
elseif (choice == 3)
    run('Microphone_Calibration.m');
elseif (choice == 4)
    disp(['There are other programs to run, but some are buggy: ']);
    prompt_2 = ['Directivity_Measurement_Farina_Method (1), Single_Tone_Test (2)'];
    choice_2 = input(prompt_2);
    if (choice_2 == 1)
        run('Microphone_Calibration.m');
    elseif (choice_2 == 2)  
        run('Directivity_Measurement_Farina_Method.m');
    end
else
    disp(['Enter valid program']);
end






