function plot_freq_responses(Fd, HF, fsample, figure_num)
    % INPUTS:
    % Frequency
    % Data
    % Sampling Rate
    % Figure number
    % OUTPUTS:
    % Plotted FR
    figure(figure_num);
    subplot(2,1,1); plot(Fd,abs(HF)); hold on %New Changes 
    % Plot the magnitude of HF on a linear scale
    grid on
    xlabel('Digital Frequency  F (cycles/sample)')
    ylabel('Magnitude Response')
    title('Frequency Response of Filter');
    % Phase
    subplot(2,1,2); plot(Fd,angle(HF)/pi); %New Changes
    grid on
    xlabel('Digital Frequency  F (cycles/sample)')
    ylabel('Phase Response /pi')
    
    figure(figure_num+1);
    f_analog = Fd.*fsample;
    subplot(2,1,1); plot((f_analog),20*log10(abs(HF))); %New Changes
    xlim([min(f_analog) max(f_analog)]);
    % Plot the magnitude of HF on a linear scale
    grid on
    xlabel('Analog Frequency  f (Hz)')
    ylabel('Magnitude Response (dB)')
    title('Frequency Response of Filter');
    % Phase
    subplot(2,1,2); plot((f_analog),angle(HF)/pi); %New Changes
    xlim([min(f_analog) max(f_analog)]);
    grid on
    xlabel('Analog Frequency  f (Hz)')
    ylabel('Phase Response /pi')

end

