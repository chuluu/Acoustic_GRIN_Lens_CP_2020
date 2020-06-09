function [hn, HF, F] = FIR_Filter_By_Freq_Sample(HF_mag_samples, figurenum)
% hn - impulse response of filter (same length as HF_mag_samples)
% HF - complex frequency response of filter
% (estimated H(F) values found by FFT or freqz)
% F ? digital frequency values corresponding to the estimated H(F)values
% HF_mag_samples ? H[k] Magnitude response samples for desired filter
% figurenum - Figure # to plot frequency responses
    M  = length(HF_mag_samples);
    dF = 1/M;
    F  = 0:dF:(M-1)/M;
    
    if mod(M, 2) == 0
       disp('provide an odd number of samples')
       
       hn = [];
       HF = [];
       F  = [];
       
       return;
    end
    
    % compute phase at each sample
    phase = zeros(M, 0);
    
    for k = 0:M-1
        phase(k+1) = exp(-j*pi*k*(M-1)/M);
    end
    
    % combine magnitude and phase components of frequency response
    HF = HF_mag_samples .* phase;
    
    hn = real(ifft(HF));
    
    % higher resolution FFT for plot
    M_hr  = 512;
    HF_hr = fft(hn, M_hr);
    F_hr  = 0:1/M_hr:(M_hr-1)/M_hr;
    
 
    figure(figurenum);
    subplot(2, 1, 1);
    title('Frequency Response');
    ylabel('Magnitude');
    xlabel('Digital Frequency (cyc/s)');
    hold on;
    stem(F, abs(HF));
    plot(F_hr, abs(HF_hr));
    
    subplot(2, 1, 2);
    ylabel('Phase');
    xlabel('Digital Frequency (cyc/s)');
    hold on;
    stem(F, angle(HF)/pi);
    plot(F_hr, angle(HF_hr)/pi);
    
    mag_dB = 20*log10(abs(HF_hr));
    
    figure(figurenum+1);
    
    plot(F_hr, mag_dB);
    title('Magnitude Response');
    ylabel('Magnitude (dB)');
    xlabel('Digital Frequency (cyc/s)');
    
 
end
