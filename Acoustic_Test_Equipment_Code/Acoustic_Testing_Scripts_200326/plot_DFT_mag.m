function [X,F,f] = plot_DFT_mag(x,fs,figure_num)
% INPUTS:
% x:    Data array.
% fs:   Sampling Rate (Hz).
% figure_num   to define the figure

% OUTPUTS:
% X:    Shifted and normalized FT of w.  Complex.  DC component located in
% the "middle".
% f:    Array of continuous-time frequencies corresponding to each
% component of W.

%Figure Plots
%subplot(1): Digital Responses of signal input
%Subplot(2): Analog Response using fft shift to have negative and positive
%analog frequency


%Supplemental help provided by Clay McKell

%Set up fft lengths
N = length(x);      % Normalize by sample length.
X = (fft(x))./N;
X_shift = fftshift(X);    % Put DC in the middle.
F = (0:1:N-1)./N;
if rem(N,2) == 0
    % even number of samples.
    f = fs*((N/-2):(-1+N/2))/N;
else
    % odd number of samples.
    f = fs*(((N-1)/-2):((N-1)/2))/N;
end

%Plotting 
figure(figure_num);
subplot(2,1,1); plot(F,abs(X));
xlabel('Digital Frequency (cyc/sample)'); ylabel('Magnitude'); 
title('FFT Spectrum Plot');

subplot(2,1,2); plot(f,abs(X_shift));
xlabel('Analog Frequency (Hz)'); ylabel('Magnitude'); 