function hn_lp = FIR_Filter_By_Window(M,Fc,window) 
% M = the filter length (odd)
% Fc = filter cutoff digital frequency (-6dB) (cycles/sample)
% window = the Matlab window function values to multiply h[n] by
% hn_lp = windowed impulse response values for the Low-pass FIR filter

% rectwin( ) - Rectangular window.
% triang( ) - Triangular window.
% bartlett( ) - Bartlett window (Triangle Window with ?0?s at both ends)
% hann( ) - von Hann window (Hanning Window with ?0?s at both ends)
% hanning( ) - Hanning window (With non-zero end samples).
% hamming ( ) - Hamming window.
% blackman( ) - Blackman window.
% kaiser( ) - Kaiser window (with b parameter)
% tukeywin( ) - Tukey window.
% Barthannwin( ) - Modified Bartlett-Hanning window.
% bohmanwin( ) - Bohman window.
% chebwin( ) - Chebyshev window.
% flattopwin( ) - Flat Top window.
% gausswin( ) - Gaussian window.
% blackmanharris( ) - Minimum 4-term Blackman-Harris window
% nuttallwin( ) - Nuttall defined minimum 4-term Blackman-Harris window.
% parzenwin( ) - Parzen (de la Valle-Poussin) window. 

n = -(M-1)/2:1:(M-1)/2;
h = 2.*Fc.*(sin(2*pi*Fc.*n)./(2*pi*Fc.*n));
h(n==0) = 2.*Fc;
w = window;
w = w';
hn_lp = (h.*w);
end