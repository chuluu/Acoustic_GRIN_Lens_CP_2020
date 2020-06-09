function [h,t,H,f] = MyFFTDeConv(x,y,fs)
% INPUTS:
% x = input signal
% y = output signal
% fs = sampling rate
% OUTPUTS:
% h = impulse response
% t = time 
% H = Deconvolved transfer function
% f = frequency?

x = x';
N_x = length(x); 
N_y = length(y);
N_sum = N_x + N_y - 1;
N_sum_pow2 = nextpow2(N_sum);
N_sum_pow2 = 2.^N_sum_pow2;
N = N_sum_pow2 - 1;

%"Convolve" Multiply frequency domain signals
X = fft(x,N);
for a = 1:length(X)
    if X(a) <0.1
        X(a) = 0.1;
    end
end

Y = fft(y,N);
H = Y./X;


%Inverse fourier transform and obtain time domain output
h = real(ifft(H))/N;
dt = 1/fs;
t = 0:dt:(length(h)*dt)-dt;
if rem(N,2) == 0
    % even number of samples.
    f = fs*((N/-2):(-1+N/2))/N;
else
    % odd number of samples.
    f = fs*(((N-1)/-2):((N-1)/2))/N;
end
end