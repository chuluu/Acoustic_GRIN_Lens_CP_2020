function [Cxy,t] = NormCrossCorrelate(xn, yn, fs)
% Convolvles 2 signals together using FFT algorithm
%Inputs:
%x:    Data Array of first time domain signal
%h:    Data Array of second time domain signal
%fs:   Sampling rate
%figure_num:   label figures generated
%Outputs:
%y:    Convolved time domain signal
%t:    Time samples based off of sampling rate
%Figure Outputs:
%figure_num: time domain plots of x,h, and y
%figure_num + 1:  Digital Frequency Spectrum of X, H, and Y
%figure_num + 2:  Analog frequency Spectrum of X, H, and Y

%Set up length of FFT to zero pad such that we get an even signal
%throughout
N_x = length(xn); 
N_h = length(yn);
N_sum = N_x + N_h - 1;
N_sum_pow2 = nextpow2(N_sum);
N_sum_pow2 = 2.^N_sum_pow2;
N = N_sum_pow2 - 1;

%"Convolve" Multiply frequency domain signals
X = fft(xn,N);
Y = fft(yn,N);
Y = conj(Y);
[X_m,X_n] = size(X); %Used to check for inccorect matrix multiplication
[Y_m,Y_n] = size(Y); %Used to check for inccorect matrix multiplication

if (Y_n == X_n)  %Used to check for inccorect matrix multiplication
    Y = Y;
else
    Y = Y';
end

% AutoCorrelation Section:
XX = X.*conj(X);
YY = Y.*conj(Y);

CXY = Y.*X;
[CXY_m,CXY_n] = size(CXY);

if (CXY_n > 1)  %Used to properly have the right array length
    CXY = CXY';
end

%Inverse fourier transform and obtain time domain output
Cxy = real(ifft(CXY));
xx = real(ifft(XX));
yy = real(ifft(YY));
Cxy = Cxy./(sqrt(xx(1)*yy(1)));
n_Cxy = 0:1:length(Cxy)-1;
dt = 1/fs;
t = (0:1:length(Cxy)-1).*dt;  %Time domain from sampling rate
end
