function [X,F] = MyDTFT(x,n)
F = linspace(-0.5,0.5,1e4);
basis = exp(-1i*2*pi*(n)'*F); %' Transposes and it outer multiplies
X = x*basis;

