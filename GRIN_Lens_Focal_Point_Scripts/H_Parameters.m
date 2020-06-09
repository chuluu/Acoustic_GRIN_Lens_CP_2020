function [H_a, H_f, Hd_a, Hd_f] = H_Parameters(x,a)
    H_a = sin(a.*x)./a;   
    H_f = cos(a.*x);
    Hd_a = cos(a.*x);
    Hd_f = -a.*sin(a.*x);
end