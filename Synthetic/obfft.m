function [k , ft] = obfft(x,f)
% Fast Fourier Transform of the pair (x,f) into (k,ft).  The length of 
% x and f must be an even number, preferably a power of two.  The index of
% the zero mode is inull=N/2.
N         = max(size(x));
k         = k_of_x(x);
Period    = N*(x(2)-x(1));
inull     = N/2;
ft        = (Period/N)*circshift(fft(f),[0 inull-1]);
