function [k , ft] = obffta(x,f)
% Zero-padded with 3/2 rule
% Fast Fourier Transform of the pair (x,f) into (k,ft).  The length of 
% x and f must be an even number, preferably a power of two.  The index of
% the zero mode is inull=N/2.
N         = max(size(x));
Nb        = 3*N/2;
dxb       = (x(2)-x(1))*2/3;
xb        = x(1) + dxb*(0:Nb-1);
k         = k_of_x(xb);
Period    = N*(x(2)-x(1));
inullb     = Nb/2;
ft        = (Period/Nb)*circshift(fft(f,Nb),[0 inullb-1]);
