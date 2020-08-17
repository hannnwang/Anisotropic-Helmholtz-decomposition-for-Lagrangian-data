function [x , f] = obifft(k,ft,rflag)
% Fast Fourier INVERSE Transform of the pair (k,ft) into (x,f).  The length of 
% x and f must be an even number, preferably a power of two.  
% x starts with zero.
% If rflag is set then only the real part of f is returned.
N         = max(size(k));
x         = x_of_k(k);
Period    = N*(x(2)-x(1));
inull     = N/2;
f         = (N/Period)*ifft(circshift(ft,-[0 inull-1]));
if (nargin == 3)
    if (rflag == '1')
        f = real(f);
    end
end
