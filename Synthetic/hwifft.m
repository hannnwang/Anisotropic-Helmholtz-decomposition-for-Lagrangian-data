function [f] = hwifft(xr,k,ft,rflag)
%[f] = hwifft(xr,k,ft,rflag)
% Fast Fourier INVERSE Transform of the pair (k,ft) into (x,f).  The length of 
% x and f must be an even number, preferably a power of two.  
% If rflag is set then only the real part of f is returned.

%20171204: IMPORTANT: corrected the circshift operation.
N         = max(size(k));
Period    = N*(xr(2)-xr(1));
inull     = N/2;

f         = (N/Period)*circshift(ifft(circshift(ft,-(inull-1))),(inull-1));

%f         = (N/Period)*ifftshift(ifft(fftshift(ft)));
%f         = (N/Period)*circshift(ifft(circshift(ft,-[0 inull-1])),[0 inull-1]);

%f         = (N/Period)*ifftshift(ifft(fftshift(ft)));
if (nargin == 4)
    if (rflag == 1)
        if (max(abs(imag(f)))==0)
            f = real(f);
        else
            f = real(f);
%             warning('the inverse transform has nonzero imaginary parts. The imaginary part IS TRUNCATED.')
        end
    end
end
