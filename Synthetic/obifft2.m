function [x , y, f] = obifft2(k,l,ft,rflag)
% Fast Fourier INVERSE Transform of (k,l,ft) into (x,y,f).  The length of 
% x,y and f must be an even number, preferably a power of two.  
% x,y start with zero.
% If rflag is set then only the real part of f is returned.
Nx         = max(size(k));
Ny         = max(size(l));
x         = x_of_k(k);
y         = x_of_k(l);
Periodx    = Nx*(x(2)-x(1));
Periody    = Ny*(y(2)-y(1));
inull     = Nx/2;
jnull     = Ny/2;
f         = (Nx/Periodx)*(Ny/Periody)*ifft2(circshift(ft,-[inull-1 jnull-1]));
if (nargin == 4)
    if (rflag == '1')
        f = real(f);
    end
end
