function [x , y, z, f] = obifft3(k,l,m,ft,rflag)
% Fast Fourier INVERSE Transform of (k,l,m,ft) into (x,y,z,f).  The length of 
% x,y,z and f must be an even number, preferably a power of two.  
% x,y,z start with zero.
% If rflag is set then only the real part of f is returned.
Nx         = max(size(k));
Ny         = max(size(l));
Nz         = max(size(m));
x         = x_of_k(k);
y         = x_of_k(l);
z         = x_of_k(m);
Periodx    = Nx*(x(2)-x(1));
Periody    = Ny*(y(2)-y(1));
Periodz    = Nz*(z(2)-z(1));
inull     = Nx/2;
jnull     = Ny/2;
knull     = Nz/2;
f         = (Nx/Periodx)*(Ny/Periody)*(Nz/Periodz)*...
    ifftn(circshift(ft,-[inull-1 jnull-1 knull-1]));
if (nargin == 5)
    if (rflag == '1')
        f = real(f);
    end
end
