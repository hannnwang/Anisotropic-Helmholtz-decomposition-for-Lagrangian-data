function [f] = hwifft3(xr,yr,zr,k,l,m,ft,rflag)
% Fast Fourier INVERSE Transform of (k,l,m,ft) into (x,y,z,f).  The length of 
% x,y,z and f must be an even number, preferably a power of two.  
% x,y,z start with zero.
% If rflag is set then only the real part of f is returned.
Nx         = max(size(k));
Ny         = max(size(l));
Nz=max(size(m));
Periodx    = Nx*(xr(2)-xr(1));
Periody    = Ny*(yr(2)-yr(1));
Periodz=Nz*(zr(2)-zr(1));
inull     = Nx/2;
jnull     = Ny/2;
knull     = Nz/2;

f         = (Nx/Periodx)*(Ny/Periody)*(Nz/Periodz)*circshift(ifftn(circshift(ft,-[inull-1 jnull-1 knull-1])),[inull-1 jnull-1 knull-1]);
if (nargin == 8)
    if (rflag == '1')
        if (max(max(max(abs(imag(f)))))==0)
            f = real(f);
        else
            warning('the inverse transform has nonzero imaginary parts. The imaginary part is then NOT truncated.')
        end
    end
end