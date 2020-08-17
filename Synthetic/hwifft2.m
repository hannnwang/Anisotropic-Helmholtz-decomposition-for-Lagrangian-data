function [f] = hwifft2(xr,yr,k,l,ft,rflag)
%[f] = hwifft2(xr,yr,k,l,ft,rflag)
% compared with obifft2: 
% (1) x,y doesn't change within the main program any more.
% (2) input ft is CENTERED at (Nx/2,Ny/2);
%     output f is CENTERED at (Nx/2,Ny/2) 
%     (replacing circumshift with fftshift/ifftshift)
%
% Fast Fourier INVERSE Transform of (k,l,ft) into f.  
% x=x_of_k(k),y=y_of_l(l).
% The length of 
% x,y and f must be an even number, preferably a power of two.  
% x,y start with zero.
% If rflag is set then only the real part of f is returned.
Nx         = max(size(k));
Ny         = max(size(l));
Periodx    = Nx*(xr(2)-xr(1));
Periody    = Ny*(yr(2)-yr(1));
inull     = Nx/2;
jnull     = Ny/2;

%f         = (Nx/Periodx)*(Ny/Periody)*ifftshift(ifft2(fftshift(ft)));

f         = (Nx/Periodx)*(Ny/Periody)*circshift(ifft2(circshift(ft,-[inull-1 jnull-1])),[inull-1 jnull-1]);
if (nargin == 6)
    if (rflag == '1')
        if (max(max(abs(imag(f))))==0)
            f = real(f);
        else
            warning('the inverse transform has nonzero imaginary parts. The imaginary part is then NOT truncated.')
        end
    end
end


% 2D Fast Fourier Transform of (x,y,f) into ft; k=k_of_x(x),l=l_of_y(y). 