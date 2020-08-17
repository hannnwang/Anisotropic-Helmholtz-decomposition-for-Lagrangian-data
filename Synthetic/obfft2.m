function [k , l, ft] = obfft2(x,y,f)
% 2D Fast Fourier Transform of (x,y,f) into (k,l,ft).  The length of 
% x,y  and f must be an even number, preferably a power of two.  The index of
% the zero mode is inull=jnull=N/2.
Nx         = max(size(x));
Ny         = max(size(y));
k          = k_of_x(x);
l          = k_of_x(y);
Periodx    = Nx*(x(2)-x(1));
Periody    = Ny*(y(2)-y(1));
inull      = Nx/2;
jnull      = Ny/2;
ft         = (Periodx/Nx)*(Periody/Ny)*circshift(fft2(f),[inull-1 jnull-1]);
