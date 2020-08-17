function [k , l, m, ft] = obfft3(x,y,z,f)
% 3D Fast Fourier Transform of (x,y,z,f) into (k,l,m,ft).  The length of 
% x,y,z  and f must be an even number, preferably a power of two.  The index of
% the zero mode is inull=jnull=N/2.
Nx         = max(size(x));
Ny         = max(size(y));
Nz         = max(size(z));
k          = k_of_x(x);
l          = k_of_x(y);
m          = k_of_x(z);
Periodx    = Nx*(x(2)-x(1));
Periody    = Ny*(y(2)-y(1));
Periodz    = Nz*(z(2)-z(1));
inull      = Nx/2;
jnull      = Ny/2;
knull      = Nz/2;
ft         = (Periodx/Nx)*(Periody/Ny)*(Periodz/Nz)*...
    circshift(fftn(f),[inull-1 jnull-1 knull-1]);
