function [ft] = hwfft3(xr,yr,zr,k,l,m,f)
Nx         = max(size(xr));
Ny         = max(size(yr));
Nz=max(size(zr));
Periodx    = Nx*(xr(2)-xr(1));
Periody    = Ny*(yr(2)-yr(1));
Periodz    = Nz*(zr(2)-zr(1));
inull      = Nx/2;
jnull      = Ny/2;
knull      = Nz/2;
ft         = (Periodx/Nx)*(Periody/Ny)*(Periodz/Nz)*circshift(fftn(circshift(f,-[inull-1 jnull-1 knull-1])),[inull-1 jnull-1 knull-1]);


