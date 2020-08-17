function [ft] = hwfft2(xr,yr,k,l,f)
%[ft] = hwfft2(xr,yr,k,l,f)
% compared with obfft2: 
% (1) k,l doesn't change within the main program any more.
% (2) input ft is CENTERED at (Nx/2,Ny/2);
%     output f is CENTERED at (Nx/2,Ny/2) 
%     (replacing circumshift with fftshift/ifftshift)
%
% Fast Fourier Transform of (k,l,ft) into f.  
% 2D Fast Fourier Transform of (x,y,f) into ft; k=k_of_x(x),l=l_of_y(y). 
% The length of 
% x,y  and f must be an even number, preferably a power of two.  The index of
% the zero mode is inull=jnull=N/2.
Nx         = max(size(xr));
Ny         = max(size(yr));
Periodx    = Nx*(xr(2)-xr(1));
Periody    = Ny*(yr(2)-yr(1));
inull      = Nx/2;
jnull      = Ny/2;
ft         = (Periodx/Nx)*(Periody/Ny)*circshift(fft2(circshift(f,-[inull-1 jnull-1])),[inull-1 jnull-1]);
%ft         = (Periodx/Nx)*(Periody/Ny)*fftshift(fft2(ifftshift(f)));


