function [ft] = hwfft(xr,k,f)
% Fast Fourier Transform of the pair (x,f) into (k,ft).  The length of 
% x and f must be an even number, preferably a power of two.  The index of
% the zero mode is inull=N/2.

%20171204: IMPORTANT: corrected the circshift operation.

N         = max(size(xr));
Period    = N*(xr(2)-xr(1));
inull     = N/2;
ft        = (Period/N)*circshift(fft(circshift(f,-(inull-1))),(inull-1));


%ft        = (Period/N)*fftshift(fft(ifftshift(f)));


%ft        = (Period/N)*circshift(circshift(fft(f),-[0 inull-1]),[0 inull-1]);



% ft         = (Periodx/Nx)*(Periody/Ny)*circshift(fft2(circshift(f,-[inull-1 jnull-1])),[inull-1 jnull-1]);

