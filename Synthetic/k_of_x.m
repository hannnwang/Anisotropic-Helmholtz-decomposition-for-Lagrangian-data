function k = k_of_x(x)
% Wave vector from x-vector
%
N      = max(size(x));
dx     = x(2)-x(1);
dk     = (2*pi)/(N*dx);
inull  = N/2;
k      = dk*(linspace(1,N,N)-inull);
