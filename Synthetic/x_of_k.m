function x = x_of_k(k)
% x-vector with leading zero from k-vector 
%
N      = max(size(k));
dk     = k(2)-k(1);
dx     = 2*pi/(N*dk);
x      = dx*(linspace(1,N,N)-1);
