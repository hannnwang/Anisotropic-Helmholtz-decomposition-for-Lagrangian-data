function [C_0] = hwC0(C,xr,yr,k,l,K,L,Kappa,N,dk,dl)
%Calculate the c0 mode of a 2D function
%Notation used in this code:
%C=C_0+(C_2sin)*sin(2*(alpha-theta))+(C_4sin)*sin(4*(alpha-theta))+.......
%C is input: real, even, and N-by-N; 
%C_0 starts from r=0

Ct=real(hwfft2(xr,yr,k,l,C));%Note that the real part is taken since C is assumed real and even.

for i=N/2:N
    ri=xr(i);
    J0kr=besselj(0,Kappa*ri);    
    C_0(i-N/2+1)=sum(sum(Ct.*J0kr))*dk*dl/(2*pi)^2;
end
end

