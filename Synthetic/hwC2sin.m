function [C_2sin] = hwC2sin(C,theta,xr,yr,k,l,K,L,Kappa,N,dk,dl)
%Calculate the s2 mode of a 2D function
%Notation used in this code:
%C=C_0+(C_2sin)*sin(2*(alpha-theta))+(C_4sin)*sin(4*(alpha-theta))+.......
%C is input: real, even, and N-by-N; 
%Output starts from r=0

Ct=real(hwfft2(xr,yr,k,l,C));%Note that the real part is taken since C is assumed real and even.

atanlk=atan2(K,L);

coe_m2_sin=-sin(2*(pi/2-atanlk-theta));
for i=(N/2+1):N
    ri=xr(i);
    %J0kr=besselj(0,Kappa*ri);

    J2kr=besselj(2,Kappa*ri);
    
 C_2sin(i-N/2+1)=sum(sum(Ct.*J2kr.*coe_m2_sin))*dk*dl/(2*pi*pi);

end
C_2sin(1)=0;
end