function [C6cos] = hwC6cos(C,theta,xr,yr,k,l,K,L,Kappa,N,dk,dl)
%Calculate the c6 mode of a 2D function
%Notation used in this code:
%C=C_0+(C_2sin)*sin(2*(alpha-theta))+(C_4sin)*sin(4*(alpha-theta))+.......
%C is input: real, even, and N-by-N; 
%Output starts from r=0

Ct=real(hwfft2(xr,yr,k,l,C));%Note that the real part is taken since C is assumed real and even.

atanlk=atan2(K,L);

coe_m6_cos=cos(6*(pi/2-atanlk-theta));

for i=(N/2+1):N%starting from the 2nd element in C4cos
    ri=xr(i);
    
    kr=Kappa*ri;
    
    J0kr=besselj(0,kr);
    
    J1kr=besselj(1,kr);
    
 intgd_6=1./(kr.^5).*(kr.*(1920-144*...
        kr.^2+kr.^4).*J0kr-6*(640-128.*kr.^2+3*kr.^4).*...
        J1kr);
    intgd_6(N/2,N/2)=0;%know this analytically
   
    C6cos(i-N/2+1)=sum(sum(Ct.*intgd_6.*...
        coe_m6_cos))*dk*dl/(2*pi*pi);
end
C6cos(1)=0;
end