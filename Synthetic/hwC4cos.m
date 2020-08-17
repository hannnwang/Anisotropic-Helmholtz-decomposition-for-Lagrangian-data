function [C4cos] = hwC4cos(C,theta,xr,yr,k,l,K,L,Kappa,N,dk,dl)
%Calculate the c4 mode of a 2D function
%Notation used in this code:
%C=C_0+(C_2sin)*sin(2*(alpha-theta))+(C_4sin)*sin(4*(alpha-theta))+.......
%C is input: real, even, and N-by-N; 
%Output starts from r=0

Ct=real(hwfft2(xr,yr,k,l,C));%Note that the real part is taken since C is assumed real and even.

atanlk=atan2(K,L);

coe_m4_cos=cos(4*(pi/2-atanlk-theta));

for i=(N/2+1):N
    ri=xr(i);
    
    kr=Kappa*ri;
    
    J0kr=besselj(0,kr);
    
    J1kr=besselj(1,kr);
    
    
    intgd_4cos=1./kr.^3.*(kr.*...
        (-24+kr.^2).*J0kr-8*(-6+kr.^2).*J1kr);
    intgd_4cos(N/2,N/2)=0;%know this analytically
    
    
    C4cos(i-N/2+1)=sum(sum(Ct.*intgd_4cos.*...
        coe_m4_cos))*dk*dl/(2*pi*pi);
  
end
C4cos(1)=0;
end