hold on
[irneg,ithetaneg]=find(Ddd_2D<0);
Xneg=raxis(irneg).*cos(thetaxis(ithetaneg))/1000;
Yneg=raxis(irneg).*sin(thetaxis(ithetaneg))/1000;
negb=boundary(Xneg,Yneg);
plot(Xneg(negb),Yneg(negb))