function [h] = hwdrawcircle(radius,ox,oy)
hold on
th = 0:pi/50:2*pi;
xunit = radius * cos(th)+ox;
yunit = radius * sin(th)+oy;
hold on 
plot(xunit,yunit,'k--')
end
