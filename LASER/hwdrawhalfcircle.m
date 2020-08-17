function [h] = hwdrawhalfcircle(radius,ox,oy)
th = 0:pi/50:pi;
xunit = radius * cos(th)+ox;
yunit = radius * sin(th)+oy;
hold on 
plot(xunit,yunit,'r--')
end

