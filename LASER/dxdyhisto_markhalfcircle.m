clear all
close all


open dxdy_histo_N1000.fig
% circlespec
figure(1)
hold all
radius=5*10^5;
hwdrawhalfcircle(radius,0,0)

radius=1*10^5;
hwdrawhalfcircle(radius,0,0)

radius=3*10^5;
hwdrawhalfcircle(radius,0,0)


savefig('dxdy_histo_Rmarked.fig')