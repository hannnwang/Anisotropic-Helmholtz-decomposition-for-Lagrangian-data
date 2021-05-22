clear all
close all

doPlotSetup

open dxdy_histo.fig
circlespec
figure(1)
hold all
hwdrawhalfcircle(radius,0,0)

savefig('dxdy_histo_Rmarked.fig')