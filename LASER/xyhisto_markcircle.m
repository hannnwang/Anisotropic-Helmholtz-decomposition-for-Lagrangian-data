clear all
close all
doPlotSetup
set(0,'DefaultLineLineWidth',1)

open xyhisto.fig
circlespec
figure(1)
hold all
hwdrawcircle(radius,center_x,center_y)

savefig('xyhisto_circlemarked.fig')