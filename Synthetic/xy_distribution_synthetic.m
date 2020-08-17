%Plot the histogram of drifter locations, corresponding to figure 2 in
%paper

clear all
close all
doPlotSetup

load ('Synthetictraj.mat')

Nbins=2000;


ndrfts=size(trajmat_X,2);

ntimesnap_total = size(trajmat_X,1);

converted_X=reshape(trajmat_X,ndrfts*ntimesnap_total,1);
converted_Y=reshape(trajmat_Y,ndrfts*ntimesnap_total,1);


figure
N = hist3([converted_X converted_Y],'Nbins',[Nbins Nbins]);
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(min(converted_X),max(converted_X),size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(min(converted_Y),max(converted_Y),size(N_pcolor,1)); % Rows of N_pcolor
N_plot=log10(N_pcolor);
N_plot(N_pcolor==0)=NaN;


h= pcolor(xl/1000,yl/1000,N_plot);
set(h, 'linestyle', 'none')

colormap((hot.^3)) % Change color scheme 

h.ZData = -max(N_pcolor(:))*ones(size(N_pcolor));
ax = gca;
ax.ZTick(ax.ZTick < 0) = [];
pbaspect([max(converted_X)-min(converted_X) max(converted_Y)-min(converted_Y) 1])
xlabel('x [km]')
ylabel('y [km]')
title('drifter numbers')

Label(1)=1;
Label(end+1)=5;
label=log10(Label);

hC=colorbar; % Display colorbar
set(hC,'Ytick',label,'YTicklabel',Label);
yticks([-80 -40 0 40 80])
xticks([-80 -40 0 40 80])
ax.FontSize = 13; 

savefig('xyhisto_nomap.fig')
