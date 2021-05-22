clear all
close all
doPlotSetup

load ('LASER_all.mat')

Nbins=500;
Rearth=6.3782*10^6;


ndrfts=size(trajmat_X,2);

ntimesnap_total = size(trajmat_X,1);

converted_X=reshape(trajmat_X.*cosd(trajmat_Y)*Rearth*pi/180,ndrfts*ntimesnap_total,1);

converted_Y=reshape(trajmat_Y*Rearth*pi/180,ndrfts*ntimesnap_total,1);

% %Cutting some parts of X and Y
% xupl=-8.49e+6;
% xlwl=-8.71e+6;
% yupl=3.265e+6;
% ylwl=3.135e+6;
% 
%  
% cutind=find(converted_Y>ylwl & converted_Y<yupl...
%      & converted_X>xlwl & converted_X<xupl);
% converted_X=converted_X(cutind);
% converted_Y=converted_Y(cutind);

figure

N = hist3([converted_X converted_Y],'Nbins',[Nbins Nbins]);
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(min(converted_X),max(converted_X),size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(min(converted_Y),max(converted_Y),size(N_pcolor,1)); % Rows of N_pcolor
N_plot=log(N_pcolor);
N_plot(N_pcolor==0)=NaN;
h = pcolor(xl,yl,N_plot);
set(h, 'linestyle', 'none')
colormap((jet)) % Change color scheme 
colorbar % Display colorbar
h.ZData = -max(N_pcolor(:))*ones(size(N_pcolor));
ax = gca;
ax.ZTick(ax.ZTick < 0) = [];
pbaspect([max(converted_X)-min(converted_X) max(converted_Y)-min(converted_Y) 1])
xlabel('x [m]')
ylabel('y [m]')
title('Log[drifter numbers]')
