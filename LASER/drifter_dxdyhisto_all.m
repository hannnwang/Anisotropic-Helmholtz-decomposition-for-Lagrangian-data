%20200104(needs to be run on HPC)
%Make the histogram for (dx,dy) (all LASER data)
clear all
close all
doPlotSetup

load ('LASER_all.mat')
Nbins=1000;
Rearth=6.3782*10^6;
% calculate time series of pair separation
ninterval=1

ndrfts=size(trajmat_X,2)

% ndrfts_sampled=ndrfts%Include all floaters

ntimesnap_total = size(trajmat_X,1);

ntimesnap_end=ntimesnap_total/3;

timevec=1:ninterval:ntimesnap_end;

n=1;
for i=1:ndrfts-1
    for j = i+1:ndrfts     
            %(X1, Y1 are in non-standard units)
    
        for it=1:ninterval:ntimesnap_end
            if ~isnan(trajmat_X(it,i)) && ~isnan(trajmat_X(it,j))
                X1(n) = trajmat_X(it,i);
                X2(n) = trajmat_X(it,j);
                Y1(n) = trajmat_Y(it,i);
                Y2(n) = trajmat_Y(it,j);
                n = n+1;
          end
        end
    end
end

disp('X1,X2,etc. are formed.')
dX=(X1.*cosd(Y1)-X2.*cosd(Y2))*Rearth*pi/180;
dY=(Y1-Y2)*Rearth*pi/180;

clearvars X1 X2 Y1 Y2
dX(dY<0)=-dX(dY<0);
dY(dY<0)=-dY(dY<0);

[sizex, sizey]=size(dX);
dX=reshape(dX,sizex*sizey,1);
dY=reshape(dY,sizex*sizey,1);

figure
N = hist3([dX dY],'Nbins',[Nbins Nbins]);
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(min(dX),max(dX),size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(min(dY),max(dY),size(N_pcolor,1)); % Rows of N_pcolor
N_plot=log(N_pcolor);
N_plot(N_pcolor==0)=NaN;
h = pcolor(xl,yl,N_plot);
set(h, 'linestyle', 'none')
colormap((jet)) % Change color scheme 
colorbar % Display colorbar
h.ZData = -max(N_pcolor(:))*ones(size(N_pcolor));
ax = gca;
ax.ZTick(ax.ZTick < 0) = [];
pbaspect([max(dX)-min(dX) max(dY)-min(dY) 1])
xlabel('\Delta x [m]')
ylabel('\Delta y [m]')
title('Log[drifter numbers]')
savefig('dxdy_histo.fig')

