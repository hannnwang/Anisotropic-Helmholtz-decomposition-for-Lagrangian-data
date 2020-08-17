%Make the histogram for separation vectors (dx,dy). Corresponding to figure
%3 in paper.
clear all
close all
doPlotSetup

load ('Synthetictraj.mat')
Nbins=1000;
% calculate time series of pair separation. Same as in dist_axis_synthetic.m
ninterval=1;

nsnapshot_total = size(trajmat_X,1);

nsnapshot_end=nsnapshot_total;

timevec=1:ninterval:nsnapshot_end;

nsnapshot_sampled = length(timevec);

n=1;

tic

for i=1:ndrfts-1
    for j = i+1:ndrfts           
        for it=1:ninterval:nsnapshot_end
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
toc
disp('X1,X2,etc. are formed.')
dX=(X1-X2);
dY=(Y1-Y2);

clearvars X1 X2 Y1 Y2
dX(dY<0)=-dX(dY<0);
dY(dY<0)=-dY(dY<0);

[sizex, sizey]=size(dX);
dX=reshape(dX,sizex*sizey,1);
dY=reshape(dY,sizex*sizey,1);


figure
N = hist3([dX/1000 dY/1000],'Nbins',[Nbins Nbins]);
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(min(dX),max(dX),size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(min(dY),max(dY),size(N_pcolor,1)); % Rows of N_pcolor
N_plot=log10(N_pcolor);
N_plot(N_pcolor==0)=NaN;
h= pcolor(xl/1000,yl/1000,N_plot);
set(h, 'linestyle', 'none')
colormap((bone.^0.8)) % Change color scheme 
hC=colorbar; % Display colorbar
h.ZData = -max(N_pcolor(:))*ones(size(N_pcolor));
ax = gca;
ax.ZTick(ax.ZTick < 0) = [];
pbaspect([max(dX)-min(dX) max(dY)-min(dY) 1])

Lmax=floor(nanmax(nanmax(N_plot)));
Label=10.^(0:1:Lmax);
label=log10(Label);
set(hC,'Ytick',label,'YTicklabel',Label);
title('pair numbers')
xlabel('\Delta x [km]')
ylabel('\Delta y [km]')

xticks([-150 -100 -50 0 50 100 150])
yticks([0 30 60 90])
ax.FontSize = 13; 

savefig('dxdy_histo.fig')

