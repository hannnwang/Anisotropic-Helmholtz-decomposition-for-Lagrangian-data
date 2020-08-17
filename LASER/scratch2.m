
N = hist3([deltaX deltaY],'Nbins',[Nbins Nbins]);
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(min(deltaX),max(deltaX),size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(min(deltaY),max(deltaY),size(N_pcolor,1)); % Rows of N_pcolor

h = pcolor(xl,yl,N_pcolor);
set(h, 'linestyle', 'none')
colormap(flipud(bone)) % Change color scheme 
colorbar % Display colorbar
h.ZData = -max(N_pcolor(:))*ones(size(N_pcolor));
ax = gca;
ax.ZTick(ax.ZTick < 0) = [];
axis square