clear all
close all
doPlotSetup

cd 'D:\OBresearch\Project_aniso\LASER\20200214\all\trapezoidal'
  Moment=0;
mattitle=sprintf('LASER_all_ensemble_Moment_%d.mat',Moment);
load(mattitle);

dist_axis=dist_axis';

xcutoff=3*10^5;
iplot=find(dist_axis<xcutoff);

rv=dist_axis(iplot);
raxis=rv;
%raxis(2:end+1)=rv;raxis(1)=0;
thetaxis=0:pi/180:pi;
[R,Theta]=ndgrid(raxis,thetaxis);

Drr_2D=0.*R;
Ddd_2D=0.*R;

Drrc0_2D=0.*R;
Dddc0_2D=0.*R;

%moment
intgrd2M_cos_ensemble=(-Moment.*dlt2M_sin_ensemble-dll2M_cos_ensemble+...
    dtt2M_cos_ensemble)./dist_axis;
intgrd2M_sin_ensemble=(+Moment.*dlt2M_cos_ensemble-dll2M_sin_ensemble+...
    dtt2M_sin_ensemble)./dist_axis;

Drr2M_cos_ensemble=dtt2M_cos_ensemble+...
    intRHS_trapz_nobin(dist_axis,intgrd2M_cos_ensemble);

Ddd2M_cos_ensemble=dll2M_cos_ensemble-...
    intRHS_trapz_nobin(dist_axis,intgrd2M_cos_ensemble);

Drr2M_sin_ensemble=dtt2M_sin_ensemble+...
    intRHS_trapz_nobin(dist_axis,intgrd2M_sin_ensemble);

Ddd2M_sin_ensemble=dll2M_sin_ensemble-...
    intRHS_trapz_nobin(dist_axis,intgrd2M_sin_ensemble);

Drr_c0_ensemble=Drr2M_cos_ensemble;
Ddd_c0_ensemble=Ddd2M_cos_ensemble;

%M=2
 Moment=2;
mattitle=sprintf('LASER_all_ensemble_Moment_%d.mat',Moment);
load(mattitle);
dist_axis=dist_axis';

%moment
intgrd2M_cos_ensemble=(-Moment.*dlt2M_sin_ensemble-dll2M_cos_ensemble+...
    dtt2M_cos_ensemble)./dist_axis;
intgrd2M_sin_ensemble=(+Moment.*dlt2M_cos_ensemble-dll2M_sin_ensemble+...
    dtt2M_sin_ensemble)./dist_axis;

Drr2M_cos_ensemble=dtt2M_cos_ensemble+...
    intRHS_trapz_nobin(dist_axis,intgrd2M_cos_ensemble);

Ddd2M_cos_ensemble=dll2M_cos_ensemble-...
    intRHS_trapz_nobin(dist_axis,intgrd2M_cos_ensemble);

Drr2M_sin_ensemble=dtt2M_sin_ensemble+...
    intRHS_trapz_nobin(dist_axis,intgrd2M_sin_ensemble);

Ddd2M_sin_ensemble=dll2M_sin_ensemble-...
    intRHS_trapz_nobin(dist_axis,intgrd2M_sin_ensemble);

Drr_c2_ensemble=Drr2M_cos_ensemble;
Ddd_c2_ensemble=Ddd2M_cos_ensemble;
Drr_s2_ensemble=Drr2M_sin_ensemble;
Ddd_s2_ensemble=Ddd2M_sin_ensemble;

%M=4
 Moment=4;
mattitle=sprintf('LASER_all_ensemble_Moment_%d.mat',Moment);
load(mattitle);
dist_axis=dist_axis';

%moment
intgrd2M_cos_ensemble=(-Moment.*dlt2M_sin_ensemble-dll2M_cos_ensemble+...
    dtt2M_cos_ensemble)./dist_axis;
intgrd2M_sin_ensemble=(+Moment.*dlt2M_cos_ensemble-dll2M_sin_ensemble+...
    dtt2M_sin_ensemble)./dist_axis;

Drr2M_cos_ensemble=dtt2M_cos_ensemble+...
    intRHS_trapz_nobin(dist_axis,intgrd2M_cos_ensemble);

Ddd2M_cos_ensemble=dll2M_cos_ensemble-...
    intRHS_trapz_nobin(dist_axis,intgrd2M_cos_ensemble);

Drr2M_sin_ensemble=dtt2M_sin_ensemble+...
    intRHS_trapz_nobin(dist_axis,intgrd2M_sin_ensemble);

Ddd2M_sin_ensemble=dll2M_sin_ensemble-...
    intRHS_trapz_nobin(dist_axis,intgrd2M_sin_ensemble);

Drr_c4_ensemble=Drr2M_cos_ensemble;
Ddd_c4_ensemble=Ddd2M_cos_ensemble;
Drr_s4_ensemble=Drr2M_sin_ensemble;
Ddd_s4_ensemble=Ddd2M_sin_ensemble;

%Reconstruct the 2D functions

for itheta=1:length(thetaxis)
        cos2theta=cos(2*thetaxis(itheta));
        sin2theta=sin(2*thetaxis(itheta));
        cos4theta=cos(4*thetaxis(itheta));
        sin4theta=sin(4*thetaxis(itheta));
   for ir=1:length(raxis)
        Drr_2D(ir,itheta)=Drr_c0_ensemble(ir)+...
            Drr_c2_ensemble(ir)*cos2theta+...
            Drr_s2_ensemble(ir)*sin2theta+...
            Drr_c4_ensemble(ir)*cos4theta+...
            Drr_s4_ensemble(ir)*sin4theta;
        
        Drrc0_2D(ir,itheta)=Drr_c0_ensemble(ir);
        
          Ddd_2D(ir,itheta)=Ddd_c0_ensemble(ir)+...
            Ddd_c2_ensemble(ir)*cos2theta+...
            Ddd_s2_ensemble(ir)*sin2theta+...
            Ddd_c4_ensemble(ir)*cos4theta+...
            Ddd_s4_ensemble(ir)*sin4theta;
        
                Dddc0_2D(ir,itheta)=Ddd_c0_ensemble(ir);

   end
end

Drr_2D_plot=Drr_2D;
%Drr_2D_plot(Drr_2D<0)=NaN;
Ddd_2D_plot=Ddd_2D;
% Ddd_2D_plot(Dddc0_2D<0)=NaN;
ratio_2D=abs(Drr_2D_plot)./(abs(Drr_2D_plot)+abs(Ddd_2D_plot));
% ratio_2D(Ddd_2D<0 & Drr_2D./abs(Ddd_2D)>=0)=NaN;
% ratio_2D(Drr_2D<0 & Ddd_2D./abs(Drr_2D)>=0)=NaN;

% 
% %Plot Drr
% figure(1)
% subplot(1,2,1)
% X=R.*cos(Theta)/1000;
% Y=R.*sin(Theta)/1000;
% hold on;
% contourf(X,Y,Drr_2D,length(dist_axis)/5,'edgecolor','none')
% ylim([0 (raxis(end)/1000)])
% xlim([-(raxis(end)/1000) (raxis(end)/1000)])
% pbaspect([2 1 1]);
% colorbar
% ax1 = gca;                   % gca = get current axis
% ax1.YAxis.Visible = 'off';   % remove y-axis
% xticks([-300 -250 -200 -150 -100 -50 0 50 100 150 200 250 300])
% set(gca,'TickDir','out');
% xlabel('r [km]')
% subplot(1,2,2)
% base=5;%log base
% Xlog=log(R)/log(base).*cos(Theta);
% Ylog=log(R)/log(base).*sin(Theta);
% hold on;
% contourf(Xlog,Ylog,log(Drr_2D)/log(base),length(dist_axis)/5,'edgecolor','none')
% ylim([0 abs(log(raxis(end))/log(base))])
% xlim([-log(raxis(end))/log(base) log(raxis(end))/log(base)])
% pbaspect([2 1 1]);
% ax1 = gca;                   % gca = get current axis
% ax1.YAxis.Visible = 'off';   % remove y-axis
% 
% xticks([-log(10^5)/log(base) -log(10^4)/log(base) -log(10^3)/log(base)...
%     -log(10^2)/log(base) -log(10)/log(base)...
%     0 log(10)/log(base) log(100)/log(base) log(1000)/log(base) log(10^4)/log(base)... ...
%      log(10^5)/log(base)])
% xticklabels({-100,-10,-1,-0.1,-0.01,0,0.01,0.1,1,10,100})
% set(gca,'TickDir','out');
% xlabel('r [km]')
% 
% colormap((hot)) % Change color scheme 
% hC=colorbar; % Display colorbar
% 
% L(1)=10^(-5);
% L(end+1)=L(end)*10;
% L(end+1)=L(end)*10;
% L(end+1)=L(end)*10;
% L(end+1)=L(end)*10;
% 
% l=log(L)/log(base);
% set(hC,'Ytick',l,'YTicklabel',L);
% 

figure(2)
%Plot ratio
base=5;%log base
Xlog=log(R)/log(base).*cos(Theta);
Ylog=log(R)/log(base).*sin(Theta);
hold on;
contourf(Xlog,Ylog,ratio_2D,10,'edgecolor','none')
ylim([0 abs(log(raxis(end))/log(base))])
xlim([-log(raxis(end))/log(base) log(raxis(end))/log(base)])
pbaspect([2 1 1]);
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis

xticks([-log(10^5)/log(base) -log(10^4)/log(base) -log(10^3)/log(base)...
    -log(10^2)/log(base) -log(10)/log(base)...
    0 log(10)/log(base) log(100)/log(base) log(1000)/log(base) log(10^4)/log(base)... ...
     log(10^5)/log(base)])
xticklabels({-100,-10,-1,-0.1,-0.01,0,0.01,0.1,1,10,100})
set(gca,'TickDir','out');
xlabel('r [km]')

colormap(copper);
hC=colorbar; % Display colorbar
caxis([0 1])


