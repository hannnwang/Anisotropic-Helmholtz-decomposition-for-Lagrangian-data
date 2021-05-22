clear all
close all
rstart=0;


%Check dLT
doPlotSetup
Moment=0;
mattitle=sprintf('LASER_all_ensemble_Moment_%d.mat',Moment);
load(mattitle);

rend=radius;
rend_sub=rend/10;

dlt_m0_ensemble=dlt2M_cos_ensemble;
dll_m0_ensemble=dll2M_cos_ensemble;
dtt_m0_ensemble=dtt2M_cos_ensemble;

Moment=2;
mattitle=sprintf('LASER_all_ensemble_Moment_%d.mat',Moment);
load(mattitle);
dll2M_ensemble=sqrt(dll2M_cos_ensemble.^2+dll2M_sin_ensemble.^2);
dtt2M_ensemble=sqrt(dtt2M_cos_ensemble.^2+dtt2M_sin_ensemble.^2);
dlt2M_ensemble=sqrt(dlt2M_cos_ensemble.^2+dlt2M_sin_ensemble.^2);

dll_m2_mag=dll2M_ensemble;
dtt_m2_mag=dtt2M_ensemble;
dlt_m2_mag=dlt2M_ensemble;

iplot_submeso=find(dist_axis<rend_sub & dist_axis>rstart);
iplot_meso=find(dist_axis<=rend);

%Mark the limiting behavior 
load LASER_CuCvCuv.mat
figure(2);subplot(3,1,2) 
SLL0=0.*dist_axis(iplot_meso)+Cu0+Cv0;
SLL2=0.*dist_axis(iplot_meso)+sqrt((Cu0-Cv0)^2+(2*Cuv0)^2);
% SLL0(1:round(length(dist_bin)*0.3))=NaN;
% SLL2(1:round(length(dist_bin)*0.3))=NaN;

plot(dist_axis(iplot_meso),SLL0,'Color','[0.9290 0.6940 0.1250]','Linestyle','-.')
hold on
plot(dist_axis(iplot_meso),SLL2,'Color','[0 0.4470 0.7410]','Linestyle','-.')
hold on
legend('$S_{c0}^{LL}(\infty)$',...
    '$S_{a2}^{LL}(\infty)$',...
    'Location','northoutside',...
    'Interpreter','latex','AutoUpdate','off','Orientation','horizontal') % R2018b and later

figure(2);subplot(3,1,1) 
SLT2=0.*dist_axis(iplot_meso)+sqrt((Cu0-Cv0)^2+(2*Cuv0)^2);

plot(dist_axis(iplot_meso),SLT2,'Color','[0 0.4470 0.7410]','Linestyle','-.')
hold on
legend('$S_{a2}^{LT}(\infty)$',...
    'Location','northoutside',...
    'Interpreter','latex','AutoUpdate','off','Orientation','horizontal') % R2018b and later

figure(2);subplot(3,1,3) 
STT0=0.*dist_axis(iplot_meso)+Cu0+Cv0;
STT2=0.*dist_axis(iplot_meso)+sqrt((Cu0-Cv0)^2+(2*Cuv0)^2);

plot(dist_axis(iplot_meso),STT0,'Color','[0.9290 0.6940 0.1250]','Linestyle','-.')
hold on
plot(dist_axis(iplot_meso),STT2,'Color','[0 0.4470 0.7410]','Linestyle','-.')
hold on
legend('$S_{c0}^{TT}(\infty)$',...
    '$S_{a2}^{TT}(\infty)$',...
    'Location','northoutside',...
    'Interpreter','latex','AutoUpdate','off','Orientation','horizontal') % R2018b and later


% dump=dlt_m0_ensemble;dump(dump>0)=NaN;
dump=-abs(dlt_m0_ensemble);
figure(1);subplot(1,3,2)
loglog(dist_axis(iplot_meso),-dump(iplot_meso),...
    'Color','[0.8500 0.3250 0.0980]','Linestyle','-','Linewidth',3.5)
hold on
loglog(dist_axis(iplot_meso),dll_m0_ensemble(iplot_meso),'Color','[0.9290 0.6940 0.1250]')
hold on
loglog(dist_axis(iplot_meso),dll_m2_mag(iplot_meso),'Color','[0 0.4470 0.7410]')
hold on

figure(1);subplot(1,3,3)
loglog(dist_axis(iplot_meso),-dump(iplot_meso),...
    'Color','[0.8500 0.3250 0.0980]','Linestyle','-','Linewidth',3.5)
hold on
loglog(dist_axis(iplot_meso),dtt_m0_ensemble(iplot_meso),'Color','[0.9290 0.6940 0.1250]')
hold on
loglog(dist_axis(iplot_meso),dtt_m2_mag(iplot_meso),'Color','[0 0.4470 0.7410]')
hold on

figure(1);subplot(1,3,1)
loglog(dist_axis(iplot_meso),-dump(iplot_meso),...
    'Color','[0.8500 0.3250 0.0980]','Linestyle','-','Linewidth',3.5)
hold on
loglog(dist_axis(iplot_meso),dlt_m2_mag(iplot_meso),'Color','[0 0.4470 0.7410]')
hold on



figure(2);subplot(3,1,1) 
plot(dist_axis(iplot_meso),-dump(iplot_meso),...
    'Color','[0.8500 0.3250 0.0980]','Linestyle','-','Linewidth',3.5)
hold on
plot(dist_axis(iplot_meso),dlt_m2_mag(iplot_meso),'Color','[0 0.4470 0.7410]')
hold on



figure(2);subplot(3,1,2) 
plot(dist_axis(iplot_meso),-dump(iplot_meso),...
    'Color','[0.8500 0.3250 0.0980]','Linestyle','-','Linewidth',3.5)
hold on
plot(dist_axis(iplot_meso),dll_m0_ensemble(iplot_meso),'Color','[0.9290 0.6940 0.1250]')
hold on
plot(dist_axis(iplot_meso),dll_m2_mag(iplot_meso),'Color','[0 0.4470 0.7410]')
hold on




figure(2);subplot(3,1,3) 
plot(dist_axis(iplot_meso),-dump(iplot_meso),...
    'Color','[0.8500 0.3250 0.0980]','Linestyle','-','Linewidth',3.5)
hold on
plot(dist_axis(iplot_meso),dtt_m0_ensemble(iplot_meso),'Color','[0.9290 0.6940 0.1250]')
hold on
plot(dist_axis(iplot_meso),dtt_m2_mag(iplot_meso),'Color','[0 0.4470 0.7410]')
hold on


Moment=4;
mattitle=sprintf('LASER_all_ensemble_Moment_%d.mat',Moment);
load(mattitle);

dll2M_ensemble=sqrt(dll2M_cos_ensemble.^2+dll2M_sin_ensemble.^2);
dtt2M_ensemble=sqrt(dtt2M_cos_ensemble.^2+dtt2M_sin_ensemble.^2);
dlt2M_ensemble=sqrt(dlt2M_cos_ensemble.^2+dlt2M_sin_ensemble.^2);

dll_2M_mag=dll2M_ensemble;
dtt_2M_mag=dtt2M_ensemble;
dlt_2M_mag=dlt2M_ensemble;


figure(2);subplot(3,1,1) 
plot(dist_axis(iplot_meso),dlt_2M_mag(iplot_meso),'Color','[0.4660 0.6740 0.1880]')
hold on

figure(2);subplot(3,1,2) 
plot(dist_axis(iplot_meso),dll_2M_mag(iplot_meso),'Color','[0.4660 0.6740 0.1880]')
hold on

figure(2);subplot(3,1,3) 
plot(dist_axis(iplot_meso),dtt_2M_mag(iplot_meso),'Color','[0.4660 0.6740 0.1880]')
hold on


Moment=100;
mattitle=sprintf('LASER_all_ensemble_Moment_%d.mat',Moment);
load(mattitle);

dll2M_ensemble=sqrt(dll2M_cos_ensemble.^2+dll2M_sin_ensemble.^2);
dtt2M_ensemble=sqrt(dtt2M_cos_ensemble.^2+dtt2M_sin_ensemble.^2);
dlt2M_ensemble=sqrt(dlt2M_cos_ensemble.^2+dlt2M_sin_ensemble.^2);

dll_2M_mag=dll2M_ensemble;
dtt_2M_mag=dtt2M_ensemble;
dlt_2M_mag=dlt2M_ensemble;

figure(1);subplot(1,3,1)
loglog(dist_axis(iplot_meso),dlt_2M_mag(iplot_meso),'Color','[0.4940 0.1840 0.5560]')
hold on

legend('$|S_{c0}^{LT}|$',...
    '$S_{a2}^{LT}$',...
    '$S_{a100}^{LT}$',...
        'Location','northoutside',...
    'Interpreter','latex','AutoUpdate','off',...
    'Location','southwest','FontSize',12,'Orientation','horizontal') % R2018b and later


xLimits = get(gca,'XLim');  %# Get the range of the x axis
yLimits = get(gca,'YLim');
% AREA = area([4000 xLimits(2)], [yLimits(2) yLimits(2)],yLimits(1),'FaceAlpha',.3);
% set(AREA,'EdgeColor', 'none'); % make area's edge invisible


xlabel('m')
ylabel('m^2/s^2')
pbaspect([1.2 1 1])

ax = gca;
ax.FontSize = 13; 

figure(1);subplot(1,3,3)


legend('$|S_{c0}^{LT}|$',...
    '$S_{c0}^{TT}$',...
    '$S_{a2}^{TT}$',...
        'Location','northoutside',...
    'Interpreter','latex','AutoUpdate','off',...
    'Location','southwest','FontSize',12,'Orientation','horizontal') % R2018b and later


xLimits = get(gca,'XLim');  %# Get the range of the x axis
yLimits = get(gca,'YLim');
% AREA = area([4000 xLimits(2)], [yLimits(2) yLimits(2)],yLimits(1),'FaceAlpha',.3);
% set(AREA,'EdgeColor', 'none'); % make area's edge invisible


xlabel('m')
ylabel('m^2/s^2')
pbaspect([1.2 1 1])

ax = gca;
ax.FontSize = 13; 

figure(1);subplot(1,3,2)


legend('$|S_{c0}^{LT}|$',...
    '$S_{c0}^{LL}$',...
    '$S_{a2}^{LL}$',...
        'Location','northoutside',...
    'Interpreter','latex','AutoUpdate','off',...
    'Location','southwest','FontSize',12,'Orientation','horizontal') % R2018b and later


xLimits = get(gca,'XLim');  %# Get the range of the x axis
yLimits = get(gca,'YLim');
% AREA = area([4000 xLimits(2)], [yLimits(2) yLimits(2)],yLimits(1),'FaceAlpha',.3);
% set(AREA,'EdgeColor', 'none'); % make area's edge invisible


xlabel('m')
ylabel('m^2/s^2')
pbaspect([1.2 1 1])

ax = gca;
ax.FontSize = 13; 


figure(2);subplot(3,1,1) 


xlabel('m')
ylabel('m^2/s^2')
pbaspect([1.6 1 1])

ax = gca;
ax.FontSize = 13; 


figure(2);subplot(3,1,2) 


xlabel('m')
ylabel('m^2/s^2')
pbaspect([1.6 1 1])

ax = gca;
ax.FontSize = 13; 




figure(2);subplot(3,1,3) 

xlabel('m')
ylabel('m^2/s^2')
pbaspect([1.6 1 1])

ax = gca;
ax.FontSize = 13; 
load('LASER_all_alpha.mat');

angleva=0.*dist_axis;
for i=1:length(dist_axis)
    mm=mean(nanglecount(i,:));
    anglevar(i)=mean((nanglecount(i,:)-mm).^2);
end
alphafactor_var=anglevar/mean(mm)^2;


fg2=figure(2);
% subplot(3,1,1)
% markalert
% subplot(3,1,2)
% markalert
% subplot(3,1,3)
% markalert
sameaxes([],[fg2])
savefig('figure2_3panel.fig')


