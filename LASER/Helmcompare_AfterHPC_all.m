clear all
close all
doPlotSetup

cd 'D:\OBresearch\Project_aniso\LASER\20200214\all\trapezoidal'
  Moment=0;
mattitle=sprintf('LASER_all_ensemble_Moment_%d.mat',Moment);
load(mattitle);

fg2=figure(2);

    
dist_axis=dist_axis';

xcutoff=3*10^5;
iplot=find(dist_axis<xcutoff);


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

if Moment>0
Drr2M_mag=sqrt(Drr2M_sin_ensemble.^2+Drr2M_cos_ensemble.^2);
Ddd2M_mag=sqrt(Ddd2M_sin_ensemble.^2+Ddd2M_cos_ensemble.^2);
else 
Drr2M_mag=Drr2M_sin_ensemble;
Ddd2M_mag=Ddd2M_sin_ensemble;
end

Drr_m0_ensemble=Drr2M_mag;
Ddd_m0_ensemble=Ddd2M_mag;

figure(2)
subplot(1,3,1) 
loglog(dist_axis(iplot)/1000,Drr_m0_ensemble(iplot),'Color','[0.0313, 0.3164, 0.6094 ]');
hold on
loglog(dist_axis(iplot)/1000,Ddd_m0_ensemble(iplot),'Color','[0.4157, 0.3176, 0.6392 ]');

pbaspect([1.2 1 1]); axis tight  

%title('0th')
%xlabel('r [km]')
lgd.FontSize = 16;

Momentvec_plot(1)=2;
Momentvec_plot(2)=4;

for imoment=1:length(Momentvec_plot)
    Moment=Momentvec_plot(imoment);
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

if Moment>0
Drr2M_mag=sqrt(Drr2M_sin_ensemble.^2+Drr2M_cos_ensemble.^2);
Ddd2M_mag=sqrt(Ddd2M_sin_ensemble.^2+Ddd2M_cos_ensemble.^2);
else 
Drr2M_mag=Drr2M_sin_ensemble;
Ddd2M_mag=Ddd2M_sin_ensemble;
end

figure(2)
subplot(1,3,imoment+1)
loglog(dist_axis(iplot)/1000,Drr2M_mag(iplot),'Color','[0.4980, 0.8039, 0.7333 ]')
hold on
loglog(dist_axis(iplot)/1000,Ddd2M_mag(iplot),'Color','[0.9804, 0.6235, 0.7098 ]')

fm021=loglog(dist_axis(iplot)/1000,Drr_m0_ensemble(iplot),'Color','[0.0313, 0.3164, 0.6094 ]',...
    'Linestyle','-','Linewidth',1);
hold on
fm022=loglog(dist_axis(iplot)/1000,Ddd_m0_ensemble(iplot),'Color','[0.4157, 0.3176, 0.6392 ]',...
    'Linestyle','-','Linewidth',1);
% fm021.Color(4) = 0.5;
% fm022.Color(4) = 0.5;

pbaspect([1.2 1 1]); axis tight  
if Moment==2
    figtitle=sprintf('2nd mode',Moment);
else
figtitle=sprintf('%d th mode',Moment);
end
%title(figtitle)
%xlabel('r [km]')
lgd.FontSize = 16;

end

cd 'D:\OBresearch\Project_aniso\LASER\20200214\all\noweighting'
Moment=0;
mattitle=sprintf('LASER_all_ensemble_Moment_%d.mat',Moment);
load(mattitle);

dist_axis=dist_axis';

xcutoff=3*10^5;
iplot=find(dist_axis<xcutoff);


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

if Moment>0
Drr2M_mag=sqrt(Drr2M_sin_ensemble.^2+Drr2M_cos_ensemble.^2);
Ddd2M_mag=sqrt(Ddd2M_sin_ensemble.^2+Ddd2M_cos_ensemble.^2);
else 
Drr2M_mag=Drr2M_sin_ensemble;
Ddd2M_mag=Ddd2M_sin_ensemble;
end

Drr_m0_ensemble=Drr2M_mag;
Ddd_m0_ensemble=Ddd2M_mag;

figure(2)
subplot(1,3,1) 
loglog(dist_axis(iplot)/1000,Drr_m0_ensemble(iplot),'Color','[0.6196, 0.7922, 0.8824 ]',...
    'Linestyle','-');
hold on
loglog(dist_axis(iplot)/1000,Ddd_m0_ensemble(iplot),'Color','[0.7373, 0.7412, 0.8627 ]',...
    'Linestyle','-');


%plot negative Ddd
subplot(1,3,1) 
negplot_trapezoidal(1)
subplot(1,3,2) 
negplot_trapezoidal(1)

subplot(1,3,3) 
negplot_trapezoidal(1)



subplot(1,3,1)
ylabel('[m^2/s^2]')
xlabel('r [km]')
axis tight
ax = gca;
ax.FontSize = 13; 
legend('$D^{c0}_{RR}$','$D^{c0}_{DD}$',...
    '$D^{c0}_{RR}$, unw.','$D^{c0}_{DD}$, unw.',...
    '$D_{DD}^{c0}$, neg.',...
    'Interpreter','latex','Location','northwest',...
    'Autoupdate','Off',...
    'Fontsize',16)

subplot(1,3,2)
legend('$D^{a2}_{RR}$','$D^{a2}_{DD}$',...
    '$D^{c0}_{RR}$','$D^{c0}_{DD}$',...
    '$D_{DD}^{c0}$, neg.',...
    'Interpreter','latex','Location','northwest','Autoupdate','Off',...
    'Fontsize',16)
xlabel('r [km]')
axis tight
ax = gca;
ax.FontSize = 13; 
set(gca,'YTickLabel',[]);
subplot(1,3,3)
xlabel('r [km]')
axis tight
ax = gca;
ax.FontSize = 13; 
legend('$D^{a4}_{RR}$','$D^{a4}_{DD}$',...
    '$D^{c0}_{RR}$','$D^{c0}_{DD}$',...
     '$D_{DD}^{c0}$, neg.',...
    'Interpreter','latex','Location','northwest','Autoupdate','Off',...
    'Fontsize',16)
set(gca,'YTickLabel',[]);

f2=figure(2);
sameaxes([],[f2]);

cd 'D:\OBresearch\Project_aniso\LASER\20200214'

savefig('Helm.fig')

function a=negplot_noweighting()
cd 'D:\OBresearch\Project_aniso\LASER\20200214\all\noweighting'
Moment=0;
mattitle=sprintf('LASER_all_ensemble_Moment_%d.mat',Moment);
load(mattitle);

dist_axis=dist_axis';

xcutoff=3*10^5;
iplot=find(dist_axis<xcutoff);


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

if Moment>0
Drr2M_mag=sqrt(Drr2M_sin_ensemble.^2+Drr2M_cos_ensemble.^2);
Ddd2M_mag=sqrt(Ddd2M_sin_ensemble.^2+Ddd2M_cos_ensemble.^2);
else 
Drr2M_mag=Drr2M_sin_ensemble;
Ddd2M_mag=Ddd2M_sin_ensemble;
end

Drr_m0_ensemble=Drr2M_mag;
Ddd_m0_ensemble=Ddd2M_mag;

Ddd_m0_ensemble_neg=-Ddd_m0_ensemble;

% loglog(dist_axis(iplot)/1000,Ddd_m0_ensemble_neg(iplot),'Color','[0.4157, 0.3176, 0.6392 ]',...
%     'Linestyle','--');
% hold on
% a=loglog(dist_axis(iplot)/1000,Ddd_m0_ensemble_neg(iplot),'Color','[0.4157, 0.3176, 0.6392 ]',...
%     'Linewidth',6);
% a.Color(4)=0.3;
loglog(dist_axis(iplot)/1000,Ddd_m0_ensemble_neg(iplot),'Color','[0.5, 0, 0.1484 ]');


end
function a=negplot_trapezoidal(flag)
cd 'D:\OBresearch\Project_aniso\LASER\20200214\all\trapezoidal'
Moment=0;
mattitle=sprintf('LASER_all_ensemble_Moment_%d.mat',Moment);
load(mattitle);

dist_axis=dist_axis';

xcutoff=3*10^5;
iplot=find(dist_axis<xcutoff);


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

if Moment>0
Drr2M_mag=sqrt(Drr2M_sin_ensemble.^2+Drr2M_cos_ensemble.^2);
Ddd2M_mag=sqrt(Ddd2M_sin_ensemble.^2+Ddd2M_cos_ensemble.^2);
else 
Drr2M_mag=Drr2M_sin_ensemble;
Ddd2M_mag=Ddd2M_sin_ensemble;
end

Drr_m0_ensemble=Drr2M_mag;
Ddd_m0_ensemble=Ddd2M_mag;

Ddd_m0_ensemble_neg=-Ddd_m0_ensemble;
if(flag==0)
% loglog(dist_axis(iplot)/1000,Ddd_m0_ensemble_neg(iplot),'Color','[0.4157, 0.3176, 0.6392 ]');
% 
% hold on
% a=loglog(dist_axis(iplot)/1000,Ddd_m0_ensemble_neg(iplot),'Color','[0.4157, 0.3176, 0.6392 ]',...
%     'Linewidth',6);
% a.Color(4)=0.3;
loglog(dist_axis(iplot)/1000,Ddd_m0_ensemble_neg(iplot),'Color','[0.5, 0, 0.1484 ]','LineStyle','--');

elseif(flag==1)
%     loglog(dist_axis(iplot)/1000,Ddd_m0_ensemble_neg(iplot),'Color','[0.4157, 0.3176, 0.6392 ]',...
%         'Linewidth',1);
% hold on
loglog(dist_axis(iplot)/1000,Ddd_m0_ensemble_neg(iplot),'Color','[0.5, 0, 0.1484 ]',...
    'Linewidth',1,'LineStyle','--');


end
%a.Color(4)=0.3;
end
