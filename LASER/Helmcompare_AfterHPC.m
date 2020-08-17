clear all
close all
doPlotSetup
  Moment=0;
mattitle=sprintf('LASER_all_ensemble_Moment_%d.mat',Moment);
load(mattitle);

fg1=figure(1);
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


figure(1)
subplot(1,3,1)
plot(dist_axis(iplot),Drr_m0_ensemble(iplot),'Color','[0.4940, 0.1840, 0.5560]');
hold on
plot(dist_axis(iplot),Ddd_m0_ensemble(iplot),'Color','[0.4660, 0.6740, 0.1880]');

legend('Rotational','Divergent',...
    'Interpreter','latex','Location','northwest','Autoupdate','Off')
axis square  

title('0th')
xlabel('r(m)')
lgd.FontSize = 16;

figure(2)
subplot(1,3,1)
loglog(dist_axis(iplot),Drr_m0_ensemble(iplot),'Color','[0.4940, 0.1840, 0.5560]');
hold on
loglog(dist_axis(iplot),Ddd_m0_ensemble(iplot),'Color','[0.4660, 0.6740, 0.1880]');

legend('Rotational','Divergent',...
    'Interpreter','latex','Location','northwest','Autoupdate','Off')
axis square  

title('0th')
xlabel('r(m)')
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
figure(1)
subplot(1,3,imoment+1)
plot(dist_axis(iplot),Drr2M_mag(iplot),'Color','[0.4940, 0.1840, 0.5560]')
hold on
plot(dist_axis(iplot),Ddd2M_mag(iplot),'Color','[0.4660, 0.6740, 0.1880]')

legend('Rotational','Divergent',...
    'Interpreter','latex','Location','northwest','Autoupdate','Off')
    fm01=plot(dist_axis(iplot),Drr_m0_ensemble(iplot),'Color','[0.4940, 0.1840, 0.5560]',...
    'Linestyle',':','Linewidth',1);
hold on
fm02=plot(dist_axis(iplot),Ddd_m0_ensemble(iplot),'Color','[0.4660, 0.6740, 0.1880]',...
    'Linestyle',':','Linewidth',1);
% fm01.Color(4) = 0.5;
% fm02.Color(4) = 0.5;

axis square  

if Moment==2
    figtitle=sprintf('2nd',Moment);
else
figtitle=sprintf('%d th',Moment);
end
title(figtitle)
xlabel('r(m)')
lgd.FontSize = 16;

figure(2)
subplot(1,3,imoment+1)
loglog(dist_axis(iplot),Drr2M_mag(iplot),'Color','[0.4940, 0.1840, 0.5560]')
hold on
loglog(dist_axis(iplot),Ddd2M_mag(iplot),'Color','[0.4660, 0.6740, 0.1880]')

legend('Rotational','Divergent',...
    'Interpreter','latex','Location','northwest','Autoupdate','Off')
fm021=loglog(dist_axis(iplot),Drr_m0_ensemble(iplot),'Color','[0.4940, 0.1840, 0.5560]',...
    'Linestyle',':','Linewidth',1);
hold on
fm022=loglog(dist_axis(iplot),Ddd_m0_ensemble(iplot),'Color','[0.4660, 0.6740, 0.1880]',...
    'Linestyle',':','Linewidth',1);
% fm021.Color(4) = 0.5;
% fm022.Color(4) = 0.5;

axis square  
if Moment==2
    figtitle=sprintf('2nd mode',Moment);
else
figtitle=sprintf('%d th mode',Moment);
end
title(figtitle)
xlabel('r(m)')
lgd.FontSize = 16;

% figure(3)
% subplot(1,2,imoment)
% plot(dist_axis(iplot),Drr2M_mag(iplot)./Drr_m0_ensemble(iplot),'Color','[0.4940, 0.1840, 0.5560]')
% hold on
% plot(dist_axis(iplot),Ddd2M_mag(iplot)./Ddd_m0_ensemble(iplot),'Color','[0.4660, 0.6740, 0.1880]')
% 
% legend('$<u_r^2> ratio$','$<u_d^2> ratio$',...
%     'Interpreter','latex','Location','northwest','Autoupdate','Off')
%   
% axis square  
% 
% xlabel('r(m)')
% figtitle=sprintf('%d th mode',Moment);
% title(figtitle)
% lgd.FontSize = 16;

end
f1=figure(1);
sameaxes([],[f1]);

f2=figure(2);
sameaxes([],[f2]);

figure(1)
savefig('Helm1.fig')


figure(2)
savefig('Helm2.fig')

