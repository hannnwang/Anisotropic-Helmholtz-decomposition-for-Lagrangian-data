%UNRESOLVED: BUG WITH ALPHA DISCOVERED IN OCTOBER, 2019 HASN'T BEEN CHECKED
%IN THIS FILE YET.

%Comparing the Helm decomp, the important change made here is:
%For all anisotropic components, scales smaller than xsubmeso are deemed
%unreliable and removed. 

clear all
close all
doPlotSetup
  Moment=0;
mattitle=sprintf('LASER_all_ensemble_Moment_%d.mat',Moment);
load(mattitle);

fg1=figure(1);
fg2=figure(2);

dist_axis=dist_axis';

xcutoff=2*10^5;
iplot=find(dist_axis<xcutoff);

%moment
TwoM=2*M;
intgrd2M_mu_ensemble=(-TwoM.*dlt2M_mu_ensemble-dll2M_mu_ensemble+...
    dtt2M_mu_ensemble)./dist_axis;
intgrd2M_nu_ensemble=(-TwoM.*dlt2M_nu_ensemble-dll2M_nu_ensemble+...
    dtt2M_nu_ensemble)./dist_axis;

Drr2M_cos_mu_ensemble=dtt2M_mu_ensemble+...
    intRHS_trapz_nobin(dist_axis,intgrd2M_mu_ensemble);

Ddd2M_cos_mu_ensemble=dll2M_mu_ensemble-...
    intRHS_trapz_nobin(dist_axis,intgrd2M_mu_ensemble);

Drr2M_cos_nu_ensemble=dtt2M_nu_ensemble+...
    intRHS_trapz_nobin(dist_axis,intgrd2M_nu_ensemble);

Ddd2M_cos_nu_ensemble=dll2M_nu_ensemble-...
    intRHS_trapz_nobin(dist_axis,intgrd2M_nu_ensemble);

Drr2M_sin_mu_ensemble=Drr2M_cos_nu_ensemble;

Ddd2M_sin_mu_ensemble=Ddd2M_cos_nu_ensemble;
if Moment>0
Drr2M_mag=sqrt(Drr2M_sin_mu_ensemble.^2+Drr2M_cos_mu_ensemble.^2);
Ddd2M_mag=sqrt(Ddd2M_sin_mu_ensemble.^2+Ddd2M_cos_mu_ensemble.^2);
else 
Drr2M_mag=Drr2M_sin_mu_ensemble;
Ddd2M_mag=Ddd2M_sin_mu_ensemble;
end

Drr_m0_ensemble=Drr2M_mag;
Ddd_m0_ensemble=Ddd2M_mag;


Momentvec(1)=2;
Momentvec(2)=4;
for imoment=1:length(Momentvec)
    Moment=Momentvec(imoment);
mattitle=sprintf('LASER_all_ensemble_Moment_%d.mat',Moment);
load(mattitle);

dist_axis=dist_axis';

xcutoff=2*10^5;
iplot=find(dist_axis<xcutoff);

xsubmeso=50000;
iinvalid=find(dist_axis<=xsubmeso);
dll2M_mu_ensemble(iinvalid)=0;
dll2M_nu_ensemble(iinvalid)=0;
dtt2M_mu_ensemble(iinvalid)=0;
dtt2M_nu_ensemble(iinvalid)=0;
dlt2M_mu_ensemble(iinvalid)=0;
dlt2M_nu_ensemble(iinvalid)=0;

%moment
TwoM=2*M;
intgrd2M_mu_ensemble=(-TwoM.*dlt2M_mu_ensemble-dll2M_mu_ensemble+...
    dtt2M_mu_ensemble)./dist_axis;
intgrd2M_nu_ensemble=(-TwoM.*dlt2M_nu_ensemble-dll2M_nu_ensemble+...
    dtt2M_nu_ensemble)./dist_axis;

Drr2M_cos_mu_ensemble=dtt2M_mu_ensemble+...
    intRHS_trapz_nobin(dist_axis,intgrd2M_mu_ensemble);

Ddd2M_cos_mu_ensemble=dll2M_mu_ensemble-...
    intRHS_trapz_nobin(dist_axis,intgrd2M_mu_ensemble);

Drr2M_cos_nu_ensemble=dtt2M_nu_ensemble+...
    intRHS_trapz_nobin(dist_axis,intgrd2M_nu_ensemble);

Ddd2M_cos_nu_ensemble=dll2M_nu_ensemble-...
    intRHS_trapz_nobin(dist_axis,intgrd2M_nu_ensemble);

Drr2M_sin_mu_ensemble=Drr2M_cos_nu_ensemble;

Ddd2M_sin_mu_ensemble=Ddd2M_cos_nu_ensemble;
if Moment>0
Drr2M_mag=sqrt(Drr2M_sin_mu_ensemble.^2+Drr2M_cos_mu_ensemble.^2);
Ddd2M_mag=sqrt(Ddd2M_sin_mu_ensemble.^2+Ddd2M_cos_mu_ensemble.^2);
else 
Drr2M_mag=Drr2M_sin_mu_ensemble;
Ddd2M_mag=Ddd2M_sin_mu_ensemble;
end
figure(1)
subplot(1,2,imoment)
plot(dist_axis(iplot),Drr2M_mag(iplot),'r')
hold on
plot(dist_axis(iplot),Ddd2M_mag(iplot),'g')

legend('Rotational','Divergent',...
    'Interpreter','latex','Location','northwest','Autoupdate','Off')
    fm01=plot(dist_axis(iplot),Drr_m0_ensemble(iplot),'r:');
hold on
fm02=plot(dist_axis(iplot),Ddd_m0_ensemble(iplot),'g:');
fm01.Color(4) = 0.5;
fm02.Color(4) = 0.5;

axis square  

figtitle=sprintf('%d th mode',Moment);
title(figtitle)
xlabel('r(m)')
lgd.FontSize = 16;

figure(2)
subplot(1,2,imoment)
loglog(dist_axis(iplot),Drr2M_mag(iplot),'r')
hold on
loglog(dist_axis(iplot),Ddd2M_mag(iplot),'g')

legend('Rotational','Divergent',...
    'Interpreter','latex','Location','northwest','Autoupdate','Off')
fm021=loglog(dist_axis(iplot),Drr_m0_ensemble(iplot),'r:');
hold on
fm022=loglog(dist_axis(iplot),Ddd_m0_ensemble(iplot),'g:');
fm021.Color(4) = 0.5;
fm022.Color(4) = 0.5;

axis square  

figtitle=sprintf('%d th mode',Moment);
title(figtitle)
xlabel('r(m)')
lgd.FontSize = 16;

% figure(3)
% subplot(1,2,imoment)
% plot(dist_axis(iplot),Drr2M_mag(iplot)./Drr_m0_ensemble(iplot),'r')
% hold on
% plot(dist_axis(iplot),Ddd2M_mag(iplot)./Ddd_m0_ensemble(iplot),'g')
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

