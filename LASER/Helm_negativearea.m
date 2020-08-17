
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



    Moment=2;
mattitle=sprintf('LASER_all_ensemble_Moment_%d.mat',Moment);
load(mattitle);

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
figure

loglog(dist_axis(iplot),Ddd2M_mag(iplot)/2,'g')
hold on
loglog(dist_axis(iplot),Ddd_m0_ensemble(iplot),'g:');


axis square  

xlabel('r(m)')
lgd.FontSize = 16;

inegative=find(Ddd2M_mag/2>Ddd_m0_ensemble);

ybars(1)=dist_axis(min(inegative));
ybars(2)=dist_axis(max(inegative));

patch([ybars(1) ybars(1), ybars(2) ybars(2)],[min(ylim) max(ylim) max(ylim) min(ylim)], [0.8 0.8 0.8])
alpha(.5)