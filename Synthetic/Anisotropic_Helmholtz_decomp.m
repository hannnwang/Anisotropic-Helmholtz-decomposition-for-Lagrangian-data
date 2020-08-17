%Helmholtz decomposition applying eq.(45) in paper 
clear all
%close all
doPlotSetup

%The true answers directly from Datagen.m
load('Synthetictraj.mat')
Drr_m2_test=sqrt(Drr_m2_cos_test.^2+Drr_m2_sin_test.^2);
Drr_m4_test=sqrt(Drr_m4_cos_test.^2+Drr_m4_sin_test.^2);
Ddd_m2_test=sqrt(Ddd_m2_cos_test.^2+Ddd_m2_sin_test.^2);
Ddd_m4_test=sqrt(Ddd_m4_cos_test.^2+Ddd_m4_sin_test.^2);
figure(2)
subplot(1,3,1)
plot(xr(N/2+1:end)/1000,Drr_m0_test(2:end),'Color','[0.0313, 0.3164, 0.6094 ]')
hold on 
plot(xr(N/2+1:end)/1000,Ddd_m0_test(2:end),'Color','[0.4157, 0.3176, 0.6392 ]')
hold on
subplot(1,3,2)
plot(xr(N/2+1:end)/1000,Drr_m2_test(2:end),'Color','[0.4980, 0.8039, 0.7333 ]')
hold on 
plot(xr(N/2+1:end)/1000,Ddd_m2_test(2:end),'Color','[0.9804, 0.6235, 0.7098 ]')
hold on
subplot(1,3,3)
plot(xr(N/2+1:end)/1000,Drr_m4_test(2:end),'Color','[0.4980, 0.8039, 0.7333 ]')
hold on 
plot(xr(N/2+1:end)/1000,Ddd_m4_test(2:end),'Color','[0.9804, 0.6235, 0.7098 ]')
hold on

%Plot Mode 9 first
 Mode=0;
mattitle=sprintf('Synthetic_ensemble_Mode_%d.mat',Mode);
load(mattitle);

dist_axis=dist_axis';

%Directly applying eq.(45) in paper
intgrd_M_cos_ensemble=(-Mode.*dlt_M_sin_ensemble-dll_M_cos_ensemble+...
    dtt_M_cos_ensemble)./dist_axis;
intgrd_M_sin_ensemble=(+Mode.*dlt_M_cos_ensemble-dll_M_sin_ensemble+...
    dtt_M_sin_ensemble)./dist_axis;

Drr_M_cos_ensemble=dtt_M_cos_ensemble+...
    intRHS_trapz_nobin(dist_axis,intgrd_M_cos_ensemble);

Ddd_M_cos_ensemble=dll_M_cos_ensemble-...
    intRHS_trapz_nobin(dist_axis,intgrd_M_cos_ensemble);

Drr_M_sin_ensemble=dtt_M_sin_ensemble+...
    intRHS_trapz_nobin(dist_axis,intgrd_M_sin_ensemble);

Ddd_M_sin_ensemble=dll_M_sin_ensemble-...
    intRHS_trapz_nobin(dist_axis,intgrd_M_sin_ensemble);

if Mode>0
Drr_M_amp=sqrt(Drr_M_sin_ensemble.^2+Drr_M_cos_ensemble.^2);
Ddd_M_amp=sqrt(Ddd_M_sin_ensemble.^2+Ddd_M_cos_ensemble.^2);
else 
Drr_M_amp=Drr_M_sin_ensemble;
Ddd_M_amp=Ddd_M_sin_ensemble;
end

Drr_m0_ensemble=Drr_M_amp;
Ddd_m0_ensemble=Ddd_M_amp;

irecon=1:length(dist_axis);
% irecon=find(dist_axis<0.9*rcutoff); %Uncaption this if you want to plot up
% to rcutoff, which is ~84 km. As explained in paper section 5(a), or as
% can be deduced from Datagen.m at line 200, beyond rcutoff, there will be
% signifant angle gaps in the distribution of separation vectors.

figure(2)
subplot(1,3,1) 
plot(dist_axis(irecon)/1000,Drr_m0_ensemble(irecon),'Color','[0.0313, 0.3164, 0.6094 ]','Marker','o','LineStyle', 'none');
hold on
plot(dist_axis(irecon)/1000,Ddd_m0_ensemble(irecon),'Color','[0.4157, 0.3176, 0.6392 ]','Marker','o','LineStyle', 'none');

pbaspect([1.6 1 1]); axis tight  

lgd.FontSize = 16;

%Repeat the same thing for other modes.
Modevec_plot(1)=2;
Modevec_plot(2)=4;

for imoment=1:length(Modevec_plot)
    Mode=Modevec_plot(imoment);
mattitle=sprintf('Synthetic_ensemble_Mode_%d.mat',Mode);
load(mattitle);

dist_axis=dist_axis';

intgrd_M_cos_ensemble=(-Mode.*dlt_M_sin_ensemble-dll_M_cos_ensemble+...
    dtt_M_cos_ensemble)./dist_axis;
intgrd_M_sin_ensemble=(+Mode.*dlt_M_cos_ensemble-dll_M_sin_ensemble+...
    dtt_M_sin_ensemble)./dist_axis;

Drr_M_cos_ensemble=dtt_M_cos_ensemble+...
    intRHS_trapz_nobin(dist_axis,intgrd_M_cos_ensemble);

Ddd_M_cos_ensemble=dll_M_cos_ensemble-...
    intRHS_trapz_nobin(dist_axis,intgrd_M_cos_ensemble);

Drr_M_sin_ensemble=dtt_M_sin_ensemble+...
    intRHS_trapz_nobin(dist_axis,intgrd_M_sin_ensemble);

Ddd_M_sin_ensemble=dll_M_sin_ensemble-...
    intRHS_trapz_nobin(dist_axis,intgrd_M_sin_ensemble);

if Mode>0
Drr_M_amp=sqrt(Drr_M_sin_ensemble.^2+Drr_M_cos_ensemble.^2);
Ddd_M_amp=sqrt(Ddd_M_sin_ensemble.^2+Ddd_M_cos_ensemble.^2);
else 
Drr_M_amp=Drr_M_sin_ensemble;
Ddd_M_amp=Ddd_M_sin_ensemble;
end

figure(2)
subplot(1,3,imoment+1)
plot(dist_axis(irecon)/1000,Drr_M_amp(irecon),'Color','[0.4980, 0.8039, 0.7333 ]','Marker','o','LineStyle', 'none')
hold on
plot(dist_axis(irecon)/1000,Ddd_M_amp(irecon),'Color','[0.9804, 0.6235, 0.7098 ]','Marker','o','LineStyle', 'none')

pbaspect([1.6 1 1]); axis tight  
if Mode==2
    figtitle=sprintf('2nd mode',Mode);
else
figtitle=sprintf('%d th mode',Mode);
end
lgd.FontSize = 16;

end



irecon=1:length(dist_axis);

figure(2)
subplot(1,3,1) 
pbaspect([1.6 1 1]); axis tight  

subplot(1,3,1)
xlabel('r [km]');ylabel('[m^2/s^2]')  
legend('$D^{c0}_{RR}$,true','$D^{c0}_{DD}$,true','$D^{c0}_{RR}$,estimated',...
    '$D^{c0}_{DD}$,estimated','Interpreter','latex','Location','northwest',...
    'Autoupdate','Off',...
    'Fontsize',16)

subplot(1,3,2)
legend('$D^{a2}_{RR}$,true','$D^{a2}_{DD}$,true','$D^{a2}_{RR}$,estimated',...
    '$D^{a2}_{DD}$,estimated','Interpreter','latex','Location','northwest',...
    'Autoupdate','Off',...
    'Fontsize',16)

xlabel('r [km]');
set(gca,'YTickLabel',[]);
subplot(1,3,3)
xlabel('r [km]');
legend('$D^{a4}_{RR}$,true','$D^{a4}_{DD}$,true','$D^{a4}_{RR}$,estimated',...
    '$D^{a4}_{DD}$,estimated','Interpreter','latex','Location','northwest',...
    'Autoupdate','Off',...
    'Fontsize',16)


set(gca,'YTickLabel',[]);

f2=figure(2);
sameaxes([],[f2]);


