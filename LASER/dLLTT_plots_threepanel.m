clear all
close all
rstart=0;


%Check dLT
doPlotSetup
Moment=0;
mattitle=sprintf('LASER_all_ensemble_Moment_%d.mat',Moment);
load(mattitle);

rend=3*10^5;

%Check dLT
doPlotSetup
Moment=0;
mattitle=sprintf('LASER_all_ensemble_Moment_%d.mat',Moment);
load(mattitle);
dlt_m0_ensemble=dlt2M_sin_ensemble;
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

Moment=4;
mattitle=sprintf('LASER_all_ensemble_Moment_%d.mat',Moment);
load(mattitle);

dll2M_ensemble=sqrt(dll2M_cos_ensemble.^2+dll2M_sin_ensemble.^2);
dtt2M_ensemble=sqrt(dtt2M_cos_ensemble.^2+dtt2M_sin_ensemble.^2);
dlt2M_ensemble=sqrt(dlt2M_cos_ensemble.^2+dlt2M_sin_ensemble.^2);

dll_m4_mag=dll2M_ensemble;
dtt_m4_mag=dtt2M_ensemble;
dlt_m4_mag=dlt2M_ensemble;


irecon=find(dist_axis<rend);

figure(1)
subplot(1,3,1)
loglog(dist_axis(irecon)/1000,dll_m0_ensemble(irecon),'Color','[0.5, 0, 0.1484 ]')

hold on
subplot(1,3,1)
loglog(dist_axis(irecon)/1000,dll_m2_mag(irecon),'Color','[0.2549, 0.6706, 0.3607]')
hold on

subplot(1,3,1)
loglog(dist_axis(irecon)/1000,dll_m4_mag(irecon),'Color','[0.5, 0.4883, 0.7294 ]')

hold on


%figure(2)
subplot(1,3,2)
loglog(dist_axis(irecon)/1000,dtt_m0_ensemble(irecon),'Color','[0.5, 0, 0.1484 ]')
hold on

subplot(1,3,2)
loglog(dist_axis(irecon)/1000,dtt_m2_mag(irecon),'Color','[0.2549, 0.6706, 0.3607]')
hold on

subplot(1,3,2)
loglog(dist_axis(irecon)/1000,dtt_m4_mag(irecon),'Color','[0.5, 0.4883, 0.7294 ]')
hold on

%figure(3)
subplot(1,3,3)
loglog(dist_axis(irecon)/1000,abs(dlt_m0_ensemble(irecon)),'Color','[0.5, 0, 0.1484 ]')
hold on

subplot(1,3,3)
loglog(dist_axis(irecon)/1000,dlt_m2_mag(irecon),'Color','[0.2549, 0.6706, 0.3607]')
hold on

subplot(1,3,3)
loglog(dist_axis(irecon)/1000,dlt_m4_mag(irecon),'Color','[0.5, 0.4883, 0.7294 ]')
hold on

%Mark the limiting behavior 
load LASER_CuCvCuv.mat

SLL0=0.*dist_axis+Cu0+Cv0;
SLL2=0.*dist_axis+sqrt((Cu0-Cv0)^2+(2*Cuv0)^2);

figure(1);
subplot(1,3,1)
plot(dist_axis/1000,SLL0,'Color','[0.5, 0, 0.1484 ]','Linestyle','-.')
hold on
subplot(1,3,1)
plot(dist_axis/1000,SLL2,'Color','[0.2549, 0.6706, 0.3607]','Linestyle','-.')
hold on

% savefig('figure1.fig')

figure(1);subplot(1,3,3) 
SLT2=0.*dist_axis+sqrt((Cu0-Cv0)^2+(2*Cuv0)^2);

plot(dist_axis/1000,SLT2,'Color','[0.2549, 0.6706, 0.3607]','Linestyle','-.')
% savefig('figure3.fig')


figure(1);
STT0=0.*dist_axis+Cu0+Cv0;
STT2=0.*dist_axis+sqrt((Cu0-Cv0)^2+(2*Cuv0)^2);
subplot(1,3,2)
plot(dist_axis/1000,STT0,'Color','[0.5, 0, 0.1484 ]','Linestyle','-.')
hold on
subplot(1,3,2)
plot(dist_axis/1000,STT2,'Color','[0.2549, 0.6706, 0.3607]','Linestyle','-.')
hold on
% savefig('figure2.fig')

figure(1)
subplot(1,3,1)
pbaspect([1.2 1 1])
axis tight
ax = gca;
ax.FontSize = 13; 
xlabel('r [km]');ylabel('[m^2/s^2]') 
legend('$D^{c0}_{LL}$',...
    '$D^{a2}_{LL}$',...
    '$D^{a4}_{LL}$',...
        'Location','northoutside',...
    'Interpreter','latex','AutoUpdate','off',...
    'Location','northwest','Fontsize',16) % R2018b and later
chi=get(gca, 'Children');
%Reverse the stacking order
set(gca, 'Children',flipud(chi))

subplot(1,3,2)
set(gca,'YTickLabel',[]);
pbaspect([1.2 1 1])
axis tight
ax = gca;
ax.FontSize = 13; 
xlabel('r [km]'); 
legend('$D^{c0}_{TT}$',...
    '$D^{a2}_{TT}$',...
    '$D^{a4}_{TT}$',...
        'Location','northoutside',...
    'Interpreter','latex','AutoUpdate','off',...
    'Location','northwest','Fontsize',16) % R2018b and later
chi=get(gca, 'Children');
%Reverse the stacking order
set(gca, 'Children',flipud(chi))

subplot(1,3,3)
set(gca,'YTickLabel',[]);
pbaspect([1.2 1 1])
axis tight
ax = gca;
ax.FontSize = 13; 
xlabel('r [km]');
legend('$|D^{c0}_{LT}|$',...
    '$D^{a2}_{LT}$',...
    '$D^{a4}_{LT}$',...
        'Location','northoutside',...
    'Interpreter','latex','AutoUpdate','off',...
    'Location','northwest','Fontsize',16) % R2018b and later
chi=get(gca, 'Children');
%Reverse the stacking order
set(gca, 'Children',flipud(chi))

fg=figure(1);
sameaxes([], fg);
subplot(1,3,1)
yticks([10^(-6) 10^(-5) 10^(-4) 10^(-3) 10^(-2) 10^(-1)]);
xticks([10^(-2) 10^(-1) 10^(0) 10^(1) 10^(2)])
axis tight
yl=ylim;
ylim([yl(1) yl(2)*1.5]);
subplot(1,3,2)
xticks([10^(-2) 10^(-1) 10^(0) 10^(1) 10^(2)])
axis tight
yl=ylim;
ylim([yl(1) yl(2)*1.5]);
subplot(1,3,3)
xticks([10^(-2) 10^(-1) 10^(0) 10^(1) 10^(2)])
axis tight
yl=ylim;
ylim([yl(1) yl(2)*1.5]);