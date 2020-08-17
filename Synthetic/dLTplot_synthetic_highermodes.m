 %plotting the reconstructed azimuthal Fourier modes of
%velocity structure functions, and check if they agree with the true
%answers as well as the decorrelation limits. The figure is not excatly
%formatted as in the paper. The legends etc. are close, but there is one
%major difference: we are presenting the outcomes at the largest r that is
%spanned in the data here, and do NOT cut it off at around r= 84 km, which
%is the correct one argued in paper in section 5.(a). This is mainly to
%demonstrate that the outcomes are indeed erroneous at r>84 km in this
%example.

%clear all
close all
rstart=0;

load('Synthetictraj.mat')

%Check dLT
doPlotSetup
Mode=0;
mattitle=sprintf('Synthetic_ensemble_Mode_%d.mat',Mode);
load(mattitle);
dlt_m0_ensemble=dlt_M_sin_ensemble;
dll_m0_ensemble=dll_M_cos_ensemble;
dtt_m0_ensemble=dtt_M_cos_ensemble;

Mode=2;
mattitle=sprintf('Synthetic_ensemble_Mode_%d.mat',Mode);
load(mattitle);

dll_M_ensemble=sqrt(dll_M_cos_ensemble.^2+dll_M_sin_ensemble.^2);
dtt_M_ensemble=sqrt(dtt_M_cos_ensemble.^2+dtt_M_sin_ensemble.^2);
dlt_M_ensemble=sqrt(dlt_M_cos_ensemble.^2+dlt_M_sin_ensemble.^2);

dll_m2_amp=dll_M_ensemble;
dtt_m2_amp=dtt_M_ensemble;
dlt_m2_amp=dlt_M_ensemble;

Mode=4;
mattitle=sprintf('Synthetic_ensemble_Mode_%d.mat',Mode);
load(mattitle);

dll_M_ensemble=sqrt(dll_M_cos_ensemble.^2+dll_M_sin_ensemble.^2);
dtt_M_ensemble=sqrt(dtt_M_cos_ensemble.^2+dtt_M_sin_ensemble.^2);
dlt_M_ensemble=sqrt(dlt_M_cos_ensemble.^2+dlt_M_sin_ensemble.^2);

dll_m4_amp=dll_M_ensemble;
dtt_m4_amp=dtt_M_ensemble;
dlt_m4_amp=dlt_M_ensemble;


iplot=(N/2):N;

irecon=1:length(dist_axis);
% irecon=find(dist_axis<0.9*rcutoff); %Uncaption this if you want to plot up
% to rcutoff, which is ~84 km. As explained in paper section 5(a), or as
% can be deduced from Datagen.m at line 200, beyond rcutoff, there will be
% signifant angle gaps in the distribution of separation vectors.

figure(5)
subplot(1,3,1)
plot(xr(iplot)/1000,dL2_m0_test,'Color','[0.5, 0, 0.1484 ]')
hold on
plot(dist_axis(irecon)/1000,dll_m0_ensemble(irecon),'Color','[0.5, 0, 0.1484 ]','Marker','o','LineStyle', 'none' )


subplot(1,3,1)
plot(xr(iplot)/1000,sqrt(dL2_m2_cos_test.^2+dL2_m2_sin_test.^2),'Color','[0.2549, 0.6706, 0.3607]')
hold on
plot(dist_axis(irecon)/1000,dll_m2_amp(irecon),'Color','[0.2549, 0.6706, 0.3607]','Marker','o','LineStyle', 'none' )

subplot(1,3,1)

plot(xr(iplot)/1000,sqrt(dL2_m4_cos_test.^2+dL2_m4_sin_test.^2),'Color','[0.5, 0.4883, 0.7294 ]')
hold on
plot(dist_axis(irecon)/1000,dll_m4_amp(irecon),'Color','[0.5, 0.4883, 0.7294 ]','Marker','o','LineStyle', 'none' )



%figure(2)
subplot(1,3,2)
plot(xr(iplot)/1000,dT2_m0_test,'Color','[0.5, 0, 0.1484 ]')
hold on
plot(dist_axis(irecon)/1000,dtt_m0_ensemble(irecon),'Color','[0.5, 0, 0.1484 ]','Marker','o','LineStyle', 'none' )

subplot(1,3,2)
plot(xr(iplot)/1000,sqrt(dT2_m2_cos_test.^2+dT2_m2_sin_test.^2),'Color','[0.2549, 0.6706, 0.3607]')
hold on
plot(dist_axis(irecon)/1000,dtt_m2_amp(irecon),'Color','[0.2549, 0.6706, 0.3607]','Marker','o','LineStyle', 'none' )

subplot(1,3,2)

plot(xr(iplot)/1000,sqrt(dT2_m4_cos_test.^2+dT2_m4_sin_test.^2),'Color','[0.5, 0.4883, 0.7294 ]')
hold on
plot(dist_axis(irecon)/1000,dtt_m4_amp(irecon),'Color','[0.5, 0.4883, 0.7294 ]','Marker','o','LineStyle', 'none' )

%figure(3)
subplot(1,3,3)
plot(xr(iplot)/1000,abs(dLT_m0_test),'Color','[0.5, 0, 0.1484 ]')
hold on
plot(dist_axis(irecon)/1000,abs(dlt_m0_ensemble(irecon)),'Color','[0.5, 0, 0.1484 ]','Marker','o','LineStyle', 'none' )

subplot(1,3,3)
plot(xr(iplot)/1000,sqrt(dLT_m2_sin_test.^2+dLT_m2_cos_test.^2),'Color','[0.2549, 0.6706, 0.3607]')
hold on
plot(dist_axis(irecon)/1000,dlt_m2_amp(irecon),'Color','[0.2549, 0.6706, 0.3607]','Marker','o','LineStyle', 'none' )

subplot(1,3,3)

plot(xr(iplot)/1000,sqrt(dLT_m4_sin_test.^2+dLT_m4_cos_test.^2),'Color','[0.5, 0.4883, 0.7294 ]')
hold on
plot(dist_axis(irecon)/1000,dlt_m4_amp(irecon),'Color','[0.5, 0.4883, 0.7294 ]','Marker','o','LineStyle', 'none' )


%Mark the limiting behavior 
load Synthetic_CuCvCuv.mat

SLL0=0.*dist_axis+Cu0+Cv0;
SLL2=0.*dist_axis+sqrt((Cu0-Cv0)^2+(2*Cuv0)^2);

figure(5);
subplot(1,3,1)
plot(dist_axis/1000,SLL0,'Color','[0.5, 0, 0.1484 ]','Linestyle','-.')
hold on
subplot(1,3,1)
plot(dist_axis/1000,SLL2,'Color','[0.2549, 0.6706, 0.3607]','Linestyle','-.')
hold on

% savefig('figure1.fig')

figure(5);subplot(1,3,3) 
SLT2=0.*dist_axis+sqrt((Cu0-Cv0)^2+(2*Cuv0)^2);

plot(dist_axis/1000,SLT2,'Color','[0.2549, 0.6706, 0.3607]','Linestyle','-.')
% savefig('figure3.fig')


figure(5);
STT0=0.*dist_axis+Cu0+Cv0;
STT2=0.*dist_axis+sqrt((Cu0-Cv0)^2+(2*Cuv0)^2);
subplot(1,3,2)
plot(dist_axis/1000,STT0,'Color','[0.5, 0, 0.1484 ]','Linestyle','-.')
hold on
subplot(1,3,2)
plot(dist_axis/1000,STT2,'Color','[0.2549, 0.6706, 0.3607]','Linestyle','-.')
hold on
% savefig('figure2.fig')

figure(5)
subplot(1,3,1)
pbaspect([1.6 1 1]) 
axis tight
ax = gca;
ax.FontSize = 13; 
xlabel('r [km]');ylabel('[m^2/s^2]') 
legend('$D^{c0}_{LL}$,true',...  
    '$D^{c0}_{LL}$,estimated',...
    '$D^{a2}_{LL}$,true',...    
    '$D^{a2}_{LL}$,estimated',...
    '$D^{a4}_{LL}$,true',...  
    '$D^{a4}_{LL}$,estimated',...
        'Location','northoutside',...
    'Interpreter','latex','AutoUpdate','off',...
    'Location','northwest','Fontsize',16) % R2018b and later
chi=get(gca, 'Children');
%Reverse the stacking order
set(gca, 'Children',flipud(chi))

subplot(1,3,2)
set(gca,'YTickLabel',[]);
pbaspect([1.6 1 1]) 
axis tight
ax = gca;
ax.FontSize = 13; 
xlabel('r [km]'); 
legend('$D^{c0}_{TT}$,true',...  
    '$D^{c0}_{TT}$,estimated',...
    '$D^{a2}_{TT}$,true',...    
    '$D^{a2}_{TT}$,estimated',...
    '$D^{a4}_{TT}$,true',...  
    '$D^{a4}_{TT}$,estimated',...
        'Location','northoutside',...
    'Interpreter','latex','AutoUpdate','off',...
    'Location','northwest','Fontsize',16)  % R2018b and later
chi=get(gca, 'Children');
%Reverse the stacking order
set(gca, 'Children',flipud(chi))

subplot(1,3,3)
set(gca,'YTickLabel',[]);
pbaspect([1.6 1 1]) 
axis tight
ax = gca;
ax.FontSize = 13; 
xlabel('r [km]');
legend('$D^{c0}_{LT}$,true',...  
    '$D^{c0}_{LT}$,estimated',...
    '$D^{a2}_{LT}$,true',...    
    '$D^{a2}_{LT}$,estimated',...
    '$D^{a4}_{LT}$,true',...  
    '$D^{a4}_{LT}$,estimated',...
        'Location','northoutside',...
    'Interpreter','latex','AutoUpdate','off',...
    'Location','northwest','Fontsize',16)  % R2018b and later
chi=get(gca, 'Children');
%Reverse the stacking order
set(gca, 'Children',flipud(chi))

fg1=figure(5);
sameaxes([],[fg1])


