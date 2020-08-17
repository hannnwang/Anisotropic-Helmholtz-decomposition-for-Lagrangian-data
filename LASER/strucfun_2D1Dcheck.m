%UNRESOLVED: BUG WITH ALPHA DISCOVERED IN OCTOBER, 2019 HASN'T BEEN CHECKED
%IN THIS FILE YET.
% 
% %20190920
% %Check if my 2D structure functions are consistent with 1D structure
% %functions. The main motivation for this is that dlt_2D seems positive in the 2D plot
% clear all
% %close all
% doPlotSetup
% load LASER_2Dstrucfun.mat
% 
% %Check the zeroth mode
% for i=1:length(dist_axis)
% dll_m0_from2D(i)=mean(dll_2D(i,:));
% dtt_m0_from2D(i)=mean(dtt_2D(i,:));
% dlt_m0_from2D(i)=mean(dlt_2D(i,:));
% end
% 
%   Moment=0;
% mattitle=sprintf('LASER_all_ensemble_Moment_%d.mat',Moment);
% load(mattitle);
% 
% %20190920:IMPORTANT: All the higher mode terms calculated in 
% %Helm_arbitrarymode_angleweighted_LASER.m are MISSING a factor of 2
% %here! I have not changed that program yet because that's already run on HPC
% %by the time I find this out. Please NOTE THAT THE RESULTS NEED TO BE
% %CORRECTED IN ALL THE LATER TREATMENTS!
% %Belo
% if M~=0
%     dll2M_mu_ensemble=dll2M_mu_ensemble*2;
%     dll2M_nu_ensemble=dll2M_nu_ensemble*2;
%     dtt2M_mu_ensemble=dtt2M_mu_ensemble*2;
%     dtt2M_nu_ensemble=dtt2M_nu_ensemble*2;
%     dlt2M_mu_ensemble=dlt2M_mu_ensemble*2;
%     dlt2M_nu_ensemble=dlt2M_nu_ensemble*2;
% end
% dll_m0_1D=dll2M_mu_ensemble;    
% dtt_m0_1D=dtt2M_mu_ensemble;    
% dlt_m0_1D=dlt2M_mu_ensemble;  
% figure
% subplot(2,2,1)
% semilogx(dist_axis,dll_m0_1D)
% hold on
% semilogx(dist_axis,dll_m0_from2D)
% legend('dll_0, from 1D','dll_0, from 2D')
% subplot(2,2,2)
% semilogx(dist_axis,dtt_m0_1D)
% hold on
% semilogx(dist_axis,dtt_m0_from2D)
% legend('dtt_0, from 1D','dtt_0, from 2D')
% subplot(2,2,3)
% semilogx(dist_axis,dlt_m0_1D)
% hold on
% semilogx(dist_axis,dlt_m0_from2D)
% legend('dlt_0, from 1D','dlt_0, from 2D')
% 
% 
% %Check the second mode
% theta_all=((1:180)-0.5)*pi/180;
% for i=1:length(dist_axis)
% dll_m2_cos_from2D(i)=2*mean(dll_2D(i,:).*cos(2*theta_all));
% dtt_m2_cos_from2D(i)=2*mean(dtt_2D(i,:).*cos(2*theta_all));
% dlt_m2_cos_from2D(i)=2*mean(dlt_2D(i,:).*cos(2*theta_all));
% 
% dll_m2_sin_from2D(i)=2*mean(dll_2D(i,:).*sin(2*theta_all));
% dtt_m2_sin_from2D(i)=2*mean(dtt_2D(i,:).*sin(2*theta_all));
% dlt_m2_sin_from2D(i)=2*mean(dlt_2D(i,:).*sin(2*theta_all));
% 
% end
% 
% dll_m2_from2D=sqrt(dll_m2_cos_from2D.^2+dll_m2_sin_from2D.^2);
% dtt_m2_from2D=sqrt(dtt_m2_cos_from2D.^2+dtt_m2_sin_from2D.^2);
% dlt_m2_from2D=sqrt(dlt_m2_cos_from2D.^2+dlt_m2_sin_from2D.^2);
% Moment=2;
% mattitle=sprintf('LASER_all_ensemble_Moment_%d.mat',Moment);
% load(mattitle);
% 
% %20190920:IMPORTANT: All the higher mode terms calculated in 
% %Helm_arbitrarymode_angleweighted_LASER.m are MISSING a factor of 2
% %here! I have not changed that program yet because that's already run on HPC
% %by the time I find this out. Please NOTE THAT THE RESULTS NEED TO BE
% %CORRECTED IN ALL THE LATER TREATMENTS!
% %Belo
% if M~=0
%     dll2M_mu_ensemble=dll2M_mu_ensemble*2;
%     dll2M_nu_ensemble=dll2M_nu_ensemble*2;
%     dtt2M_mu_ensemble=dtt2M_mu_ensemble*2;
%     dtt2M_nu_ensemble=dtt2M_nu_ensemble*2;
%     dlt2M_mu_ensemble=dlt2M_mu_ensemble*2;
%     dlt2M_nu_ensemble=dlt2M_nu_ensemble*2;
% end
% 
% dll2M_ensemble=sqrt(dll2M_mu_ensemble.^2+dll2M_nu_ensemble.^2);
% dtt2M_ensemble=sqrt(dtt2M_mu_ensemble.^2+dtt2M_nu_ensemble.^2);
% dlt2M_ensemble=sqrt(dlt2M_mu_ensemble.^2+dlt2M_nu_ensemble.^2);
% 
% dll_m2_1D=dll2M_ensemble;
% dtt_m2_1D=dtt2M_ensemble;
% dlt_m2_1D=dlt2M_ensemble;
% 
% figure
% subplot(2,2,1)
% semilogx(dist_axis,dll_m2_1D)
% hold on
% semilogx(dist_axis,dll_m2_from2D)
% legend('dll_m2, from 1D','dll_m2, from 2D')
% subplot(2,2,2)
% semilogx(dist_axis,dtt_m2_1D)
% hold on
% semilogx(dist_axis,dtt_m2_from2D)
% legend('dtt_m2, from 1D','dtt_m2, from 2D')
% subplot(2,2,3)
% semilogx(dist_axis,dlt_m2_1D)
% hold on
% semilogx(dist_axis,dlt_m2_from2D)
% legend('dlt_m2, from 1D','dlt_m2, from 2D')
% 
% 
