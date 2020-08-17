clear all
load ('LASER_all.mat')

% %Only keep the drifters whose positions fall into a circle
% circlespec

Rearth=6.3782*10^6;
% 
% converted_X=trajmat_X.*cosd(trajmat_Y)*Rearth*pi/180;
% converted_Y=trajmat_Y*Rearth*pi/180;
% 
% Radius_test=sqrt((converted_X-center_x).^2+(converted_Y-center_y).^2);
% 
% trajmat_U(Radius_test>radius)=NaN;
% trajmat_V(Radius_test>radius)=NaN;

clearvars trajmat_X trajmat_Y laser_all

trajmat_U=reshape(trajmat_U,9729*1432,1);
trajmat_V=reshape(trajmat_V,9729*1432,1);

Cu0=nanmean(trajmat_U.^2);
Cv0=nanmean(trajmat_V.^2);
Cuv0=nanmean(trajmat_U.*trajmat_V);

save('LASER_CuCvCuv.mat','Cu0','Cv0','Cuv0')