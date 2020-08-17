%Calculate the average of u^2, v^2 and uv from all drifter
%observations. This is useful in estimating the decorrelation limits, as
%shown in eq.(40)(41) in paper.
clear all
load ('Synthetictraj.mat')

clearvars trajmat_X trajmat_Y 
trajmat_U=reshape(trajmat_U,size(trajmat_U,1)*size(trajmat_U,2),1);
trajmat_V=reshape(trajmat_V,size(trajmat_U,1)*size(trajmat_U,2),1);

Cu0=nanmean(trajmat_U.^2);
Cv0=nanmean(trajmat_V.^2);
Cuv0=nanmean(trajmat_U.*trajmat_V);

save('Synthetic_CuCvCuv.mat','Cu0','Cv0','Cuv0')