%UNRESOLVED: BUG WITH ALPHA DISCOVERED IN OCTOBER, 2019 HASN'T BEEN CHECKED
%IN THIS FILE YET.

% clear all
% %close all
% load LASER_2Dstrucfun.mat
% 
% 
% rcut=5000;
% theta_all=((1:180)-0.5)*pi/180;
% [Theta_all,Axis_all]=meshgrid(theta_all,dist_axis);
% idist_axis=find(dist_axis<rcut);
% len=length(idist_axis);
% iscatter=zeros(1,len*180)*NaN;
% for i=1:180
%     ind=(i-1)*len;
%     iscatter(ind+1:ind+len)=idist_axis+(i-1)*1000;
% end
% iscatter=iscatter';
% th=Theta_all(:);
% r=Axis_all(:);
% sz=5;
% dll_scaller=dll_2D(:);
% figure
% polarscatter(th(iscatter),r(iscatter),sz,dll_scaller(iscatter))
% colorbar
% title('Dll(r,\alpha)')
% rticks([dist_axis(1) (rcut-dist_axis(1))/4+dist_axis(1) ....
%     (rcut-dist_axis(1))/2+dist_axis(1) (rcut-dist_axis(1))*3/4+dist_axis(1) rcut])
% 
% figure
% polarscatter(th(iscatter),log(r(iscatter)),sz,dll_scaller(iscatter))
% colorbar
% title('Dll(r,\alpha)')
% rticks([dist_axis(1) (rcut-dist_axis(1))/4+dist_axis(1) ....
%     (rcut-dist_axis(1))/2+dist_axis(1) (rcut-dist_axis(1))*3/4+dist_axis(1) rcut])