%UNRESOLVED: BUG WITH ALPHA DISCOVERED IN OCTOBER, 2019 HASN'T BEEN CHECKED
%IN THIS FILE YET.
% 
% %20190920
% %See how many modes are required to reconstruct 2D dll
% clear all
% close all
% load LASER_2Dstrucfun.mat
% 
% 
% rcut=5000;
% 
% 
% vtitle=sprintf('struct2Dreconstr_r%d.avi',rcut);
% 
% v = VideoWriter(vtitle);
% v.Quality = 100;
% open(v);
% 
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
% dll_stretch_true=dll_2D(:);
% figure(1)
% subplot(1,2,1)
% polarscatter(th(iscatter),r(iscatter),sz,dll_stretch_true(iscatter))
% colorbar
% title('Dll(r,\alpha), true answer')
% 
% %Reconstruct from the zeroth mode
% for i=1:length(dist_axis)
% dll_m0_from2D(i)=mean(dll_2D(i,:));
% dtt_m0_from2D(i)=mean(dtt_2D(i,:));
% dlt_m0_from2D(i)=mean(dlt_2D(i,:));
% end
% dll_m0_2D=zeros(length(dist_axis),nangle)*NaN;
% for i=1:length(dist_axis)
%     dll_m0_2D(i,:)=dll_m0_from2D(i);
% end
% 
% dll_2Dreconstr=dll_m0_2D;
% dll_stretch=dll_2Dreconstr(:);
% 
% Error_dll(1)=rms(dll_stretch_true(iscatter)-dll_stretch(iscatter));
% 
% figure(1)
% subplot(1,2,2)
% polarscatter(th(iscatter),r(iscatter),sz,dll_stretch(iscatter))
% colorbar
% title('Dll(r,\alpha), reconstructed')
% sameaxes()
% 
% hold on
%  frame = getframe(gcf);
%  writeVideo(v,frame);
% 
% %Reconstruct more modes
% Modevec=2:2:200;
% for ifig=1:length(Modevec)
%     Mode=Modevec(ifig);
%   theta_all=((1:180)-0.5)*pi/180;
% for i=1:length(dist_axis)
% dll_mM_cos_from2D(i)=2*mean(dll_2D(i,:).*cos(Mode*theta_all));
% dtt_mM_cos_from2D(i)=2*mean(dtt_2D(i,:).*cos(Mode*theta_all));
% dlt_mM_cos_from2D(i)=2*mean(dlt_2D(i,:).*cos(Mode*theta_all));
% 
% dll_mM_sin_from2D(i)=2*mean(dll_2D(i,:).*sin(Mode*theta_all));
% dtt_mM_sin_from2D(i)=2*mean(dtt_2D(i,:).*sin(Mode*theta_all));
% dlt_mM_sin_from2D(i)=2*mean(dlt_2D(i,:).*sin(Mode*theta_all));
% end
% 
% dll_mM_2D=zeros(length(dist_axis),nangle)*NaN;
% for i=1:length(dist_axis)
%     for j=1:nangle
%     dll_mM_2D(i,j)=dll_mM_cos_from2D(i)*cos(Mode*theta_all(j))+...
%         dll_mM_sin_from2D(i)*sin(Mode*theta_all(j));
%     end
% end
% 
% dll_2Dreconstr=dll_2Dreconstr+dll_mM_2D;
% dll_stretch=dll_2Dreconstr(:);
% 
% Error_dll(end+1)=rms(dll_stretch_true(iscatter)-dll_stretch(iscatter));
% 
% figure(1)
% subplot(1,2,2)
% polarscatter(th(iscatter),r(iscatter),sz,dll_stretch(iscatter))
% colorbar
% figtitle=sprintf('Reconstructed up to mode %d', Mode);
% 
% 
% title(figtitle)
% hold on
% 
% sameaxes()
% 
%  frame = getframe(gcf);
%  writeVideo(v,frame);
% 
% end
% hold off
% close(v)
% 
% mattitle=sprintf('Strucfun_2D_Error_dll_rcut%d.mat', rcut);
% 
% save(mattitle,'Error_dll')