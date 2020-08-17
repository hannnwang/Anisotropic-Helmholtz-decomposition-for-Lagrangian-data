
%  domain=polyshape([-92.6035 -85.1646 -85.6233 -93.3796],...
%     [25.6275 26.4706 30.9037 30.0657]);
% 
% figure
% plot(domain)
% hold on
% xplot=reshape(trajmat_X,[1,5373*311]);
% yplot=reshape(trajmat_Y,[1,5373*311]);
% scatter(xplot,yplot);


xv=[-92.6035 -85.1646 -85.6233 -93.3796];
yv=[25.6275 26.4706 30.9037 30.0657];

iindomain=inpolygon(trajmat_X,trajmat_Y,xv,yv);

trajmat_U(iindomain~=true)=NaN;
trajmat_V(iindomain~=true)=NaN;
trajmat_X(iindomain~=true)=NaN;
trajmat_Y(iindomain~=true)=NaN;