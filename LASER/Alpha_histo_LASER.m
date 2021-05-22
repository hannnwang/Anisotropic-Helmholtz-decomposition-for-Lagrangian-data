%20190916:modified from Alpha_histo.m in GLAD/20190619
%(Change it so that less storage is demanded, but it takes long to run)

clear all
close all
doPlotSetup

load ('LASER_all.mat')


% calculate time series of pair separation
ninterval=1

ndrfts=size(trajmat_X,2)

% ndrfts_sampled=ndrfts%Include all floaters

ntimesnap_total = size(trajmat_X,1);

ntimesnap_end=ntimesnap_total;

timevec=1:ninterval:ntimesnap_end;

ntimesnap_sampled = length(timevec);

% npairs = factorial(ndrfts)/factorial(ndrfts-2)/factorial(2);
npairs = ndrfts*(ndrfts-1)/2;
% npairs_sampled = ndrfts_sampled*(ndrfts_sampled-1)/2;

n=1;
%I have to dump NaN values at this
%step to save memory. Specifically, in X1,Y1,U1,V1 etc.,
%pair number and time are treated equally. They are all 1D arrays now. This
%makes it easier to dump NaN values.



% %Only keep the drifters whose positions fall into a circle
circlespec

Rearth=6.3782*10^6;

converted_X=trajmat_X.*cosd(trajmat_Y)*Rearth*pi/180;
converted_Y=trajmat_Y*Rearth*pi/180;

Radius_test=sqrt((converted_X-center_x).^2+(converted_Y-center_y).^2);

trajmat_X(Radius_test>radius)=NaN;
trajmat_Y(Radius_test>radius)=NaN;
trajmat_U(Radius_test>radius)=NaN;
trajmat_V(Radius_test>radius)=NaN;
% 
for i=1:ndrfts-1
    for j = i+1:ndrfts
            
            %(X1, Y1 are in non-standard units)
    
        for it=1:ninterval:ntimesnap_end
            if ~isnan(trajmat_X(it,i)) && ~isnan(trajmat_X(it,j))
               X1(n) = trajmat_X(it,i);
                X2(n) = trajmat_X(it,j);
                Y1(n) = trajmat_Y(it,i);
                Y2(n) = trajmat_Y(it,j);
                n = n+1;
          end
        end
    end
end

disp('X1,X2,etc. are formed.')


dr=sqrt((X1.*cosd(Y1)-X2.*cosd(Y2)).^2+...
    (Y1-Y2).^2)*Rearth*pi/180;

Cosine=X1.*cosd(Y1)-X2.*cosd(Y2);
Sine = Y1 - Y2;
clearvars X1 X2 Y1 Y2
Fnormalization= sqrt(Cosine.^2+Sine.^2);
% normalize to unit vectors
Cosine = Cosine./Fnormalization;
Sine = Sine./Fnormalization;
clearvars Fnormalization

%Added on 20191105
% THIS IS WRONG!!!!
% Alpha=acos(Cosine);
% Alpha(Alpha==pi)=0;%There is one Alpha that happens to be pi!
Alpha=atan2(Sine,Cosine);
Alpha(Alpha<0)=Alpha(Alpha<0)+pi;%Using that the function is even

% alpha_apair=acos(Cosine);

%load dist_axis from HPC outcomes
load('distaxis_LASER.mat','dist_axis','dist_bin');
%(This is generated in dist_axis_comp.m)
nangleaxis=32;
	nanglewall=nangleaxis+1;%Number of angle bin walls	
	angle_bin=linspace(0,pi,nanglewall);	
	angle_axis = 0.5*(angle_bin(1:end-1) + angle_bin(2:end));
    
nanglecount=zeros(length(dist_axis),nangleaxis);

%dr is 1D array. Probably shouldn't search in 2D.
%Doesn't cause apparent harm so far.-20191113
for ix=1:length(dist_axis)
    [ivalid_t,ivalid_pair] = find(dr<dist_bin(ix+1) & ...
        dr>=dist_bin(ix));
    npair_atx=length(ivalid_pair);%npair_atx is already ensured to be larger than 1 from dist_axis_comp.m   
        %Get valid angles at this ix
        for it=1:npair_atx
            alpha_cur=Alpha(ivalid_t(it),ivalid_pair(it));
            
	            if alpha_cur==pi	
	                ialpha=1;	
	            else	
	         ialpha=find(angle_bin<=alpha_cur,1,'last');	
	            end	
	            	
	        nanglecount(ix,ialpha)=nanglecount(ix,ialpha)+1;
            
        end
end

disp('Number of angle bins:')
nangleaxis

%Save everything except for the really big variables
clearvars Alpha Cosine dr Sine

save('LASER_all_alpha.mat')

%check if there are "gaps"
ixvec(1)=1;
ixvec(end+1)=round(length(dist_axis)/2);
ixvec(end+1)=round(length(dist_axis));
    figure(1)
for ifig=1:length(ixvec)
    subplot(2,2,ifig)
    plot(1:nanglewall-1,nanglecount(ixvec(ifig),:))
titname=sprintf('at r=%.1f',dist_axis(ixvec(ifig)));
title(titname)
hold on
end

