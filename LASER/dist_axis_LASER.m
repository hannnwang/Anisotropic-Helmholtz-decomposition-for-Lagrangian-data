%20200120: enforce that the smallest bin center is larger than 15 m.

%20190914:modified from dist_axis_comp.m in GLAD/20190619
%(Change it so that less storage is demanded, but it takes long to run)

%20190615:Calculate the dist_axis according to the distribution of floaters
%so that the number of floaters would be the same in each bin.

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

tic
%Only keep the drifters whose positions fall into a circle
circlespec

Rearth=6.3782*10^6;

converted_X=trajmat_X.*cosd(trajmat_Y)*Rearth*pi/180;
converted_Y=trajmat_Y*Rearth*pi/180;

Radius_test=sqrt((converted_X-center_x).^2+(converted_Y-center_y).^2);

trajmat_X(Radius_test>radius)=NaN;
trajmat_Y(Radius_test>radius)=NaN;
trajmat_U(Radius_test>radius)=NaN;
trajmat_V(Radius_test>radius)=NaN;

n=1;
% I have to dump NaN values at this
%step to save memory.Specifically, in X1,Y1,U1,V1 etc.,
%pair number and time are treated equally. They are all 1D arrays now. This
%makes it easier to dump NaN values.

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
toc
disp('X1,X2,etc. are formed.')


dr=sqrt((X1.*cosd(Y1)-X2.*cosd(Y2)).^2+...
    (Y1-Y2).^2)*Rearth*pi/180;

disp('Hopefully all the big-momory variables are done now.')

clearvars X1 X2 Y1 Y2

 %Dump the ones at too large a scale
xcutoff=2*10^5;

dr=dr(dr<=xcutoff);

% %Dump the ones less than 0.1 meter
dr=dr(dr>=0.1);

dr=sort(dr);


nallvalid=sum(~isnan(dr));
nBin=min(1000,round((nallvalid/(5000))));

nallpairs=length(dr);

neachBin_r=(nallpairs+nBin-1)/nBin-1; %nBin*(x+1)-(nBin-1)=nallpairs
%don't round neachBin in order for fBin_index to reach nallpairs.

%round fBin_index instead
fBin_r_index=1:neachBin_r:nallpairs;
fBin_r_index=round(fBin_r_index);

dist_bin=dr(fBin_r_index);

%Prevent bins centered at less than 15 m.
idelete=find(dist_bin<30,1,'last');
if idelete>1
    dist_bin(2:idelete)=[];
end

%Bin faces and centers used for actual structure function calculations
dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));


if length(dist_axis)~=length(unique(dist_axis))||dist_axis(1)==0
    error('it turns out there are empty bins')
end


toc
disp('dist_axis calculated.')
%check if the floaters are evenly distributed in each bin
ncount_check=histcounts(dr,dist_bin);

clearvars dr
mattitle=sprintf('distaxis_LASER.mat');

save(mattitle)%Saving this for debugging; would be replaced if the next steps run without error warning.
