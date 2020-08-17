%Calculate the binning in r
clear all
close all
doPlotSetup

load ('Synthetictraj.mat')


% calculate time series of pair separation
ninterval=1

nsnapshot_total = size(trajmat_X,1);

nsnapshot_end=nsnapshot_total;

timevec=1:ninterval:nsnapshot_end;%sampled snapshots -- can modify ninterval 
%or nsnapshot_end to include less snapshots. 

nsnapshot_sampled = length(timevec);

n=1;
%Form the separation vectors.
%Note: this code is intended to be
%fast and memory-saving. 
%The trick is to stretch all the inputs into 1D arrays first, and try to
%take advantage of fast matrix manipulations built in Matlab.

tic

for i=1:ndrfts-1
    for j = i+1:ndrfts
            
            %(X1, Y1 are in non-standard units)
    
        for it=1:ninterval:nsnapshot_end
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

%These are all matrix manipulations, which are fast in Matlab!
dr=sqrt((X1-X2).^2+(Y1-Y2).^2);

disp('Hopefully all the big-memory variables are done now.')

clearvars X1 X2 Y1 Y2

tic
xcutoff=xr(end);%Dump the ones larger than xcutoff

xstart=(x(2)-x(1))/100;

dr=dr(dr<=xcutoff);

dr=dr(dr>0);%Deleting "self pairs"

dr=sort(dr);

nallpairs=length(dr);

nBin=length(xr)/2;%number of bins
%Now, set the bin faces so that approximately the number of drifters
%falling into each bin are the same.
neachBin_r=(nallpairs+nBin-1)/nBin-1; 

fBin_r_index=1:neachBin_r:nallpairs;
fBin_r_index=round(fBin_r_index);

dist_bin=dr(fBin_r_index);%bin faces

%Bin centers 
dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));

if length(dist_axis)~=length(unique(dist_axis))||dist_axis(1)==0
    error('it turns out there are bins with no drifters')
end

toc
disp('dist_axis calculated.')
%check if the floaters are evenly distributed in each bin
ncount_check=histcounts(dr,dist_bin);% each element in ncount_check shows how many pairs are there in that bin.

clearvars dr
mattitle=sprintf('Synthetic_distaxis.mat');

save(mattitle) 