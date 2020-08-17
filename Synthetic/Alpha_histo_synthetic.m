%Count how many drifters are there in each bin in
%alpha at each r. This program is only useful if we use the alternative
%angle weighting approach corresponding to Appendix A in paper, or to the
%program "Arbitrarymode_angleweighted_equalbins_synthetic". You could
%ignore this program if you are not running that approach.

%clear all
%close all
doPlotSetup

load ('Synthetictraj.mat') %generated from dist_axis_synthetic.m

% calculate time series of pair separation. Same as in dist_axis_synthetic.m
ninterval=1

ndrfts=size(trajmat_X,2)

nsnapshot_total = size(trajmat_X,1);

nsnapshot_end=nsnapshot_total;

timevec=1:ninterval:nsnapshot_end;

nsnapshot_sampled = length(timevec);

n=1;

tic

for i=1:ndrfts-1
    for j = i+1:ndrfts           
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

dr=sqrt((X1-X2).^2+(Y1-Y2).^2);

%Calculate alpha, the separation angle, for each pair.
Cosine=X1-X2;
Sine = Y1 - Y2;
clearvars X1 X2 Y1 Y2
Fnormalization= sqrt(Cosine.^2+Sine.^2);
% normalize to unit vectors
Cosine = Cosine./Fnormalization;
Sine = Sine./Fnormalization;
clearvars Fnormalization

Alpha=atan2(Sine,Cosine); 
Alpha(Alpha<0)=Alpha(Alpha<0)+pi;%using the fact that structure 
%functions considered are even

load('Synthetic_distaxis.mat','dist_axis','dist_bin');

nangleaxis=16; %Taking 16 equally-spaced bins in alpha. You can modify this number
%But there is not much point doing it, as explained in appendix A in paper.

nanglewall=nangleaxis+1;%Number of angle bin walls
angle_bin=linspace(0,pi,nanglewall);
angle_axis = 0.5*(angle_bin(1:end-1) + angle_bin(2:end));

%Count how many pairs fall into each bin in (r,alpha)
nanglecount=zeros(length(dist_axis),length(angle_axis));

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
toc


%Save everything except for the really big variables
clearvars Alpha Cosine dr Sine

save('Synthetic_alpha.mat')

