
%20190916:modified from Helm_arbitrarymode_angleweighted.m in GLAD/20190619
%(Change it so that less storage is demanded, but it takes long to run)

%20190620
%Calculate the Helmholtz decomp with this. Include the weighting by angles.

%Can evaluate one mode with this.
%The only variable that needs to be specified is "Moment".
%"Moment" needs to be specified
%before the program.
%e.g,Moment=4: calculate the cos(4*(alpha-theta)) mode


%
% clear all
% close all
% doPlotSetup
tic

mu=0;
% Moment=2
M=Moment/2;
if M<0
    error('Moment needs to be a nonnegative even integer')
end

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

dr=sqrt((X1.*cosd(Y1)-X2.*cosd(Y2)).^2+...
    (Y1-Y2).^2)*Rearth*pi/180;
% dr=sqrt(((X1-X2).*cosd(0.5*(Y1+Y2))).^2+...
%     (Y1-Y2).^2)*Rearth*pi/180;%Dhruv's formula 
Cosine=X1.*cosd(Y1)-X2.*cosd(Y2);
Sine = Y1 - Y2;
clearvars X1 X2 Y1 Y2
Fnormalization= sqrt(Cosine.^2+Sine.^2);
% normalize to unit vectors
Cosine = Cosine./Fnormalization;
Sine = Sine./Fnormalization;
clearvars Fnormalization

%Corrected in 20191105:
% THIS IS WRONG!!!!
% Alpha=acos(Cosine);
% Alpha(Alpha==pi)=0;%There is one Alpha that happens to be pi!
Alpha_pairs=atan2(Sine,Cosine);

Alpha_histo=Alpha_pairs;
Alpha_histo(Alpha_histo<0)=Alpha_histo(Alpha_histo<0)+pi;%Using that the function is even
%just for convenience in debugging: 
Alpha=Alpha_histo;

n=1;

for i=1:ndrfts-1
    for j = i+1:ndrfts
        for it=1:ninterval:ntimesnap_end
            if ~isnan(trajmat_X(it,i)) && ~isnan(trajmat_X(it,j))
                dU(n)=trajmat_U(it,i)-trajmat_U(it,j);
                dV(n)=trajmat_V(it,i)-trajmat_V(it,j);
                
                n = n+1;
            end
        end
    end
end

% convert to longitudnal and transverse structure functions
dl = dU.*Cosine + dV.*Sine;
dt = dV.*Cosine - dU.*Sine;

clearvars dU dV Cosine Sine

dll=dl.^2;
dtt=dt.^2;
dlt=dl.*dt;
clearvars dl dt

if Moment>0
    dll2M_cos=dll.*cos(Moment*Alpha_pairs);
    dll2M_sin=dll.*sin(Moment*Alpha_pairs);
    
    dtt2M_cos=dtt.*cos(Moment*Alpha_pairs);
    dtt2M_sin=dtt.*sin(Moment*Alpha_pairs);
    
    dlt2M_cos=dlt.*cos(Moment*Alpha_pairs);
    dlt2M_sin=dlt.*sin(Moment*Alpha_pairs);
elseif Moment==0
    dll2M_cos=dll;
    dll2M_sin=dll;
    
    dtt2M_cos=dtt;
    dtt2M_sin=dtt;
    
    dlt2M_cos=dlt;
    dlt2M_sin=dlt;
end

clearvars dll dtt dlt

%So now we only have dr,dll2M_mu etc.

%load dist_axis and angle counts from HPC outcomes
load('LASER_all_alpha.mat','dist_axis','dist_bin','nanglecount','dist_angle','nangle');
%(This is generated in Alpha.histo.m)

%structure functions

dll2M_cos_ensemble= zeros(length(dist_axis),1);

dtt2M_cos_ensemble= zeros(length(dist_axis),1);

dlt2M_cos_ensemble= zeros(length(dist_axis),1);

dll2M_sin_ensemble= zeros(length(dist_axis),1);

dtt2M_sin_ensemble= zeros(length(dist_axis),1);

dlt2M_sin_ensemble= zeros(length(dist_axis),1);

for ix=1:length(dist_axis)
    ivalid_pair = find(dr<dist_bin(ix+1) & ...
        dr>=dist_bin(ix));
    npair_atx=length(ivalid_pair);%npair_atx is guaranteed from dist_axis_comp.m to be larger than zero
    
    %Calculate the angle weighting factor
    Alpha_atx=nanglecount(ix,:);
    
    Alpha_weight=1./Alpha_atx;
    
    %Get the ensemble expectation
    dll2M_cos_v=zeros(npair_atx,1);
    dll2M_sin_v=zeros(npair_atx,1);
    dtt2M_cos_v=zeros(npair_atx,1);
    dtt2M_sin_v=zeros(npair_atx,1);
    dlt2M_cos_v=zeros(npair_atx,1);
    dlt2M_sin_v=zeros(npair_atx,1);
    
    weight=zeros(npair_atx,1);
    
    for it=1:npair_atx
        dll2M_cos_v(it)=dll2M_cos(ivalid_pair(it));
        dll2M_sin_v(it)=dll2M_sin(ivalid_pair(it));
        
        dtt2M_cos_v(it)=dtt2M_cos(ivalid_pair(it));
        dtt2M_sin_v(it)=dtt2M_sin(ivalid_pair(it));
        
        dlt2M_cos_v(it)=dlt2M_cos(ivalid_pair(it));
        dlt2M_sin_v(it)=dlt2M_sin(ivalid_pair(it));
        
        alpha_cur=Alpha_histo(ivalid_pair(it));
        
        ialpha=find(dist_angle<=alpha_cur & ...
            dist_angle>alpha_cur-1/nangle*pi);
        
        weight(it)=Alpha_weight(ialpha);
    end
    %structure functions
    dll2M_cos_ensemble(ix)=sum(dll2M_cos_v.*weight)/nangle;
    dll2M_sin_ensemble(ix)=sum(dll2M_sin_v.*weight)/nangle;
    dtt2M_cos_ensemble(ix)=sum(dtt2M_cos_v.*weight)/nangle;
    dtt2M_sin_ensemble(ix)=sum(dtt2M_sin_v.*weight)/nangle;
    dlt2M_cos_ensemble(ix)=sum(dlt2M_cos_v.*weight)/nangle;
    dlt2M_sin_ensemble(ix)=sum(dlt2M_sin_v.*weight)/nangle;

end
toc

%added on 20191110: by definition all the modes at n>0 need to be doubled, but the mode at n=0 shouldn't
if Moment~=0
    dll2M_cos_ensemble=dll2M_cos_ensemble*2;
    dll2M_sin_ensemble=dll2M_sin_ensemble*2;
    dtt2M_cos_ensemble=dtt2M_cos_ensemble*2;
    dtt2M_sin_ensemble=dtt2M_sin_ensemble*2;
    dlt2M_cos_ensemble=dlt2M_cos_ensemble*2;
    dlt2M_sin_ensemble=dlt2M_sin_ensemble*2;
end

    

%Save everything except for the really big variables
clearvars Alpha_hist Alpha_pairs dll2M_cos dll2M_sin dlt2M_cos dlt2M_sin dr dtt2M_cos dtt2M_sin
clearvars trajmat_U trajmat_V trajmat_X trajmat_Y
savetitle=sprintf('LASER_all_ensemble_Moment_%d.mat',Moment);
save(savetitle);
