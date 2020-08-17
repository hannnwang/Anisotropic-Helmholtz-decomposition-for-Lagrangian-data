%20200120
%Modified from Arbitrarymode_axisweighted_interpolate_synthetic.m
tic

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

% %Only keep the drifters whose positions fall into a circle
	% circlespec
	
	Rearth=6.3782*10^6;
	% 
	% converted_X=trajmat_X.*cosd(trajmat_Y)*Rearth*pi/180;
	% converted_Y=trajmat_Y*Rearth*pi/180;
	% 
	% Radius_test=sqrt((converted_X-center_x).^2+(converted_Y-center_y).^2);
	% 
	% trajmat_X(Radius_test>radius)=NaN;
	% trajmat_Y(Radius_test>radius)=NaN;
	% trajmat_U(Radius_test>radius)=NaN;
	% trajmat_V(Radius_test>radius)=NaN;
	
    
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
dr=sqrt(((X1-X2).*(cosd(Y1)+cosd(Y2))/2).^2+...
    (Y1-Y2).^2)*Rearth*pi/180;
%dr=sqrt((cosave*(X1-X2)).^2+...
	 %   (Y1-Y2).^2)*Rearth*pi/180;
    
% 	dr=sqrt((X1.*cosd(Y1)-X2.*cosd(Y2)).^2+...
% 	    (Y1-Y2).^2)*Rearth*pi/180;
	% dr=sqrt(((X1-X2).*cosd(0.5*(Y1+Y2))).^2+...
	%     (Y1-Y2).^2)*Rearth*pi/180;%Dhruv's formula -- I think it's less accurate

Cosine=(X1-X2).*(cosd(Y1)+cosd(Y2))/2;
% 	Cosine=X1.*cosd(Y1)-X2.*cosd(Y2);
    
Sine = Y1 - Y2;
clearvars X1 X2 Y1 Y2
Fnormalization= sqrt(Cosine.^2+Sine.^2);
% normalize to unit vectors
Cosine = Cosine./Fnormalization;
Sine = Sine./Fnormalization;
clearvars Fnormalization

% THIS IS WRONG!!!!
% Alpha=acos(Cosine);
% Alpha(Alpha==pi)=0;%There is one Alpha that happens to be pi!
Alpha_pairs=atan2(Sine,Cosine);

Alpha_histo=Alpha_pairs;
Alpha_histo(Alpha_histo<0)=Alpha_histo(Alpha_histo<0)+pi;%Using that the function is even

%Alpha that correspond to bin axis:
Alpha_axis=Alpha_histo*NaN;
for i=1:length(Alpha_histo)
    alpha_cur=Alpha_histo(i);
    if alpha_cur==pi
        ialpha=1;
    else
        ialpha=find(angle_bin<=alpha_cur,1,'last');
    end
    
    Alpha_axis(i)=angle_axis(ialpha);
end

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
    dll2M_cos=dll.*cos(Moment*Alpha_axis);
    dll2M_sin=dll.*sin(Moment*Alpha_axis);
    
    dtt2M_cos=dtt.*cos(Moment*Alpha_axis);
    dtt2M_sin=dtt.*sin(Moment*Alpha_axis);
    
    dlt2M_cos=dlt.*cos(Moment*Alpha_axis);
    dlt2M_sin=dlt.*sin(Moment*Alpha_axis);
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
load('LASER_all_alpha.mat','dist_axis','dist_bin','nanglecount','angle_bin','nangleaxis');
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
    nanglecount_ix=nanglecount(ix,:);
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
    
    
    %20191217 locate all observations into angle bins
    
    dll2M_cos_anglebins=zeros(1,nangleaxis);
    dll2M_sin_anglebins=zeros(1,nangleaxis);
    dtt2M_cos_anglebins=zeros(1,nangleaxis);
    dtt2M_sin_anglebins=zeros(1,nangleaxis);
    dlt2M_cos_anglebins=zeros(1,nangleaxis);
    dlt2M_sin_anglebins=zeros(1,nangleaxis);
    
    for it=1:npair_atx
        dll2M_cos_v=dll2M_cos(ivalid_pair(it));
        dll2M_sin_v=dll2M_sin(ivalid_pair(it));
        
        dtt2M_cos_v=dtt2M_cos(ivalid_pair(it));
        dtt2M_sin_v=dtt2M_sin(ivalid_pair(it));
        
        dlt2M_cos_v=dlt2M_cos(ivalid_pair(it));
        dlt2M_sin_v=dlt2M_sin(ivalid_pair(it));
        
        alpha_cur=Alpha_histo(ivalid_pair(it));
        
        if alpha_cur==pi
            ialpha=1;
        else
            ialpha=find(angle_bin<=alpha_cur,1,'last');
        end
        
        weight=Alpha_weight(ialpha);
        
        dll2M_cos_anglebins(ialpha)=dll2M_cos_anglebins(ialpha)+...
            dll2M_cos_v*weight;
        dll2M_sin_anglebins(ialpha)=dll2M_sin_anglebins(ialpha)+...
            dll2M_sin_v*weight;
        
        
        dtt2M_cos_anglebins(ialpha)=dtt2M_cos_anglebins(ialpha)+...
            dtt2M_cos_v*weight;
        dtt2M_sin_anglebins(ialpha)=dtt2M_sin_anglebins(ialpha)+...
            dtt2M_sin_v*weight;
        
        
        dlt2M_cos_anglebins(ialpha)=dlt2M_cos_anglebins(ialpha)+...
            dlt2M_cos_v*weight;
        dlt2M_sin_anglebins(ialpha)=dlt2M_sin_anglebins(ialpha)+...
            dlt2M_sin_v*weight;
    end
    dll2M_cos_anglebins(nanglecount_ix==0)=NaN;
    dll2M_sin_anglebins(nanglecount_ix==0)=NaN;
    dtt2M_cos_anglebins(nanglecount_ix==0)=NaN;
    dtt2M_sin_anglebins(nanglecount_ix==0)=NaN;
    dlt2M_cos_anglebins(nanglecount_ix==0)=NaN;
    dlt2M_sin_anglebins(nanglecount_ix==0)=NaN;
    
    %Interpolate the empty bins
    %Case where the empty bins are not at the start or the end
    if nanglecount_ix(1)>0 && nanglecount_ix(end)>0
        %just fillmising directly
        dll2M_cos_anglebins=fillmissing(dll2M_cos_anglebins,'linear');
        dll2M_sin_anglebins=fillmissing(dll2M_sin_anglebins,'linear');
        dtt2M_cos_anglebins=fillmissing(dtt2M_cos_anglebins,'linear');
        dtt2M_sin_anglebins=fillmissing(dtt2M_sin_anglebins,'linear');
        dlt2M_cos_anglebins=fillmissing(dlt2M_cos_anglebins,'linear');
        dlt2M_sin_anglebins=fillmissing(dlt2M_sin_anglebins,'linear');
    else
        if length(nanglecount(nanglecount>0))<3
            warning('SAMPLES ARE TOO ANISOTROPIC AT THIS R!')
        else
            %copy and paste to create perodic conditions
            dll2M_cos_anglebins_expanded=repmat(dll2M_cos_anglebins,1,3);
            dll2M_sin_anglebins_expanded=repmat(dll2M_sin_anglebins,1,3);
            dtt2M_cos_anglebins_expanded=repmat(dtt2M_cos_anglebins,1,3);
            dtt2M_sin_anglebins_expanded=repmat(dtt2M_sin_anglebins,1,3);
            dlt2M_cos_anglebins_expanded=repmat(dlt2M_cos_anglebins,1,3);
            dlt2M_sin_anglebins_expanded=repmat(dlt2M_sin_anglebins,1,3);
            %get the fillmissing for the expanded
            dll2M_cos_anglebins_expanded=fillmissing(dll2M_cos_anglebins_expanded,'linear');
            dll2M_sin_anglebins_expanded=fillmissing(dll2M_sin_anglebins_expanded,'linear');
            dtt2M_cos_anglebins_expanded=fillmissing(dtt2M_cos_anglebins_expanded,'linear');
            dtt2M_sin_anglebins_expanded=fillmissing(dtt2M_sin_anglebins_expanded,'linear');
            dlt2M_cos_anglebins_expanded=fillmissing(dlt2M_cos_anglebins_expanded,'linear');
            dlt2M_sin_anglebins_expanded=fillmissing(dlt2M_sin_anglebins_expanded,'linear');
            
            %get back
            dll2M_cos_anglebins=dll2M_cos_anglebins_expanded(nangleaxis+1:2*nangleaxis);
            dll2M_sin_anglebins=dll2M_sin_anglebins_expanded(nangleaxis+1:2*nangleaxis);
            dtt2M_cos_anglebins=dtt2M_cos_anglebins_expanded(nangleaxis+1:2*nangleaxis);
            dtt2M_sin_anglebins=dtt2M_sin_anglebins_expanded(nangleaxis+1:2*nangleaxis);
            dlt2M_cos_anglebins=dlt2M_cos_anglebins_expanded(nangleaxis+1:2*nangleaxis);
            dlt2M_sin_anglebins=dlt2M_sin_anglebins_expanded(nangleaxis+1:2*nangleaxis);
            
        end
    end
    
    %structure functions
    dll2M_cos_ensemble(ix)=mean(dll2M_cos_anglebins);
    dll2M_sin_ensemble(ix)=mean(dll2M_sin_anglebins);
    dtt2M_cos_ensemble(ix)=mean(dtt2M_cos_anglebins);
    dtt2M_sin_ensemble(ix)=mean(dtt2M_sin_anglebins);
    dlt2M_cos_ensemble(ix)=mean(dlt2M_cos_anglebins);
    dlt2M_sin_ensemble(ix)=mean(dlt2M_sin_anglebins);
    
end
toc

%added on 20191110 to correct a mistake before. Should therefore delete the
%same steps in the plotting/Helm decomp. codes.
if Moment~=0
    dll2M_cos_ensemble=dll2M_cos_ensemble*2;
    dll2M_sin_ensemble=dll2M_sin_ensemble*2;
    dtt2M_cos_ensemble=dtt2M_cos_ensemble*2;
    dtt2M_sin_ensemble=dtt2M_sin_ensemble*2;
    dlt2M_cos_ensemble=dlt2M_cos_ensemble*2;
    dlt2M_sin_ensemble=dlt2M_sin_ensemble*2;
end

%Save everything except for the really big variables
clearvars Alpha_hist Alpha_pairs Alpha_axis dll2M_cos dll2M_sin dlt2M_cos dlt2M_sin dr dtt2M_cos dtt2M_sin
clearvars trajmat_U trajmat_V trajmat_X trajmat_Y
savetitle=sprintf('LASER_all_ensemble_Moment_%d.mat',Moment);
save(savetitle);
