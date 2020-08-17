%This corresponds to the
%alternative angle weighting approach described in Appendix A in paper.

%Count how many drifters are there in each bin in
%alpha at each r. I could have merged this program into this code, and
%currently they are two separate programs due to laziness and historical reasons...
Alpha_histo_synthetic 

tic

load ('Synthetic_alpha.mat')

% calculate time series of pair separation. Same as in dist_axis_synthetic.m
ninterval=1;

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

dr=sqrt((X1-X2).^2+...
    (Y1-Y2).^2);

%Calculate alpha, the separation angle, for each pair.
Cosine=X1-X2;
Sine = Y1 - Y2;
clearvars X1 X2 Y1 Y2
Fnormalization= sqrt(Cosine.^2+Sine.^2);
% normalize to unit vectors
Cosine = Cosine./Fnormalization;
Sine = Sine./Fnormalization;
clearvars Fnormalization

Alpha_pairs=atan2(Sine,Cosine);

Alpha_histo=Alpha_pairs;
Alpha_histo(Alpha_histo<0)=Alpha_histo(Alpha_histo<0)+pi;%using the fact that structure 
%functions considered are even

%Each bin axis in alpha; this will be called later. 
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

%Calculate velocity differences in u and v
n=1;
for i=1:ndrfts-1
    for j = i+1:ndrfts
        for it=1:ninterval:nsnapshot_end
            if ~isnan(trajmat_X(it,i)) && ~isnan(trajmat_X(it,j))
                dU(n)=trajmat_U(it,i)-trajmat_U(it,j);
                dV(n)=trajmat_V(it,i)-trajmat_V(it,j);   
                n = n+1;
            end
        end
    end
end

% convert to longitudnal and transverse velocity components
dl = dU.*Cosine + dV.*Sine;
dt = dV.*Cosine - dU.*Sine;

clearvars dU dV Cosine Sine

%Calculate structure functions.
%Note that this is intended to be memory-saving and fast -- we are using Matlab matrix
%manipulations as much as we can here.

%This is the zeroth modes:
dll=dl.^2;
dtt=dt.^2;
dlt=dl.*dt;
clearvars dl dt


if Mode>0  %Higher modes. The variable Mode is defined in the main program
%Synthetic_run.m.
    dll_M_cos=dll.*cos(Mode*Alpha_axis);%exactly as what is described in Appendix A.
    dll_M_sin=dll.*sin(Mode*Alpha_axis);
    
    dtt_M_cos=dtt.*cos(Mode*Alpha_axis);
    dtt_M_sin=dtt.*sin(Mode*Alpha_axis);
    
    dlt_M_cos=dlt.*cos(Mode*Alpha_axis);
    dlt_M_sin=dlt.*sin(Mode*Alpha_axis);
elseif Mode==0 %If Mode is 0, then no change is needed.
    dll_M_cos=dll;
    dll_M_sin=dll;
    
    dtt_M_cos=dtt;
    dtt_M_sin=dtt;
    
    dlt_M_cos=dlt;
    dlt_M_sin=dlt;
end

clearvars dll dtt dlt

%So now we only have dr,dll_M_mu etc.

%load dist_axis and angle counts
load('Synthetic_alpha.mat','dist_axis','dist_bin','nanglecount','angle_bin','nangleaxis');
%(This is generated in Alpha_histo_synthetic.m)

%Binning, and take averages to estimate structure functions

dll_M_cos_ensemble= zeros(length(dist_axis),1);

dtt_M_cos_ensemble= zeros(length(dist_axis),1);

dlt_M_cos_ensemble= zeros(length(dist_axis),1);

dll_M_sin_ensemble= zeros(length(dist_axis),1);

dtt_M_sin_ensemble= zeros(length(dist_axis),1);

dlt_M_sin_ensemble= zeros(length(dist_axis),1);

for ix=1:length(dist_axis)
    ivalid_pair = find(dr<dist_bin(ix+1) & ...
        dr>=dist_bin(ix));
    npair_atx=length(ivalid_pair);%npair_atx is guaranteed from dist_axis_comp.m to be larger than zero
    nanglecount_ix=nanglecount(ix,:);
    %Calculate the angle weighting factor
    Alpha_atx=nanglecount(ix,:);
    
    Alpha_weight=1./Alpha_atx;
    
    %Get the ensemble expectation
    dll_M_cos_v=zeros(npair_atx,1);
    dll_M_sin_v=zeros(npair_atx,1);
    dtt_M_cos_v=zeros(npair_atx,1);
    dtt_M_sin_v=zeros(npair_atx,1);
    dlt_M_cos_v=zeros(npair_atx,1);
    dlt_M_sin_v=zeros(npair_atx,1);
    %locate all observations into angle bins
    dll_M_cos_anglebins=zeros(1,nangleaxis);
    dll_M_sin_anglebins=zeros(1,nangleaxis);
    dtt_M_cos_anglebins=zeros(1,nangleaxis);
    dtt_M_sin_anglebins=zeros(1,nangleaxis);
    dlt_M_cos_anglebins=zeros(1,nangleaxis);
    dlt_M_sin_anglebins=zeros(1,nangleaxis);
    
    for it=1:npair_atx
        dll_M_cos_v=dll_M_cos(ivalid_pair(it));
        dll_M_sin_v=dll_M_sin(ivalid_pair(it));
        
        dtt_M_cos_v=dtt_M_cos(ivalid_pair(it));
        dtt_M_sin_v=dtt_M_sin(ivalid_pair(it));
        
        dlt_M_cos_v=dlt_M_cos(ivalid_pair(it));
        dlt_M_sin_v=dlt_M_sin(ivalid_pair(it));
        
        alpha_cur=Alpha_histo(ivalid_pair(it));
        
        if alpha_cur==pi
            ialpha=1;
        else
            ialpha=find(angle_bin<=alpha_cur,1,'last');
        end
        
        weight=Alpha_weight(ialpha);
        
        dll_M_cos_anglebins(ialpha)=dll_M_cos_anglebins(ialpha)+...
            dll_M_cos_v*weight;
        dll_M_sin_anglebins(ialpha)=dll_M_sin_anglebins(ialpha)+...
            dll_M_sin_v*weight;
        
        
        dtt_M_cos_anglebins(ialpha)=dtt_M_cos_anglebins(ialpha)+...
            dtt_M_cos_v*weight;
        dtt_M_sin_anglebins(ialpha)=dtt_M_sin_anglebins(ialpha)+...
            dtt_M_sin_v*weight;
        
        
        dlt_M_cos_anglebins(ialpha)=dlt_M_cos_anglebins(ialpha)+...
            dlt_M_cos_v*weight;
        dlt_M_sin_anglebins(ialpha)=dlt_M_sin_anglebins(ialpha)+...
            dlt_M_sin_v*weight;
    end
    dll_M_cos_anglebins(nanglecount_ix==0)=NaN;
    dll_M_sin_anglebins(nanglecount_ix==0)=NaN;
    dtt_M_cos_anglebins(nanglecount_ix==0)=NaN;
    dtt_M_sin_anglebins(nanglecount_ix==0)=NaN;
    dlt_M_cos_anglebins(nanglecount_ix==0)=NaN;
    dlt_M_sin_anglebins(nanglecount_ix==0)=NaN;
    
    %Interpolate the empty bins from their adjacent bins. This is
    %straightforward and not very important, but it's taking up many lines of the code.

    %Case where the empty bins are not at the start or the end
    if nanglecount_ix(1)>0 && nanglecount_ix(end)>0
        %just fillmising directly
        dll_M_cos_anglebins=fillmissing(dll_M_cos_anglebins,'linear');
        dll_M_sin_anglebins=fillmissing(dll_M_sin_anglebins,'linear');
        dtt_M_cos_anglebins=fillmissing(dtt_M_cos_anglebins,'linear');
        dtt_M_sin_anglebins=fillmissing(dtt_M_sin_anglebins,'linear');
        dlt_M_cos_anglebins=fillmissing(dlt_M_cos_anglebins,'linear');
        dlt_M_sin_anglebins=fillmissing(dlt_M_sin_anglebins,'linear');
    else
        if length(nanglecount(nanglecount>0))<3
            warning('There are more than three empty bins in alpha at this r! Be careful if your interpolation is justified.')
        else
            %copy and paste to create perodic conditions
            dll_M_cos_anglebins_expanded=repmat(dll_M_cos_anglebins,1,3);
            dll_M_sin_anglebins_expanded=repmat(dll_M_sin_anglebins,1,3);
            dtt_M_cos_anglebins_expanded=repmat(dtt_M_cos_anglebins,1,3);
            dtt_M_sin_anglebins_expanded=repmat(dtt_M_sin_anglebins,1,3);
            dlt_M_cos_anglebins_expanded=repmat(dlt_M_cos_anglebins,1,3);
            dlt_M_sin_anglebins_expanded=repmat(dlt_M_sin_anglebins,1,3);
            %get the fillmissing for the expanded
            dll_M_cos_anglebins_expanded=fillmissing(dll_M_cos_anglebins_expanded,'linear');
            dll_M_sin_anglebins_expanded=fillmissing(dll_M_sin_anglebins_expanded,'linear');
            dtt_M_cos_anglebins_expanded=fillmissing(dtt_M_cos_anglebins_expanded,'linear');
            dtt_M_sin_anglebins_expanded=fillmissing(dtt_M_sin_anglebins_expanded,'linear');
            dlt_M_cos_anglebins_expanded=fillmissing(dlt_M_cos_anglebins_expanded,'linear');
            dlt_M_sin_anglebins_expanded=fillmissing(dlt_M_sin_anglebins_expanded,'linear');
            
            %get back
            dll_M_cos_anglebins=dll_M_cos_anglebins_expanded(nangleaxis+1:2*nangleaxis);
            dll_M_sin_anglebins=dll_M_sin_anglebins_expanded(nangleaxis+1:2*nangleaxis);
            dtt_M_cos_anglebins=dtt_M_cos_anglebins_expanded(nangleaxis+1:2*nangleaxis);
            dtt_M_sin_anglebins=dtt_M_sin_anglebins_expanded(nangleaxis+1:2*nangleaxis);
            dlt_M_cos_anglebins=dlt_M_cos_anglebins_expanded(nangleaxis+1:2*nangleaxis);
            dlt_M_sin_anglebins=dlt_M_sin_anglebins_expanded(nangleaxis+1:2*nangleaxis);
            
        end
    end
    
    %structure functions
    dll_M_cos_ensemble(ix)=mean(dll_M_cos_anglebins);
    dll_M_sin_ensemble(ix)=mean(dll_M_sin_anglebins);
    dtt_M_cos_ensemble(ix)=mean(dtt_M_cos_anglebins);
    dtt_M_sin_ensemble(ix)=mean(dtt_M_sin_anglebins);
    dlt_M_cos_ensemble(ix)=mean(dlt_M_cos_anglebins);
    dlt_M_sin_ensemble(ix)=mean(dlt_M_sin_anglebins);
    
end
toc

%Note that the higher modes involve a factor of 2 from the definitions, as
%expressed in eq. (8)-(10) in paper.
if Mode~=0
    dll_M_cos_ensemble=dll_M_cos_ensemble*2;
    dll_M_sin_ensemble=dll_M_sin_ensemble*2;
    dtt_M_cos_ensemble=dtt_M_cos_ensemble*2;
    dtt_M_sin_ensemble=dtt_M_sin_ensemble*2;
    dlt_M_cos_ensemble=dlt_M_cos_ensemble*2;
    dlt_M_sin_ensemble=dlt_M_sin_ensemble*2;
end

%Save everything except for the really big variables
clearvars Alpha_hist Alpha_pairs Alpha_axis dll_M_cos dll_M_sin dlt_M_cos dlt_M_sin dr dtt_M_cos dtt_M_sin
clearvars trajmat_U trajmat_V trajmat_X trajmat_Y
savetitle=sprintf('Synthetic_ensemble_Mode_%d.mat',Mode);
save(savetitle);
