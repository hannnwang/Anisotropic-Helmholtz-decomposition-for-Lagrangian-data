%This corresponds to the angle weighting approach described in section 3 in paper.
%Intended to be fast and memory-saving by using built-in Matlab matrix 
%manipulations as much as possible.
load('Synthetic_distaxis.mat');
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
    dll_M_cos=dll.*cos(Mode*Alpha_pairs);
    dll_M_sin=dll.*sin(Mode*Alpha_pairs);
    
    dtt_M_cos=dtt.*cos(Mode*Alpha_pairs);
    dtt_M_sin=dtt.*sin(Mode*Alpha_pairs);
    
    dlt_M_cos=dlt.*cos(Mode*Alpha_pairs);
    dlt_M_sin=dlt.*sin(Mode*Alpha_pairs);
elseif Mode==0
    dll_M_cos=dll;
    dll_M_sin=dll;
    
    dtt_M_cos=dtt;
    dtt_M_sin=dtt;
    
    dlt_M_cos=dlt;
    dlt_M_sin=dlt;
end



clearvars dll dtt dlt

%So now we only have dr,dll_M_mu etc.

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
    
    %Get the ensemble expectation
    dll_M_cos_v=zeros(npair_atx,1);
    dll_M_sin_v=zeros(npair_atx,1);
    dtt_M_cos_v=zeros(npair_atx,1);
    dtt_M_sin_v=zeros(npair_atx,1);
    dlt_M_cos_v=zeros(npair_atx,1);
    dlt_M_sin_v=zeros(npair_atx,1);
    
    alpha_v=zeros(npair_atx,1);
    
    for it=1:npair_atx
        dll_M_cos_v(it)=dll_M_cos(ivalid_pair(it));
        dll_M_sin_v(it)=dll_M_sin(ivalid_pair(it));
        
        dtt_M_cos_v(it)=dtt_M_cos(ivalid_pair(it));
        dtt_M_sin_v(it)=dtt_M_sin(ivalid_pair(it));
        
        dlt_M_cos_v(it)=dlt_M_cos(ivalid_pair(it));
        dlt_M_sin_v(it)=dlt_M_sin(ivalid_pair(it));
        
        alpha_v(it)=Alpha_histo(ivalid_pair(it));
        
    end
    
    %Sort them according to alpha in prepartion for trapezoidal
    %integration.
    Vall=[alpha_v dll_M_cos_v dll_M_sin_v dtt_M_cos_v dtt_M_sin_v ...
        dlt_M_cos_v dlt_M_sin_v];
    
    Vall=sortrows(Vall,1);    %now they are sorted.
    alpha_v=Vall(:,1);
    
    %Merge repeating values in alpha_v, as noted in the last sentence of
    %the third last paragraph in section 3 in paper; this rarely happens in data.
    idelete=0;
    for ii=1:npair_atx
        if ii>idelete(end)
            ngap=1;
            for ij=(ii+1):npair_atx
                if alpha_v(ij)-alpha_v(ii)<10e-8
                    idelete(end+1)=ij;
                    Vall(ii,:)=...
                        Vall(ii,:)+Vall(ij,:);
                    ngap=ngap+1;
                else
                    break
                end
            end
            Vall(ii,:)=Vall(ii,:)/ngap;
        end
    end
    Vall(idelete(2:end),:)=[];
    
    alpha_v=Vall(:,1);
    dll_M_cos_v=Vall(:,2);
    dll_M_sin_v=Vall(:,3);
    dtt_M_cos_v=Vall(:,4);
    dtt_M_sin_v=Vall(:,5);
    dlt_M_cos_v=Vall(:,6);
    dlt_M_sin_v=Vall(:,7);
    
    %structure functions
    dll_M_cos_ensemble(ix)=(trapz(alpha_v,dll_M_cos_v))/pi;
    dll_M_cos_ensemble(ix)=(trapz(alpha_v,dll_M_cos_v))/pi;
    dll_M_sin_ensemble(ix)=(trapz(alpha_v,dll_M_sin_v))/pi;
    dtt_M_cos_ensemble(ix)=(trapz(alpha_v,dtt_M_cos_v))/pi;
    dtt_M_sin_ensemble(ix)=(trapz(alpha_v,dtt_M_sin_v))/pi;
    dlt_M_cos_ensemble(ix)=(trapz(alpha_v,dlt_M_cos_v))/pi;
    dlt_M_sin_ensemble(ix)=(trapz(alpha_v,dlt_M_sin_v))/pi;
    
end

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
clearvars Alpha_histo Alpha Alpha_pairs dll_M_cos dll_M_sin dlt_M_cos dlt_M_sin dr dtt_M_cos dtt_M_sin
clearvars trajmat_U trajmat_V trajmat_X trajmat_Y
clearvars alpha_v dll_M_cos_v dll_M_sin_v dlt_M_cos_v dlt_M_sin_v dtt_M_cos_v dtt_M_sin_v Vall
savetitle=sprintf('Synthetic_ensemble_Mode_%d.mat',Mode);
save(savetitle);
