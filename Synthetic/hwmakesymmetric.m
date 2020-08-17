function [ Kpsi_p ] = hwmakesymmetric(temp)
%20171223
%Kpsi_p=hwmakesymmetric(temp)
%make a symmeteric 1D function at length Nx from its values at [Nx/2,Nx]
%The input needs to be at size (1,Nx/2+1), and could contain NaN values.
%The output will be at size (1,Nx); the NaN values will be subsituted with
%zero in the output, except when the NaN values are near the zero point
%(i.e. next to the Nx/2th entry)

Mmax=(length(temp)-1)*2;
% temp=hwinterp1(strato_kk_BOTH,Kpsi_strato,k((Mmax/2):(Mmax)));
indices=find(~isnan(temp));
istart=indices(1);
ltemp=length(indices);
Kpsi_p(1:Mmax)=0;
Kpsi_p((Mmax/2+istart-1):(Mmax/2+istart-1+ltemp-1))=temp(~isnan(temp));
if (Mmax/2-(istart-1)-ltemp+1)==0
    temp2=temp(~isnan(temp));
    Kpsi_p(1:(Mmax/2-(istart-1)))=flip(temp2(1:end-1));
else
    Kpsi_p((Mmax/2-(istart-1)-ltemp+1):(Mmax/2-(istart-1)))=flip(temp(~isnan(temp)));
end
Kpsi_p((Mmax/2-(istart-1)+1):(Mmax/2+istart-1-1))=NaN;
Kpsi_p(Mmax/2)=0;%The values at the center won't affect either of the decompositions,
% thus this can be set at an arbitary value; one should expect that now
% Kpsi_p doesn't have NaN values any more
if(mean(isnan(Kpsi_p))>0)
    disp('WARNING: NaN values existent in the symmetic function, even after the midpoint is set zero!')
end

end

