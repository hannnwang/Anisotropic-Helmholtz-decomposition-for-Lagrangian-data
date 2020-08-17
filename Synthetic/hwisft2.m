function [Cpsi_sift] = hwisft2(N,dx,dy,k,l,Cpsit,xq,yq)
%[Cpsi_sift] = hwisft2(N,dx,dy,k,l,Cpsit,xq,yq)
%Slow inverse Fourier transform to evaluate inverse FT of Cpsit at (xq,yq)
%xq, yq are the query points in x and y. 
%Note that the index of Cpsit follows the meshgrid tradidions
%Cpsit is N-by-N, and k, l are 1-by-N.
        exvec=exp(1i*k'*xq);%N-by-1
        Halfway=Cpsit*exvec;
        Cpsi_sift=1/(N^2*dx*dy)*...
            sum(Halfway.*exp(1i*l'*yq));

end
