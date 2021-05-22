%run all the programs, including the generation of synthetic data,
%and the Helmholtz decomposition etc.
clear all
close all

Datagen %generate synthetic data

xy_distribution_synthetic %Histogram of drifter locations
dxdy_histo %Histogram of separation vectors

%Now, start reconstructing different azimuthal Fourier modes of velocity structure
%functions, and conduct the Helmholtz decomposition. Tighten your seat
%belt, it's a long journey!

dist_axis_synthetic %Calculate the binning in r

%Alpha_histo_synthetic %Count how many drifters are there in each bin in
%alpha at each r. This process is only useful if we use the alternative
%angle weighting approach corresponding to Appendix A in paper, or to the
%program "Arbitrarymode_axisweighted_synthetic". You could
%caption this line out if you are not running that in line 25.

CuCvCuv_synthetic %Calculate the average of u^2, v^2 and uv from all drifter
%observations. This is useful in estimating the decorrelation limits, as
%shown in eq.(36)(37) in paper.

for Mode=0:2:4
    %Arbitrarymode_axisweighted_synthetic   %This corresponds to the
    %alternative angle weighting approach described in Appendix A in paper.
    
    Arbitrarymode_angleweighted_trapezoidal_synthetic   %This corresponds to the angle weighting
    %approach described in section 3 in paper.
    
    %Arbitrarymode_noweighting_synthetic  %This corresponds to the
    %traditional unweighted approach, mentioned in section 3 in paper.
end

dLTplot_synthetic_highermodes %plotting the reconstructed azimuthal Fourier modes of
%velocity structure functions, and check if they agree with the true
%answers as well as the decorrelation limits. The figure is not excatly
%formatted as in the paper. The legends etc. are close, but there is one
%major difference: we are presenting the outcomes at the largest r that is
%spanned in the data here, and do NOT cut it off at around r= 84 km, which
%is what's shown in the paper in section 5. As expected, at r>84 km, the reconstructions are pretty bad.

%Conduct the Helmholtz decomposition, and check if the reconstructed
%azimuthal Fourier modes of the divergent or rotational structure functions
%agree with the true answers.
Anisotropic_Helmholtz_decomp
