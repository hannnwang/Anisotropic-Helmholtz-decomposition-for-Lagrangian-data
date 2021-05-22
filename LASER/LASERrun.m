%Run the algorithm at one specific sub-region in the LASER data set. 
%The codes will use "LASER_all.mat", which formats the data downloaded  from the GRIIDC website. The code to generate "LASER_all.mat" is dataformating_all.m. We are assuming that you have already generated "LASER_all.mat" prior to running this code. 

%The code specifying the location and span of the sub-region is circlespec.m. The circlespec.m will be called in dist_axis_LASER.m and so on in this program. If you want to try another sub-region, you could simply modify the values in circlespec.m and re-run this program. Note that circlespec.m is only written for a sub-region that is is circular in x-y space; if you don't want your sub-region to be circular, then you need to a little more work.

clear all 
close all

%The codes to be called here are all very similar to the ones in /Synthetic/. The only changes are about: 1) converting the latitudes and longitudes from raw data to meridional and zonal distances; 2) selecting the observations that fall into the specified subregion.


dist_axis_LASER %Compute the binning in r, so that there are approximately equal numbers of drifter pairs in each bin. The program circlespec.m, which defines the location of the sub-region, is called in this program. 

CuCvCuv %Calculate the average of u^2, v^2 and uv from all drifter
%observations. This is useful in estimating the decorrelation limits, as
%shown in eq.(36)(37) in paper.


Momentvec=0:2:4; %Compute up to n=4 in the azumuthal Fourier expansion. I refer to n as "Moment" in these codes for no good reason... In the codes in /Synthetic/., they are referred to as "Mode".

Momentvec(end+1)=100;

for Moment=Momentvec(1:end)
Arbitrarymode_trapezoidal_LASER %Angle weighting described in section 3 of the paper

%Arbitrarymode_axisweighted_LASER %This corresponds to the
    %alternative angle weighting approach described in Appendix A in paper. You may need to call Alpha_histo_LASER.m first if you are using this code.
end
dLLTT_plots %Plotting the velocity structure functions
Helmcompare_AfterHPC %Conducting the Helmholtz decomposition and plotting the divergent/rotational structure functions
