clear all 
close all

dist_axis_LASER

CuCvCuv

clear all
close all

Momentvec=0:2:4;
Momentvec(end+1)=100;

for Moment=Momentvec(1:end)
Helm_arbitrarymode_Connor_LASER
%Helm_arbitrarymode_angleweighted_interpolate_LASER
%Helm_arbitrarymode_angleweighted_LASER
end

%Alpha_histo_LASER %for plotting the variance
dLLTT_plots_threepanel
Helmcompare_AfterHPC