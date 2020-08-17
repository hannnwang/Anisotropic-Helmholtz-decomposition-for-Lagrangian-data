clear all
%close all
load LASER_2Dstrucfun.mat

ixvec=1;
ixvec(end+1)=10;
ixvec(end+1)=100;
ixvec(end+1)=1000;
    figure
for ifig=1:length(ixvec)
    subplot(2,2,ifig)
    plot(1:nangle,dtt_2D(ixvec(ifig),:))
titname=sprintf('at r=%.1f',dist_axis(ixvec(ifig)));
title(titname)
xlabel('degrees')
ylabel('m^2/s^2')
hold on
end


ax = gca;
ax.FontSize = 13; 
