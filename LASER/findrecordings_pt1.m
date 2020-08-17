%20190904
%Find the recordings of LASER that are from a specific launch
%This was run at the school laptop;  the changem function only works for
%older Matlab versions.
clear all
close all

load ID_P1.mat
load laserspotdrifterscleanv15_all.mat

nfl=length(ID_P1);

laserspotdrifterscleanv15(:,1)=erase(laserspotdrifterscleanv15(:,1),"L_");
ID_all=laserspotdrifterscleanv15(:,1);

%Find the ones from P1 launch only
ID_ifP1=ismember(ID_all,ID_P1);

laserP1=laserspotdrifterscleanv15(ID_ifP1,:);

%Simplify the IDs
IDsimple=ID_P1;
for i=1:length(ID_P1)
    IDsimple(i)=num2str(i);   
end
laserP1_ordered=changem(laserP1,IDsimple,sort(ID_P1));
% 
% 
% %Sort according to time
% %sort from date first:
% laserP1=sortrows(laserP1,2);
% %then, sort from time (tricker)
% for 

saveas('LASER_P1_pt1.mat')