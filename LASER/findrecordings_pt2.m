
%20190904
%This was run on my own laptop.
clear all
load LASER_P1_pt1.mat

%Sort according to time
%find month, day, hour and minute
date_P1=char(laserP1_ordered(:,2));
month=str2num(date_P1(:,6:7));
day=str2num(date_P1(:,end-1:end));
clearvars date_P1;
time_P1=char(laserP1_ordered(:,3));
hour=str2num(time_P1(:,1:2));
minute=str2num(time_P1(:,4:5));
clearvars time_P1;
for i=1:length(month)
    if month(i)==1
timearray_P1(i)=minute(i)+60*hour(i)+60*24*day(i)+...
    60*24*0;
    elseif month(i)==2
        timearray_P1(i)=minute(i)+60*hour(i)+60*24*day(i)+...
    60*24*31;
    elseif month(i)==3
        timearray_P1(i)=minute(i)+60*hour(i)+60*24*day(i)+...
    60*24*(31+29);
    end
end
timearray_P1=(timearray_P1-min(timearray_P1))/15+1;

laserP1_all(:,1)=str2num(char(laserP1_ordered(:,1)));%floater ID
laserP1_all(:,2)=timearray_P1;%Time array
laserP1_all(:,3)=str2num(char(laserP1_ordered(:,4)));%Lat
laserP1_all(:,4)=str2num(char(laserP1_ordered(:,5)));%lon
laserP1_all(:,5)=str2num(char(laserP1_ordered(:,7)));%U
laserP1_all(:,6)=str2num(char(laserP1_ordered(:,8)));%V

laserP1_all=sortrows(laserP1_all,2);


%Save as a single huge variable, like Dhruv did.
%(using matrix instead of struct array)
trajmat_X=nan*zeros(max(timearray_P1),nfl);
trajmat_Y=trajmat_X;
trajmat_U=trajmat_X;
trajmat_V=trajmat_X;
for i=1:length(day)
    iflt=laserP1_all(i,1);
    it=laserP1_all(i,2);
    trajmat_X(it,iflt)=laserP1_all(i,4);%Note lat and lon
    trajmat_Y(it,iflt)=laserP1_all(i,3);
    trajmat_U(it,iflt)=laserP1_all(i,5);
    trajmat_V(it,iflt)=laserP1_all(i,6);
end

%turns out that it has just ~12 day

