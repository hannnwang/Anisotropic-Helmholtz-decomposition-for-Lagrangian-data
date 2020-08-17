clear all
close all

ndrfts=40;%number of floaters at the first quadrant
ntimesnap=500;

N=256;
dk=10^(-5);dl=dk;
kmax=dk*(N/2);
Hs=10^(-8)*dk;
k=linspace(-kmax+dk,kmax,N);%length(k)=2*Nmax
k=k-k(N/2);
k=hwmakesymmetric(k(N/2:end));
k(1:N/2)=-k(1:N/2);
l=k;
x=x_of_k(k);
xr=x-x(N/2);
xr=hwmakesymmetric(xr(N/2:end));
xr(1:N/2)=-xr(1:N/2);
yr=xr;
Lx=2*max(xr);
Ly=Lx;
[K,L]=meshgrid(k,l);
[X,Y]=meshgrid(xr,yr);
Kappa=sqrt(K.^2+L.^2);
dx=xr(2)-xr(1);dy=yr(2)-yr(1);
r=sqrt(X.^2+Y.^2);

dl=dk;

Nbins=200;
%This run assumes x and y are uniformly distributed
for itimesnap=1:ntimesnap
    trajmat_X(itimesnap,:)=(rand(1,ndrfts))*(N-1)+1-N/2;%[-N/2+1,N/2]
    trajmat_Y(itimesnap,:)=(rand(1,ndrfts))*(N-1)+1-N/2;
    trajmat_X(itimesnap,:)=trajmat_X(itimesnap,:)/(N/2)*xr(end);
    trajmat_Y(itimesnap,:)=trajmat_Y(itimesnap,:)/(N/2)*xr(end);

end
dfy=zeros(ntimesnap*ndrfts*(ndrfts-1)/2,1);
dfx=zeros(ntimesnap*ndrfts*(ndrfts-1)/2,1);
icount=0;
for itimesnap=1:(ntimesnap)   
    for i=1:ndrfts
        for j=1:i-1
            icount=icount+1;
            dfy(icount)=(trajmat_Y(itimesnap,i)-trajmat_Y(itimesnap,j));
            dfx(icount)=(trajmat_X(itimesnap,i)-trajmat_X(itimesnap,j));
        end
    end
end
% figure
% scatter(dfx,dfy)
% axis square

figure
hist3([dfx dfy],'Nbins',[Nbins Nbins])
figure
N = hist3([dfx dfy],'Nbins',[Nbins Nbins]);
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(min(dfx),max(dfx),size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(min(dfy),max(dfy),size(N_pcolor,1)); % Rows of N_pcolor

h = pcolor(xl,yl,N_pcolor);
colormap('hot') % Change color scheme 
colorbar % Display colorbar
h.ZData = -max(N_pcolor(:))*ones(size(N_pcolor));
ax = gca;
ax.ZTick(ax.ZTick < 0) = [];
axis square

%This run assumes r^2 and theta are uniformly distributed
rcutoff=xr(end-1);

for itimesnap=1:ntimesnap
    theta_rand=(rand(1,ndrfts))*(2*pi);%Uniformly distributed in [0 2pi]
    r2_rand=rand(1,ndrfts)*(rcutoff^2);
    
    trajmat_X(itimesnap,:)=sqrt(r2_rand).*cos(theta_rand);
    trajmat_Y(itimesnap,:)=sqrt(r2_rand).*sin(theta_rand);
end

dfy=zeros(ntimesnap*ndrfts*(ndrfts-1)/2,1);
dfx=zeros(ntimesnap*ndrfts*(ndrfts-1)/2,1);
icount=0;
for itimesnap=1:(ntimesnap)   
    for i=1:ndrfts
        for j=1:i-1
            icount=icount+1;
            dfy(icount)=(trajmat_Y(itimesnap,i)-trajmat_Y(itimesnap,j));
            dfx(icount)=(trajmat_X(itimesnap,i)-trajmat_X(itimesnap,j));
        end
    end
end
% figure
% scatter(dfx,dfy)
% axis square

figure
hist3([dfx dfy],'Nbins',[Nbins Nbins])
figure
N = hist3([dfx dfy],'Nbins',[Nbins Nbins]);
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(min(dfx),max(dfx),size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(min(dfy),max(dfy),size(N_pcolor,1)); % Rows of N_pcolor

h = pcolor(xl,yl,N_pcolor);
colormap('hot') % Change color scheme 
colorbar % Display colorbar
h.ZData = -max(N_pcolor(:))*ones(size(N_pcolor));
ax = gca;
ax.ZTick(ax.ZTick < 0) = [];
axis square

%This run assumes that r^2 is uniformly distributed and theta is not
%uniformly distributed; it follows a triangle distribution instead.
rcutoff=xr(end-1);
for itimesnap=1:ntimesnap
 rsignal=rand;

    r2_rand=rand(1,ndrfts)*(rcutoff^2);
     %theta_rand=(rand(1,ndrfts))*(2*pi);%Uniformly distributed in [0 2pi]
    pd = makedist('Triangular','a',0,'b',pi/5,'c',pi);
    theta_rand=random(pd,size(r2_rand));
    if rsignal<0.5
    trajmat_X(itimesnap,:)=sqrt(r2_rand).*cos(theta_rand);
    trajmat_Y(itimesnap,:)=sqrt(r2_rand).*sin(theta_rand);
else
    trajmat_X(itimesnap,:)=sqrt(r2_rand).*cos(theta_rand+pi);
    trajmat_Y(itimesnap,:)=sqrt(r2_rand).*sin(theta_rand+pi);
    end
end

dfy=zeros(ntimesnap*ndrfts*(ndrfts-1)/2,1);
dfx=zeros(ntimesnap*ndrfts*(ndrfts-1)/2,1);
icount=0;
for itimesnap=1:(ntimesnap)   
    for i=1:ndrfts
        for j=1:i-1
            icount=icount+1;
            dfy(icount)=(trajmat_Y(itimesnap,i)-trajmat_Y(itimesnap,j));
            dfx(icount)=(trajmat_X(itimesnap,i)-trajmat_X(itimesnap,j));
        end
    end
end
% figure
% scatter(dfx,dfy)
% axis square


figure
hist3([dfx dfy],'Nbins',[Nbins Nbins])

figure
subplot(1,2,1)
fx=reshape(trajmat_X,ntimesnap*ndrfts,1);
fy=reshape(trajmat_Y,ntimesnap*ndrfts,1);

N = hist3([fx fy],'Nbins',[Nbins Nbins]);
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;

xl = linspace(min(fx),max(fx),size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(min(fy),max(fy),size(N_pcolor,1)); % Rows of N_pcolor

h = pcolor(xl,yl,N_pcolor);
colormap(flipud(bone)) % Change color scheme 
colorbar % Display colorbar
h.ZData = -max(N_pcolor(:))*ones(size(N_pcolor));
ax = gca;
ax.ZTick(ax.ZTick < 0) = [];
axis square

xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
axis([xr(1) xr(end) yr(1) yr(end)])

hold on
th = 0:pi/50:2*pi;
xunit = rcutoff * cos(th) ;
yunit = rcutoff * sin(th) ;
plot(xunit, yunit);


subplot(1,2,2)
N = hist3([dfx dfy],'Nbins',[Nbins Nbins]);
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(min(dfx),max(dfx),size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(min(dfy),max(dfy),size(N_pcolor,1)); % Rows of N_pcolor

h = pcolor(xl,yl,N_pcolor);
colormap(flipud(bone)) % Change color scheme 
colorbar % Display colorbar
h.ZData = -max(N_pcolor(:))*ones(size(N_pcolor));
ax = gca;
ax.ZTick(ax.ZTick < 0) = [];
axis square

xlabel('$\Delta x$','Interpreter','latex')
ylabel('$\Delta y$','Interpreter','latex')
axis([xr(1) xr(end) yr(1) yr(end)])


hold on
th = 0:pi/50:2*pi;
xunit = rcutoff * cos(th) ;
yunit = rcutoff * sin(th) ;
plot(xunit, yunit);
