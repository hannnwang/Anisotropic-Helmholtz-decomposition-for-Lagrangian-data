clear all
close all

%You can decrease ndrfts, nsnapshot and N to run faster. Decreasing ndrfts
%and nsnapshot will increase statistical error, and decreasing N will
%decrease the grid resolution.

ndrfts=400;%number of floaters at the first quadrant
nsnapshot=1000;

%Setting up the grids in real and wavenumber space.
N=512;
dk=2*pi/(100*1000)/2.5;%150 km
kmax=dk*(N/2);
Hs=10^(-8)*dk;
k=linspace(-kmax+dk,kmax,N);%length(k)=2*Nmax
k=k-k(N/2);
k=hwmakesymmetric(k(N/2:end));
k(1:N/2)=-k(1:N/2);
l=k;dl=dk;
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
Alpha=atan2(Y,X);

%Start defining underlying spectra.
%Cpsi
q1=2*pi/(100000);%q1 and q2 mark the wavenumber interval in which the
%synthetic spectra are nonzero.
q2=2*pi/(1000);
npower=-2;

Upsi=0.2;%rms velocity due to psi set to be 0.2

as=5;
bs=1;
m=(1-npower)/2;
Epsi=(as*K.^2+bs*L.^2).^(-m);
Epsi(K.^2+L.^2>=q2^2)=0;
Epsi(K.^2+L.^2<=q1^2)=0;

CS=Upsi^2/(sum(sum(Epsi))*dk*dl/(4*pi));

Epsi=CS*Epsi;

Spsi=Epsi./(K.^2+L.^2)*2;
Spsi(K.^2+L.^2==0)=0;

Cpsi_2D=real(hwifft2(xr,yr,k,l,Spsi));

%Cphi
Uphi=0.1;%rms velocity due to phi set to be 0.2

ah=1;
bh=1;
Ephi=(ah*K.^2+bh*L.^2).^(-m);
Ephi(K.^2+L.^2>=q2^2)=0;
Ephi(K.^2+L.^2<=q1^2)=0;

CH=Uphi^2/(sum(sum(Ephi))*dk*dl/(4*pi));

Ephi=CH*Ephi;

Sphi=Ephi./(K.^2+L.^2)*2;
Sphi(K.^2+L.^2==0)=0;

Cphi_2D=real(hwifft2(xr,yr,k,l,Sphi));

Bphi_2D=2*(Cphi_2D(N/2,N/2)-Cphi_2D);

%Cpsiphi
a=max([as^2+ah^2,bs^2+bh^2,as*bs+ah*bh]);%a is a constant set so that the
%Cauchy-Shwartz inequality (eq. (49) in paper) is met.
%Making sure Spsiphi contains only even components. This is one of the
%assumptions made in paper.
Spsiphi_even=sqrt(CS*CH*(a*(K.^2+L.^2).^2/2).^(-m))./(K.^2+L.^2)*2;
Spsiphi_even(K.^2+L.^2>=q2^2)=0;
Spsiphi_even(K.^2+L.^2<=q1^2)=0;

Spsiphi_odd=0.*K;

Spsiphi=(Spsiphi_even+Spsiphi_odd)*0.9;

Cpsiphi=real(hwifft2(xr,yr,k,l,Spsiphi));

Cphipsi=hwflipfunc(Cpsiphi);
Sphipsi=hwfft2(xr,yr,k,l,Cphipsi);

%Get the "true answers" for different modes of structure functions
Su=L.^2.*Spsi+K.^2.*Sphi-K.*L.*Spsiphi-K.*L.*Sphipsi;
Sv=K.^2.*Spsi+L.^2.*Sphi+K.*L.*Spsiphi+K.*L.*Sphipsi;
Suv=-K.*L.*Spsi+K.*L.*Sphi-L.^2.*Spsiphi+K.^2.*Sphipsi;

Cu=real(hwifft2(xr,yr,k,l,Su));
Cv=real(hwifft2(xr,yr,k,l,Sv));
Cu0=Cu(N/2,N/2);Cv0=Cv(N/2,N/2);

Su_psi=L.^2.*Spsi;Sv_psi=K.^2.*Spsi;
Su_phi=K.^2.*Sphi;Sv_phi=L.^2.*Sphi;
Cu_rr=real(hwifft2(xr,yr,k,l,Su_psi));
Cv_rr=real(hwifft2(xr,yr,k,l,Sv_psi));
Cu_dd=real(hwifft2(xr,yr,k,l,Su_phi));
Cv_dd=real(hwifft2(xr,yr,k,l,Sv_phi));
Cu0_rr=Cu_rr(N/2,N/2);Cv0_rr=Cv_rr(N/2,N/2);
Cu0_dd=Cu_dd(N/2,N/2);Cv0_dd=Cv_dd(N/2,N/2);

Cuv=real(hwifft2(xr,yr,k,l,Suv));
Cuv0=Cuv(N/2,N/2);

Cvu=hwflipfunc(Cuv);

%Drr and Ddd corresponding to the true answer
Drr_test=2*(Cu0_rr-Cu_rr+Cv0_rr-Cv_rr);
Ddd_test=2*(Cu0_dd-Cu_dd+Cv0_dd-Cv_dd);

%In the true answers, modes are calculated from 2D functions generally from the approach
%described in appendix B
Drr_m0_test=hwC0(Drr_test,xr,yr,k,l,K,L,Kappa,N,dk,dl);
Drr_m2_cos_test=hwC2cos(Drr_test,0,xr,yr,k,l,K,L,Kappa,N,dk,dl);
Drr_m2_sin_test=hwC2sin(Drr_test,0,xr,yr,k,l,K,L,Kappa,N,dk,dl);
Ddd_m0_test=hwC0(Ddd_test,xr,yr,k,l,K,L,Kappa,N,dk,dl);
Ddd_m2_cos_test=hwC2cos(Ddd_test,0,xr,yr,k,l,K,L,Kappa,N,dk,dl);
Ddd_m2_sin_test=hwC2sin(Ddd_test,0,xr,yr,k,l,K,L,Kappa,N,dk,dl);

Drr_m4_cos_test=hwC4cos(Drr_test,0,xr,yr,k,l,K,L,Kappa,N,dk,dl);
Drr_m4_sin_test=hwC4sin(Drr_test,0,xr,yr,k,l,K,L,Kappa,N,dk,dl);
Ddd_m4_cos_test=hwC4cos(Ddd_test,0,xr,yr,k,l,K,L,Kappa,N,dk,dl);
Ddd_m4_sin_test=hwC4sin(Ddd_test,0,xr,yr,k,l,K,L,Kappa,N,dk,dl);

dL2_2D=2*(Cu0-Cu).*cos(Alpha).^2+2*(Cv0-Cv).*sin(Alpha).^2+...
    2*(2*Cuv0-Cuv-Cvu).*sin(Alpha).*cos(Alpha);
dT2_2D=2*(Cu0-Cu).*sin(Alpha).^2+2*(Cv0-Cv).*cos(Alpha).^2-...
    2*(2*Cuv0-Cuv-Cvu).*sin(Alpha).*cos(Alpha);
dLT_2D=-2*(Cu0-Cu).*sin(Alpha).*cos(Alpha)+...
    2*(Cv0-Cv).*sin(Alpha).*cos(Alpha)+(2*Cuv0-Cuv-Cvu).*...
    (cos(Alpha).^2-sin(Alpha).^2);

%In the true answers, modes are calculated from 2D functions generally from the approach
%described in appendix B
dL2_m0_test=hwC0(dL2_2D,xr,yr,k,l,K,L,Kappa,N,dk,dl);
dT2_m0_test=hwC0(dT2_2D,xr,yr,k,l,K,L,Kappa,N,dk,dl);
dLT_m0_test=hwC0(dLT_2D,xr,yr,k,l,K,L,Kappa,N,dk,dl);

dL2_m2_cos_test=hwC2cos(dL2_2D,0,xr,yr,k,l,K,L,Kappa,N,dk,dl);
dT2_m2_cos_test=hwC2cos(dT2_2D,0,xr,yr,k,l,K,L,Kappa,N,dk,dl);
dLT_m2_cos_test=hwC2cos(dLT_2D,0,xr,yr,k,l,K,L,Kappa,N,dk,dl);

dL2_m2_sin_test=hwC2sin(dL2_2D,0,xr,yr,k,l,K,L,Kappa,N,dk,dl);
dT2_m2_sin_test=hwC2sin(dT2_2D,0,xr,yr,k,l,K,L,Kappa,N,dk,dl);
dLT_m2_sin_test=hwC2sin(dLT_2D,0,xr,yr,k,l,K,L,Kappa,N,dk,dl);

dL2_m4_cos_test=hwC4cos(dL2_2D,0,xr,yr,k,l,K,L,Kappa,N,dk,dl);
dT2_m4_cos_test=hwC4cos(dT2_2D,0,xr,yr,k,l,K,L,Kappa,N,dk,dl);
dLT_m4_cos_test=hwC4cos(dLT_2D,0,xr,yr,k,l,K,L,Kappa,N,dk,dl);

dL2_m4_sin_test=hwC4sin(dL2_2D,0,xr,yr,k,l,K,L,Kappa,N,dk,dl);
dT2_m4_sin_test=hwC4sin(dT2_2D,0,xr,yr,k,l,K,L,Kappa,N,dk,dl);
dLT_m4_sin_test=hwC4sin(dLT_2D,0,xr,yr,k,l,K,L,Kappa,N,dk,dl);

trajmat_X=NaN*zeros(nsnapshot,ndrfts);trajmat_Y=NaN*zeros(nsnapshot,ndrfts);
trajmat_U=NaN*zeros(nsnapshot,ndrfts);trajmat_V=NaN*zeros(nsnapshot,ndrfts);

rcutoff=xr(end-1)*3/4;

for itimesnap=1:nsnapshot
    % doesn't inflicit reality condition at first.
    
    %Generating samples of psit and phit as described in section 5(a).
    psit=0.*K;
    phit=0.*K;
    m=2;
    Cov=zeros(2,2);%Covariance matrix, as expressed in eq.(50) in paper
    for ik=1:N
        for il=1:N
            if Spsi(ik,il)~=0 && Sphi(ik,il)~=0
                Cov(1,1)=Spsi(ik,il);
                Cov(2,2)=Sphi(ik,il);
                Cov(1,2)=Spsiphi(ik,il);
                Cov(2,1)=Sphipsi(ik,il);
                
                n = 1;
                Chol = chol(Cov);
                temp = Chol'*(randn(m,n)+1i*randn(m,n))/sqrt(2)*sqrt(Lx*Ly);
                psit(ik,il)=temp(1);
                phit(ik,il)=temp(2);
            else %One of them having zero covariance-->uncorrelated
                psit(ik,il)=((randn+1i*(randn))*...
                    sqrt(Spsi(ik,il)*Lx*Ly))/sqrt(2);
                phit(ik,il)=((randn+1i*(randn))*...
                    sqrt(Sphi(ik,il)*Lx*Ly))/sqrt(2);
            end
        end
    end
    ut=-1i*L.*psit+1i*K.*phit;
    vt=+1i*K.*psit+1i*L.*phit;
    %Making up the drifter positions
    rsignal=rand;
    
    r2_rand=rand(1,ndrfts)*(rcutoff^2);
    %theta_rand=(rand(1,ndrfts))*(2*pi);%Uniformly distributed in [0 2pi]
    pd = makedist('Triangular','a',0,'b',pi/5,'c',pi);%Angle is unevenly distributed.
    theta_rand=random(pd,size(r2_rand));
    if rsignal<0.5
        trajmat_X(itimesnap,:)=sqrt(r2_rand).*cos(theta_rand);
        trajmat_Y(itimesnap,:)=sqrt(r2_rand).*sin(theta_rand);
    else
        trajmat_X(itimesnap,:)=sqrt(r2_rand).*cos(theta_rand+pi);
        trajmat_Y(itimesnap,:)=sqrt(r2_rand).*sin(theta_rand+pi);
    end
    
    for i=1:ndrfts
        xq=trajmat_X(itimesnap,i);
        yq=trajmat_Y(itimesnap,i);
        trajmat_U(itimesnap,i)=real(hwisft2(N,dx,dy,k,l,ut,xq,yq))*sqrt(2);
        trajmat_V(itimesnap,i)=real(hwisft2(N,dx,dy,k,l,vt,xq,yq))*sqrt(2);
    end
end
%discard the drifters lying outside rcutoff
trajmat_X(trajmat_X.^2+trajmat_Y.^2>rcutoff^2)=NaN;
trajmat_Y(trajmat_X.^2+trajmat_Y.^2>rcutoff^2)=NaN;
trajmat_U(trajmat_X.^2+trajmat_Y.^2>rcutoff^2)=NaN;
trajmat_V(trajmat_X.^2+trajmat_Y.^2>rcutoff^2)=NaN;

save('Synthetictraj.mat')