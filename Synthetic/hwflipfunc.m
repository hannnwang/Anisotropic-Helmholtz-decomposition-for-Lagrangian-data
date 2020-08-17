function [Cvu] = hwflipfunc(Cuv)
%Transform Cuv into Cvu, such that Cvu(x,y)=Cuv(-x,-y)
%(only works if the numerical grids are set up as in the main program.)
[Nx,Ny]=size(Cuv);

Cvu=Cuv;

x0=Nx/2;
y0=Ny/2;
for i=1:Nx-1
    for j=1:Ny-1        
        i_reflex=Nx-i;
        j_reflex=Ny-j;
        Cvu(i,j)=Cuv(i_reflex,j_reflex);
    end
end
end

