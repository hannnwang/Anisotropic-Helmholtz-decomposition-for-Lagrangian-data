function RHS_int=intRHS_trapz_nobin(xaxis,RHS)
%May be more sophisticated if one wants to interpolate near the zero point
if size(xaxis,1)~=size(RHS,1)
    error('Check if the dimension of axis agrees with that of structure functions!')
end
RHS_int=cumtrapz(xaxis,RHS);
end