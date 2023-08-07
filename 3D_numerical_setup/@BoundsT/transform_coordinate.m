function [B,dx] = transform_coordinate(obj,B,Boundary)
%==========================================================================
% function to transform coordinate such that the curved boundary is
% deflected back to linear.
% x - dx, dx = f(x)-xn.
% Approach: take the coordinate (A.Xpart) => compute the deflection imposed
% by the type of the boundary and correct this coordinate for filling up
% the phase.
%==========================================================================
data_boundary = obj.(Boundary);
if strcmp(data_boundary{2},'Circular')
    dx  = obj.circumference_margin(Boundary,B.Ypart(:));
end
[loc1, ~] = find(tril(data_boundary{1} == data_boundary{1}', -1)); % find the main boundary.
% The main boundary is the coordinate that has the same value in
% data_boundary. 
dx  = dx-data_boundary{1}(loc1);
if rem(loc1,2)==0
    ya = data_boundary{1}(1);
    yb = data_boundary{1}(3);
else
    ya = data_boundary{1}(2);
    yb = data_boundary{1}(4);
end
if strcmp(Boundary,'A')||strcmp(Boundary,'D') % in case of A and B, you need to increase the coordinate
    dx = +dx;
else
    dx = -dx;
end
B.Xpart(B.Ypart>=ya & B.Ypart<=yb) = B.Xpart(B.Ypart>=ya & B.Ypart<=yb)+dx(B.Ypart>=ya & B.Ypart<=yb);

end