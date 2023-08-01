function [Size,xc,yc] = modify_terrane_limits(obj,R,arc_length,Boundary)
%=========================================================================%
% Argument Radius of curvature
% R = radius of curvature
% arc_length = the length of the curved slab
%=========================================================================%
theta_ = arc_length./(R);
xa     = -R.*sin(theta_./2);
xb     = R.*sin(theta_./2);
Size = xb-xa;
if strcmp(Boundary,'A')||strcmp(Boundary,'C')
    xc     = -R.*cos(theta_./2);
    yc     = obj.c(2);

else
    xc     = obj.c(1);
    yc     = -R.*cos(theta_./2);
end

end