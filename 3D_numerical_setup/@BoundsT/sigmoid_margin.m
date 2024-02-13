function [y] = sigmoid_margin(obj,Boundary,x)
%
% For being clear: the curvature IS ALWAYS directed towards the external
% portion of the terrain NOT towards its interiour. 
%
%
%
% I do not have clever idea on how to do in a simplified and elegant way. 
if strcmp(Boundary,'A')
elseif strcmp(Boundary,'B')
%     B=obj.B{3}; %function handle boundary
%     R   = B(3);
%     coord  = obj.B{1};
%     xa = coord(1);
%     ya = coord(4);
%     xb = coord(3);
%     c = B(1); 
%     sign = +1; 
elseif strcmp(Boundary,'C')
elseif strcmp(Boundary,'D')
    D = obj.D{3};
    coord  = obj.D{1};
    xa = coord(1);
    ya = coord(4);
    xb = coord(3);%[x0,K1,K2,B,v]
    K1 = D(2)-ya;
    K2 = D(3)-ya;
    x0 = D(1);
    b = D(4);
    v = D(5);

end
y  = K1+(K2-K1)./(1+exp(-b*(x-x0))).^(1/v);
y(x>=xa & x<=xb) = nan; 
end

