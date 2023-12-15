function [obj,Phase,Temp] = find_slab_mode_1(obj,A,Weak_Slab,Phase,Temp,Boundary,theta)
% linear part of the slab: should I compute?
sl    = 1    ; % Always YES, NO, only if we are considering weak zone, which has 90 deg theta
% Convert the angle in radians

x     = A.Xpart(:); % x part coordinate
z     = A.Zpart(:); % z part coordinate
y     = A.Ypart(:); % y part coordinate
data_boundary = obj.Boundary.(Boundary);
%Select the boundary = limits of the boundary
[loc1, ~] = find(tril(data_boundary{1} == data_boundary{1}', -1)); % find the main boundary.
if rem(loc1,2)==0
    ya = data_boundary{1}(1);
    yb = data_boundary{1}(3);
    xc = data_boundary{1}(2);
else
    ya = data_boundary{1}(2);
    yb = data_boundary{1}(4);
    xc = data_boundary{1}(3);

end
% Restrict our research to the point that ARE WORTH TO LOOK. You do not
% find God in a Casino, but in church. 
ind = (x>=C(1) & x(:)<=C(1)+2.*obj.L0) & z>=-1.5*obj.L0 & z(:)<=1.0 & y>ya & y<yb;

%Prepare the main vectors
d_slab = x.*Nan; 

l_slab = x.*Nan;
%Prepare other variables 
it = 1; 

l = 0; 

dl = obj.L0./obj.nseg; 
% Spell out the surfaces
T = obj.Slab_surface.Top;

M = obj.Slab_surface.MidS;

B = obj.Slab_surface.Bottom; 

while l<obj.L0 

    CPU_A = cputime;

    ln = l+dl; 

    itn = it+1; 
    %[Spell out the coordinate of the small element]

    Ap = [T(it,1),T(it,2)]; % coordinate point A

    Bp = [T(itn,2),T(itn,2)];% coordinate point B

    Cp = [B(it,1),B(it,2)]; %coordinate point C 

    Dp = [B(itn,1),B(itn,2)]; % coordinate D; 

    Mp1  = [M(it,1),M(it,2)]; % coordinate point midplane 1

    Mp2  = [M(it,1),M(it,2)]; % coordinate point midplane 2

    % Compute slope AB and CD

    s1 = (Bp(2)-Ap(2))./(Bp(1)-Ap(1));

    s2 = (Dp(2)-Cp(2))./(Dp(1)-Cp(1));
    % Compute AC slope 

    s3 = (Cp(2)-Ap(2))./(Cp(1)-Ap(1));

    s4 = (Dp(2)-Bp(2))./(Dp(1)-Bp(1)); 
    % Compute the slope of the midplane

    sm =(Mp2(2)-Mp1(1))./(Mp2(1)-Mp1(2)); 
    % find the chosen point 
    % The chosen one, the favorite marker of the Godly Slab is the one that
    % are within ya-yb, and that feature z> or < of two pairs of straight
    % line. 

    ind_chosen = ind(:) == 1 & z<=linear_equation_two_points(s1,Ap,x) ...
        & z>=linear_equation_two_points(s2,Cp,x) & z<=linear_equation_two_points(s3,Ap,x) ...
        & z>=linear_equation_two_points(s4,Bp,x);

    % Project the chosen point to the midsurface to find the lenght for the
    % chosen particles 
    % 1) find the projection. From the projection, use l+distance along the
    % projection to find the distance l(P) = l+d(Ppr(in)Midsurface) 
    %


    it = itn; 
    l = ln; 
    CPU_B = cputime;
    disp(['Segment took ',num2str(CPU_B-CPU_A,2), 'seconds'])
end


ang_  = (theta)*pi/180; % angle of the slab
ang_2 = (theta+90)*pi/180; %needed for the linear portion of the slab

r      = obj.R; % radius upper and lower surface
r_m    = sum(r)./2; % average radius, to compute the linear integral
C      = [xc,0.0-obj.R(2)]; % center of the curvature

if strcmp(Weak_Slab,'Weak')
    r(1) = r(2);
    r(2) = r(2)+obj.tk_WZ;
    if theta == 90
        sl = 0;
    end
    %=============
end
% select the point that are worth to check if they belong to the slab
ind = (x>=C(1) & x(:)<=C(1)+2.*obj.L0) & z>=-1.5*obj.L0 & z(:)<=1.0 & y>ya & y<yb;
A_time = cputime;
% Vector version
d = x*nan; % d the distance array for the slab
Length = d; % Lenght the integrated length along the slab
arc_angleS = x.*NaN; % array that stores the angular distance
continent = arc_angleS; % array that identify where is the continental portion
% Compute the angular distance w.r.t. the center C
d(ind==1) = sqrt((x(ind==1)-C(1)).^2+(z(ind==1)- C(2)).^2);
% Compute the angle between the vertical vector stemming from the center of
% circumference and the current marker
arc_angleS(ind==1) = acos((z(ind==1)-C(2))./d(ind==1)).*180/pi;
% find the portion that are considered continent
% save the information in continent
% Exclude the points that feature a distance from the center that are less
% or higher than the actual internal and external radius of circumference
% Correct temperature 1st part: Unfourtantely, temperature is a pain to
% correct and I need to introduce ad hoc correction in specific place.
if ~strcmp(Weak_Slab,'Weak')
    Temp(d<=r(1) & x>0 & z>obj.Decoupling_depth) = obj.Thermal_information.TP;
    Phase(d<=r(1) & x>0 & z>obj.Decoupling_depth) = obj.Thermal_information.Ph_Ast;
end
d(d<r(1) | d>r(2)) = NaN;
% Exclude the points that feature an angular distance w.r.t. C that is
% higher or lower than the actual angle {TO DO detecting the bending}
d(arc_angleS>theta | arc_angleS<0)=NaN;
% Compute the effective distance from the top surface
d = d-r(2);
% Compute the length of the curved slab as a function of the mid distance
% layer. {r_m*angle} => angle in radians.
Length(~isnan(d))= r_m.*arc_angleS(~isnan(d) & arc_angleS<=theta)*pi/180; % compute the length associated with the

if sl ==1

    % Compute the point using basic trigonometry:
    % 1) => Compute the P1,P2,P3,P4 that represents the rectangular shape
    % attached to the curved portion of the slab
    p1(1)    = C(1)+r(1)*sin(ang_);
    p1(2)    = C(2)+r(1)*cos(ang_);
    p2(1)    = C(1)+r(2)*sin(ang_);
    p2(2)    = C(2)+r(2)*cos(ang_);

    p3(1)    = p1(1)+obj.L0*sin(ang_2);
    p3(2)    = p1(2)+obj.L0*cos(ang_2);
    p4(1)    = p2(1)+obj.L0*sin(ang_2);
    p4(2)    = p2(2)+obj.L0*cos(ang_2);

    %
    Py(1) = p1(1);
    Py(2) = p2(1);
    Py(3) = p4(1);
    Py(4) = p3(1);
    Pz(1) = p1(2);
    Pz(2) = p2(2);
    Pz(3) = p4(2);
    Pz(4) = p3(2);
    % Inpolygon => Select the pint that belongs to the polygon
    [in,~] = inpolygon(x,z,Py,Pz);
    P = [x(in & y>=ya & y<=yb), z(in & y>=ya & y<=yb)];
    % Find the distance from the top surface
    d(in & y>=ya & y<=yb) = -(find_distance_linear(P,p2,p4))';
    %Update the object
    % Find the continent portion {additional change to the structure}
    % I have the length and the depth, which mean that I have an ortogonal
    % system of coordinate that I can use to generate my stuff.
    if ~strcmp(Weak_Slab,'Weak')
        [Length,ind_decoupling,Phase,Temp] = find_length_slab(obj,x,z,C,r_m,d,Length,Phase,Temp,theta);
        if ~strcmp(obj.length_continent{2},'none')
            [L] = Compute_properties_along_function(ya,yb,obj.length_continent,A.Length_along);
            Points = obj.Subducted_crust_L;
            Pl(1)     = Points(1);
            Pl(2)     = Points(3);
            Pl(3)     = 1.0;
            Pd(1)     = Points(2);
            Pd(2)     = Points(4);
            Pd(3)     = Points(6);
            [in,~] = inpolygon(Length./L(:),d,Pl,Pd);
            continent(in==1 & y>=ya & y<=yb)=1.0;
        else
            Points = obj.Subducted_crust_L;
            Pl(1)     = Points(1);
            Pl(2)     = Points(3);
            Pl(3)     = Points(5);
            Pd(1)     = Points(2);
            Pd(2)     = Points(4);
            Pd(3)     = Points(6);
            [in,~] = inpolygon(Length,d,Pl,Pd);
            continent(in==1 & y>=ya & y<=yb)=1.0;
        end
        % Need to do as follow, I lose any hope to have any fixed
        % organisation
        obj.Decoupling_depth(1) = obj.Decoupling_depth;
        obj.Decoupling_depth(2) = ind_decoupling;
    else
        d(z<=obj.Decoupling_depth(1))=nan;
    end
end
obj.l_slab = Length;
obj.d_slab = reshape(d,size(A.Xpart));
obj.continent = reshape(continent,size(A.Xpart));

B_time = cputime;
disp(['Finding the slab took = ', num2str((B_time-A_time),3), ' sec']);
end


function [d] = find_distance_linear(P,p1,p2)

% Ok, wiki saves my day because I was too lazy:
A  = abs((p2(1)-p1(1)).*(p1(2)-P(:,2))-(p1(1)-P(:,1)).*(p2(2)-p1(2)));
B  = (p2(1)-p1(1)).^2+(p2(2)-p1(2)).^2;
d  = A./sqrt(B);

end

function [l,ind_decoupling,Phase,Temp] = find_length_slab(obj,x,z,C,r,d,l,Phase,Temp,theta)
slab = ~isnan(d)==1;
% Compute the two point along the midplane of the slab
% 1) After the curved part of the slab
PA(1)=C(1)+r*sin(theta*pi/180);
PA(2)=C(2)+r*cos(theta*pi/180);
%Coordinate of the midsurface at the depth of slab  (bottom of the slab)
PB(1)=PA(1)+obj.L0*cos(theta*pi/180);
PB(2)=PA(2)-obj.L0*sin(theta*pi/180);
% Compute the slope of the line:
m=(PB(2)-PA(2))./(PB(1)-PA(1));
%m = 1./m;
% project all the x points into the line (from:
% https://math.stackexchange.com/questions/62633/orthogonal-projection-of-a-point-onto-a-line)
% Project all the points that belong to the slab into the mid surface
xprojection = d.*0.0;
zprojection = d.*0.0;
b = PA(2)-PA(1).*m;
xprojection(slab==1 & isnan(l)) = (x(slab==1 & isnan(l)) + m.*z(slab==1 & isnan(l))-m.*b)./(1+m.^2);
zprojection(slab==1 &isnan(l)) = (m.*x(slab==1 & isnan(l)) + m^2.*z(slab==1 & isnan(l))+b)./(1+m.^2);
% Compute the length of the slab along the midsurface
l(slab==1 & isnan(l)) = r.*theta*pi/180+sqrt((xprojection(slab==1 & isnan(l))-PA(1)).^2+(zprojection(slab==1 & isnan(l))-PA(2)).^2);
ind_decoupling = find(zprojection>=obj.Decoupling_depth,1);
% Corret the damn Phase and Temperature
Phase(z<PA(2)+m*(x-PA(1)) & isnan(d) & z<PA(2)& x>C(1) & ~isnan(Phase(:))) = obj.Thermal_information.Ph_Ast;
Temp(z<PA(2)+m*(x-PA(1)) & isnan(d) & z<PA(2) & x>C(1) &~isnan(Phase(:))) = obj.Thermal_information.TP;
end

function [z] = linear_equation_two_points(s,A,x)
   
   z = s.*(x-A(1))+A(2);

end

function [pr] = linear_projection(s,M,P)
end