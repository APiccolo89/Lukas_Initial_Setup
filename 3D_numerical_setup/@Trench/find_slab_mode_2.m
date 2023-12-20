function [obj,Phase,Temp] = find_slab_mode_2(obj,A,Weak_Slab,Phase,Temp,Boundary,theta)
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
% Restrict our research to the point the t ARE WORTH TO LOOK. You do not
% find God in a Casino, but in church.
C      = [xc,0.0-obj.R(2)]; % center of the curvature
dl = obj.L0./obj.nseg;
% Spell out the surfaces
T = obj.Slab_surface.Top; % Top surface

M = obj.Slab_surface.MidS; % MidSurface

B = obj.Slab_surface.Bottom; % Bottom

WZ = obj.Slab_surface.WZ_surf;

% Select better the area of the chosen one
if ~strcmp(Weak_Slab,'Weak')
    ind = x>=C(1) & x(:)<=max(T(:,1)) & z(:)>=min(B(:,2)) & z(:)<=1.0 & y>ya & y<yb;
else
    ind = x>=C(1) & x(:)<=max(WZ(:,1)) & z(:)>=obj.Decoupling_depth(1) & y>ya & y<yb;
end


%Prepare the main vectors
d_slab = x.*nan;

l_slab = x.*nan;

continent = x.*nan;
%Prepare other variables

xs = x(ind==1);

zs = z(ind==1);

ds = d_slab(ind==1);

ls = l_slab(ind==1);

cs = l_slab(ind==1).*nan;
%free memory
x = [];

z = [];

y = [];

it = 1;

l = 0;
% Basically the length of the midplane that reaches the decoupling depth
ind_decoupling_1 = find(M(:,2)<=obj.Decoupling_depth,1);

found  = 0;

CPU_AF = cputime;

while l<obj.L0
    %======================================================================
    % = Marcel suggested to use scatterinterpolant and inpolygon. Now the
    % routine is faster.
    %======================================================================

    ln = l+dl;

    itn = it+1;
    %[Spell out the coordinate of the small element]
    if ~strcmp(Weak_Slab,'Weak')
        Ap = [T(it,1),T(it,2)]; % coordinate point A

        Bp = [T(itn,1),T(itn,2)];% coordinate point B

        Cp = [B(it,1),B(it,2)]; %coordinate point C

        Dp = [B(itn,1),B(itn,2)]; % coordinate D;
    else
        Ap = [WZ(it,1),WZ(it,2)]; % coordinate point A

        Bp = [WZ(itn,1),WZ(itn,2)];% coordinate point B

        Cp = [T(it,1),T(it,2)]; %coordinate point C

        Dp = [T(itn,1),T(itn,2)]; % coordinate D;

    end

    % Create the polygon

    xv = [Ap(1),Bp(1),Dp(1),Cp(1)];

    zv = [Ap(2),Bp(2),Dp(2),Cp(2)];
    % Find the chosen one

    [ind_chosen] = inpoly2([xs,zs],[xv',zv']);

    %Set up the scatter interpolant classes
    F1  = scatteredInterpolant(xv',zv',[l l+dl l+dl l]');

    F2 = scatteredInterpolant(xv',zv',[0 0 -obj.D0 -obj.D0]');

    %Compute the length and d as l with linear interpolation
    ls(ind_chosen==1) = F1(double(xs(ind_chosen==1)),double(zs(ind_chosen==1)));

    ds(ind_chosen==1) = F2(double(xs(ind_chosen==1)),double(zs(ind_chosen==1)));

    if it == ind_decoupling_1 & found ==0  % Theoretically this should be avoided
        l_dc = l;
        found = 1;
    end

    it = itn;
    l = ln;

end
CPU_BF = cputime;

% Compute the poligon of the continental crust
disp(['Slab segments identification took ',num2str(CPU_BF-CPU_AF,2), 'seconds'])

if ~strcmp(Weak_Slab,'Weak')

    if ~strcmp(obj.length_continent{2},'none')

        [L] = Compute_properties_along_function(ya,yb,obj.length_continent,A.Length_along(ind==1));
        Points = obj.Subducted_crust_L;

        Pl(1)     = Points(1);

        Pl(2)     = Points(3);

        Pl(3)     = 1.0;

        Pd(1)     = Points(2);

        Pd(2)     = Points(4);

        Pd(3)     = Points(6);

        [in,~] = inpoly2([ls./L(:),ds],[Pl',Pd']);

        cs(in==1)=1.0;
    else
        Points = obj.Subducted_crust_L;

        Pl(1)     = Points(1);

        Pl(2)     = Points(3);

        Pl(3)     = Points(5);

        Pd(1)     = Points(2);

        Pd(2)     = Points(4);

        Pd(3)     = Points(6);

        [in,~] = inpoly2([ls,ds],[Pl',Pd']);

        cs(in==1)=1.0;


    end
    % Fill up the main portion
    d_slab(ind==1) = ds;

    l_slab(ind==1) = ls;

    continent(ind==1) = cs;

    % Correct Phase and Temperature

    Temp(ind==1 & isnan(d_slab) & A.Zpart(:)<-obj.D0)= obj.Thermal_information.TP;

    Phase(ind==1 & isnan(d_slab) & A.Zpart(:)<-obj.D0) = obj.Thermal_information.Ph_Ast;

    % Find the first particles fullfilling the required dept

    ind_decoupling = find(l_slab<=l_dc,1);

    % Update
    obj.l_slab = reshape(l_slab,size(A.Xpart));

    obj.d_slab = reshape(d_slab,size(A.Xpart));

    obj.continent = reshape(continent,size(A.Xpart));

    obj.Decoupling_depth(1) = obj.Decoupling_depth;

    obj.Decoupling_depth(2) = ind_decoupling;
else

    d_slab(ind==1) = ds;

    obj.d_slab = reshape(d_slab,size(A.Xpart));
end


% Debug variable.
bla = 0;

end
