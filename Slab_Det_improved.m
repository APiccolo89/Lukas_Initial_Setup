
clear all;
close all; 

T.R        = 40;   %curvature radius 
T.theta_c  = 30;   %curvature radius ingested continental crust
T.theta_dc = 20;   % additional curvature to emulate passive margin (optional)
T.theta    = 90;   % curvature slab
T.tk_WZ    = 10;   % thickness of the weak zone
T.L0       = 100;  % length of the slab from the bottom of the lithosphere
T.D0       = 80;   % Thickness of the slab
T.C  = [0.0 -T.D0-T.R]; 
T.r  = [T.R T.R+T.D0]; 
T.r_WZ = [T.r(2), T.r(2)+T.tk_WZ];


phases.Ph_Ar  = [0,10]  ;% Air 
phases.Ph_UC  = [1,2700]  ;% Upper Crust 
phases.Ph_LC  = [2,2800]  ;% Lower Crust 
phases.Ph_Clt = [3,3300] ;% Continental Lithosphere
phases.Ph_Clt2 = [4,3300] ; 
phases.Ph_WZ  = [8,3300]  ;% Weak Zone 
phases.Ph_OLt = [6,3367.957866] ;% Ocean Lithosphere 
phases.Ph_UM  = [5,3300]  ;% Upper Mantle 
phases.Ph_OC  = [7,3367.957866]   ;%place holder
phases.Ph_LC2  = [9,2700]  ;
phases.Ph_UC2  = [10,2800] ; 
phases.Ph_sed_oc = [11,2680];
phases.Ph_cont_pr = [12,2700];

Buffer         = Terrane;
Continent1     = Terrane;
Oceanic_Plate  = Terrane;
Continent2     = Terrane;


% Organize 
%% Buffer Terranes: Terranes at the left most area of the numerical domains: 
%  it is the column of material used as reference for the lithostatic
%  calculation. And, it is a low viscosity terranes that allow convergence
%  velocity; 
%==========================================================================
Buffer.order        = 1; 
Buffer.Type         = 'Buffer';
Buffer.x_lim        = [-1500.0,-1300.0];
Buffer.Phases       = [phases.Ph_sed_oc(1),phases.Ph_OC(1),phases.Ph_OLt(1)]; 
Buffer.Stratigraphy = [0.0,-2.0,-7.0,-80.0];
Buffer.Age          = 5.0;

%% Continental Terranes 1 
%Terranes.Continent1
Continent1.order = 2; 
Continent1.Type  = 'Continent';
Continent1.x_lim = [-1300.0,-600.0];
Continent1.Phases = [phases.Ph_UC(1),phases.Ph_LC(1),phases.Ph_Clt(1)];
Continent1.Stratigraphy = [0.0,-15.0,-30.0,-100.0];
Continent1.Age          = 100.0; 
Continent1.Passive_Margin = {'left','right'};
Continent1.Passive_Margin_phase = phases.Ph_cont_pr;
%% Continent Terranes 2 
Continent2.order = 3; 
Continent2.Type  = 'Continent';
Continent2.x_lim = [0.0,1500];
Continent2.Phases = [phases.Ph_UC2(1),phases.Ph_LC2(1),phases.Ph_Clt2(1)];
Continent2.Stratigraphy = [0.0,-15.0,-30.0,-100.0];
Continent2.Age          = 100.0; 
Continent2.Accretion_prism = 'Prism'; 
Continent2.Trench_properties = T; 
Continent2.prism_phase      = phases.Ph_cont_pr;
%% Oceanic plate 
Oceanic_Plate.order = 4; 
Oceanic_Plate.Type  = 'Ocean';
Oceanic_Plate.x_lim = [-600.0,0.0]; 
Oceanic_Plate.Phases = [phases.Ph_sed_oc(1),phases.Ph_OC(1),phases.Ph_OLt(1)];
Oceanic_Plate.Stratigraphy = [0.0,-2.0,-7.0,-80.0];
Oceanic_Plate.Age          = 80.0; 
Oceanic_Plate.Trench        = 'Subduction'; 
Oceanic_Plate.Trench_properties = T;  

Terranes = struct('Buffer',Buffer,'Continent1',Continent1,'Continent2',Continent2,'Oceanic_Plate',Oceanic_Plate); 
%% Generic information numerical domain: 
Gen.T_P = 1300; 
Gen.T_S = 20; 
Gen.Ph_Air   = phases.Ph_Ar(1);
Gen.Ph_UM    = phases.Ph_UM(1);
Gen.WZ       = phases.Ph_WZ(1); 

Create_Setup(Terranes,phases,Gen);
disp('The Setup is finished')



function Create_Setup(Terranes,ph,Gen)
    addpath matlab   %Here -> Folder2LaMEM/matlab => all the function that handle matlab files are there
    addpath(genpath('../geomio/src'));
    %==========================================================================
    % OUTPUT OPTIONS
    %==========================================================================
    % See model setup in Paraview 1-YES; 0-NO
    Paraview_output        = 1;
    % Output parallel files for LaMEM, using a processor distribution file (msetup = parallel)
    LaMEM_Parallel_output  = 1;
    Parallel_partition     = 'ProcessorPartitioning_8cpu_4.1.2.bin';
    RandomNoise             =   logical(0);
    Is64BIT                 =   logical(0);
    %==========================================================================
    % LOAD MESH GRID FROM LaMEM PARTITIONING FILE
    %==========================================================================
    npart_x =   3;
    npart_y =   3;
    npart_z =   3;

    % Load grid from parallel partitioning file
    [X,Y,Z,x,y,z, Xpart,Ypart,Zpart] = FDSTAGMeshGeneratorMatlab(npart_x,npart_y,npart_z,Parallel_partition, RandomNoise, Is64BIT);
    Gr.x_g = [min(x),max(x)]; 
    Gr.z_g =[min(z),max(z)];
    Gr.y_g = [min(y),max(y)]; 
    % Update variables (size of grid in [x,y,z] direction
    nump_x  =   size(X,2);
    nump_y  =   size(X,1);
    nump_z  =   size(X,3);  
    W       =   max(x)-min(x);
    mW      =   abs(min(x));
    L       =   max(y)-min(y);
    mL      =   abs(min(y));
    H       =   max(z)-min(z);
    mH      =   abs(min(z)); 
    Xvec    =   squeeze(X(1,:,1));
    Yvec    =   squeeze(Y(:,1,1));
    Zvec    =   squeeze(Z(1,1,:));
    % Temporary save memory
    clear Z Xpart Ypart Zpart
    % Prepare data for visualization/output
    A           =   struct('W',[],'L',[],'H',[],'nump_x',[],'nump_y',[],'nump_z',[],'Phase',[],'Temp',[],'x',[],'y',[],'z',[],'npart_x',[],'npart_y',[],'npart_z',[]);
    % Linear vectors containing coords
    [A.Xpart,A.Ypart,A.Zpart] =meshgrid(single(Xvec),single(Yvec),single(Zvec));
    Phase = zeros(size(A.Xpart));
    Temp  = zeros(size(A.Xpart));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                    Initial Input                                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    % Create useful structure to avoid function input argument long and
    % tedius: 
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Set Mantle Phase
    Phase(:,:,:)  = Gen.Ph_Air;
    Temp(:,:,:)   = Gen.T_P;

    % Set up continental terranes: 
    % Temperature is computed only for the layer considered, while the rest
    % is assumed to be not affected by the half space cooling - just to
    % simplify a bit 
    % To Do => Convert this shit in Python
    % To Do => Create ad hoc function
    
    %=====================================================================%
    terranes_list = fieldnames(Terranes); 
    for it =1:length(terranes_list)
            t = Terranes.(terranes_list{it});
            disp(['Filling the ',terranes_list{it} ,' terranes '])
            [Phase,Temp] =  Set_Phase_Temperature(A,Phase,Temp,t,Xvec,Zvec,Gen);
    end
    %===========================================================================%
    % Set Air Phase
    ind = Phase == 0 & A.Zpart<0.0;
    Phase(ind)  = Gen.Ph_UM;
    Temp(ind)   = Gen.T_P;
    ind = Phase == 0 & A.Zpart>0.0;
    Temp(ind)   = Gen.T_S;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    A.Xpart  =  permute(A.Xpart,[2 1 3]);
    A.Ypart  =  permute(A.Ypart,[2 1 3]);
    A.Zpart  =  permute(A.Zpart,[2 1 3]);


    % We can still manually change all this & include, say, a lower mantle
    A.nump_x = nump_x;
    A.nump_y = nump_y;
    A.nump_z = nump_z;
    A.Phase  = double(Phase); clear Phase
    A.Temp   = double(Temp);  clear Temp
    A.Phase  = permute(A.Phase,[2 1 3]);
    A.Temp   = permute(A.Temp, [2 1 3]);


    x = Xvec;
    y = Yvec;   
    z = Zvec;   
    A.x      =  double(x(:));
    A.y      =  double(y(:));
    A.z      =  double(z(:));
    
    [A] = displace_phase_isostasy(ph,A,Gr,Gen); 
   
    A.RandomNoise = RandomNoise;

    clear Temp Phase

    % Clearing up some memory for parallel partitioning
    % clearvars -except A Paraview_output LaMEM_Parallel_output Parallel_partition Is64BIT

    % PARAVIEW VISUALIZATION
    if (Paraview_output == 1)
         FDSTAGWriteMatlab2VTK(A,'BINARY'); % default option
    end

    % SAVE PARALLEL DATA (parallel)
    if (LaMEM_Parallel_output == 1)
        FDSTAGSaveMarkersParallelMatlab(A,Parallel_partition, Is64BIT);
    end

end


%==========================================================================
% Function for constructing the initial setup [Still tuned for 2D and
% monodirectional]
%==========================================================================


function [Phase,Temp] =  Set_Phase_Temperature(A,Phase,Temp,Terranes,Xvec,Zvec,Gen)
type = Terranes.Type; 
trench = Terranes.Trench; 
passive_margin = Terranes.Passive_Margin;
accretion_prism = Terranes.Accretion_prism; 
T_Tk = [0.0, Terranes.Stratigraphy(end)];

% Select the layer and portion of the model:
indx = A.Xpart >= Terranes.x_lim(1) & A.Xpart < Terranes.x_lim(2);
% Compute the thermal profile
[Temp] = compute_temperature_profile(A,Temp,1,T_Tk,Gen,Terranes,indx);
% Set the background phase
T = T_Tk(1);

for i=1:length(Terranes.Phases)
    if (i == length(Terranes.Phases))
        B=T_Tk(2);
    else
        B = Terranes.Stratigraphy(i+1);
    end
    ind = A.Zpart < T & A.Zpart>=B & indx>0 & Phase == 0;
    Phase(ind) = Terranes.Phases(i);
    T=B;
    ind = [];
end

if strcmp(trench,'Subduction')
    Xvec = A.Xpart(1,:,1);
    Zvec = A.Zpart(1,1,:);
    T = Terranes.Trench_properties;
    C     = T.C; 
    D0    = abs(T.D0);
    ind_z = find(Zvec >= C(2)-T.L0 & Zvec<=0.0);
    ind_x = find(Xvec>=C(1) & Xvec <= C(1)+D0*4); 
    [Layout,arc_angleS] = find_slab_(A,ind_x,ind_z,T,'Slab',D0);
    % Set the temperature
    [Temp] = compute_temperature_profile(Layout,Temp,-1,T_Tk,Gen,Terranes,ind_x);%(A,Temp,Type,T_tk,Gen,Terranes,indx)
    % Correct Temperaure 
    %ind = isnan(Layout);
    %Temp(ind) = Gen.T_P; 
    id1 = min(A.Xpart(~isnan(Layout)));
    id2 = max(A.Xpart(~isnan(Layout)));
    id3 = min(A.Zpart(~isnan(Layout)));
    ind_x1 = find(Xvec>=id1,1);
    ind_x2 = find(Xvec>=id2,1);
    ind_z1 = find(Zvec>=id3,1);
    for i= ind_x1:ind_x2
        ind_L = find((Layout(1,i,:)==min(Layout(1,i,:))),1);
        Phase(:,i,ind_z1:ind_L) = Gen.Ph_Air;
        Temp(:,i,ind_z1:ind_L) = Gen.T_P;
    end


    Tk = T_Tk(1);
    for i=1:length(Terranes.Phases)
        if (i == length(Terranes.Phases))
            B=T_Tk(2);
        else
            B = Terranes.Stratigraphy(i+1);
        end
        ind = Layout < Tk & Layout>=B;
        Phase(ind) = Terranes.Phases(i);
        Tk=B;
        ind = [];
    end
% 
%     theta_c = Terranes.theta_c; 
%     r       = Terranes.r; 
%     ph_C    = [Terranes.Cont.Phases(1);Terranes.Cont.Phases(2)];
%     ph_C_D  = [Terranes.Cont.Depths(1);Terranes.Cont.Depths(2)];
%     
%     TOP     = 0.0; 
% 
%     for i=1:length(ph_C)
%         BOTTOM = TOP-ph_C_D(i); 
%         ind = Layout< TOP & Layout >= BOTTOM & arc_angleS<theta_c; 
%         Phase(ind) = ph_C(i); 
%         TOP = BOTTOM; 
%         ind = []; 
%     end
% 
%     P1 = [C(1);-D0];
%     P2 = [C(1)+(r(2)-sum(ph_C_D))*sin(theta_c*pi/180);C(2)+(r(2)-sum(ph_C_D))*cos(theta_c*pi/180)];
%     m  = (P2(2)-P1(2))/(P2(1)-P1(1));
%     z   = m*(A.Xpart-P1(1))+P1(2); 
%     Phase(Layout< -sum(ph_C_D) & (A.Zpart>z))= Terranes.Cont.Phases(3); 

    Layout = []; 
    [Layout,arc_angleS] = find_slab_(A,ind_x,ind_z,T,'Weak',D0);
    ind =(Layout<=0.0 & A.Zpart >= -D0);
    Phase(ind) = Gen.WZ;

end
if strcmp(accretion_prism,'Prism')
    % Place holder
end
if strcmp(passive_margin,'none') > 0
   % Place holder 
end

end



function [Temp]= compute_temperature_profile(A,Temp,Type,T_tk,Gen,Terranes,indx)
k = Terranes.K./(Terranes.Cp.*3300); 
T_Age = Terranes.Age.*Terranes.secMyrsyear;
T_P   = Gen.T_P;
T_S   = Gen.T_S; 

if Type == 1
    ind_z = find(A.Zpart<T_tk(1) & A.Zpart>=T_tk(2) & indx == 1);
    erf_function = erf(A.Zpart(ind_z).*1000/2/(k*T_Age)^0.5);
    Temp(ind_z) = T_S - (T_P-T_S).*erf_function;
else
    ind = A < T_tk(1) & A>=T_tk(2);
    erf_function = erf(A(ind).*1000/2/(k*T_Age)^0.5);
    Temp(ind) = T_S - (T_P-T_S).*erf_function;
    ind = A < T_tk(2); 
    Temp(ind) = T_P;

end

end




function [Layout,arc_angleS] = find_slab_(A,ind_x,ind_z,Slab,type,D0)
% 
%
% Really convoluted, but it works: 
% 1) create a polygon for the slab 
% 2) find the curvature 
% 3) save the relative distannce w.t.r the top
%
%
% find the area belonging to the curvature slab: 
x     = A.Xpart;
z     = A.Zpart; 
theta = Slab.theta;
ang_  = (theta)*pi/180; 
ang_2 = (theta+90)*pi/180;
C     = Slab.C; 
r     = Slab.r; 
L0    = Slab.L0 - Slab.R;
sl    = 1               ;
if strcmp(type,'Weak')
    r(1) = r(2); 
    r(2) = r(2)+Slab.tk_WZ;
     if Slab.theta == 90
        sl = 0;
    end
%=============
   
end
   
 
 p1    = [C(1)+r(1)*sin(ang_),C(2)+r(1)*cos(ang_)];
 p2    = [C(1)+r(2)*sin(ang_),C(2)+r(2)*cos(ang_)];
 p3    = [p1(1)+L0*sin(ang_2),p1(2)+L0*cos(ang_2)];
 p4    = [p2(1)+L0*sin(ang_2),p2(2)+L0*cos(ang_2)];


Px = [p1(1) p2(1) p4(1) p3(1)];
Pz = [p1(2) p2(2) p4(2) p3(2)];

[in,out] = inpolygon(x,z,Px,Pz);
% =========================
%  P2X           P4X
%   r2            r2
%   |             |
%  P1X           P3X
%   r1            r1
%========================
ls = size(x);
arc_angleS = ones(size(x))*nan;
d = ones(size(x))*-1000;
Layout = ones(size(x))*nan;

for i=1:ls(1)
    for j=1:length(ind_x)
        for k=1:length(ind_z)
           
            lx = ind_x(j);
            lz = ind_z(k);
            u = [x(i,lx,lz);z(i,lx,lz)];
            v = [C(1);z(i,lx,lz)];
            a = sqrt(u(1).^2+u(2).^2);
            b = sqrt(v(1).^2+v(2).^2);
            %arc_angleS(i,lx,lz) = asin(dot(v,u)./(a.*b)).*180/pi;
            d(i,lx,lz) = sqrt((u(1)-C(1))^2+(u(2)- C(2))^2); 
            arc_angleS(i,lx,lz) = acos((z(i,lx,lz)-C(2))/d(i,lx,lz)).*180/pi;

            if (d(i,lx,lz)>=r(1) & d(i,lx,lz)<=r(2)) 
                if(arc_angleS(i,lx,lz)>=0.0 & arc_angleS(i,lx,lz)<=theta)
                  
                    Layout(i,lx,lz) = d(i,lx,lz)-r(2); 
                end
            end
            if sl == 1 
                if(in(i,lx,lz)>0)
                    P = [x(i,lx,lz);z(i,lx,lz)]; 
                    Layout(i,lx,lz) = -(find_distance_linear(P,p2,p4));
                end
            end
        end
    end
end


end


function [d] = find_distance_linear(P,p1,p2)

% Ok, wiki saves my day because I was too lazy: 
A  = abs((p2(1)-p1(1))*(p1(2)-P(2))-(p1(1)-P(1))*(p2(2)-p1(2)));
B  = (p2(1)-p1(1))^2+(p2(2)-p1(2))^2;
d  = A/sqrt(B);

end



function [A] = displace_phase_isostasy(ph,A,Gr,Gen) 
% =========================================================================
% Function that compute the relative displacement of the phase if we
% consider them in Isostatic equilibrium. 
%==========================================================================
% Theoretical background: I choose the bottom of the model as reference
% depth. This reference depth represents the compensation depth (i.e.,
% where the lithostatic pressure should be equal and do not show any
% horizontal gradient). Then, using the leftest colum i compute the
% refernce density of this column. The density is computed using a DataBase
% of densities (in the later stage I introduce a simplified
% pressure-temperature dependence on the density) associated to the phase,
% and the average density is assumed to be proportional to the fraction of
% each phase along z direction. 
% Then the topography is computed w.r.t. density reference colum. After
% this simple computation, I use a moving average to smooth the topography
% and use the topography as displacement vector for each marker along z. 
% Then, as last step, I interpolate the marker population and temperature
% column wise using the nearest to the original grid, and fill the the nan
% value emerging with the bottom phase or the air phase. 
%==========================================================================
% Input argument:
% A = Structure of particles 
% ph = Phase data base 
% Gr = numerical grid information
% Output Argument: 
% A = update structure of particles 
%==========================================================================

[topo,rho] = compute_topography_(A,ph);

Ph = squeeze(A.Phase(:,1,:));
T  = squeeze(A.Temp(:,1,:));
Z  = squeeze(A.Zpart(:,1,:));
x  = squeeze(A.Xpart(:,1,:));
ilx = length(x(:,1));
ily = length(squeeze(A.Ypart(1,:,1)))
Ph2 = 0.*Ph;
T2  = 0.*Ph; 
topo_M = movmean(topo,500); 

for i = 1:ilx
    z = Z(i,:)+topo_M(i);
    Ph2(i,:)=interp1(z,Ph(i,:),Z(i,:),'nearest');
    T2(i,:)=interp1(z,T(i,:),Z(i,:),'nearest'); 
    z = [];
end
 
Ph2(isnan(Ph2)& Z >0.0)=ph.Ph_Ar(1);
Ph2(isnan(Ph2)& Z <0.0)=ph.Ph_UM(1);
T2(isnan(Ph2)& Z >0.0)=Gen.T_S; 
T2(isnan(Ph2)& Z <0.0)=Gen.T_P; 

for iy = 1:ily 
    A.Phase(:,iy,:) = Ph2; 
    A.Temp(:,iy,:) = T2; 
end


[s]=save_topography(topo_M,Gr); 
disp(s); 

end

function [topo,rho] = compute_topography_(A,ph)

 x = (squeeze(A.Xpart(:,1,1)));
 rho_x = 0.*x; 
 topo  = 0.*x;

 for i = 1:length(x)
     Ph = A.Phase(i,:,:);
     Ph = Ph(Ph>0);
     l_ph = unique(Ph);

     rho = 0;
     for ip = 1:length(l_ph)

         Per_ph = length(Ph(Ph==l_ph(ip)))/length(Ph);

         rho    =rho+ Per_ph*phase_density(ph,l_ph(ip));

     end
     if i == 1
         rho_0 = rho;
     end
     topo(i) = 1000.*(rho_0-rho)./rho;
     rho_x(i)  = (rho);
     Ph = [];
     l_ph = [];
     rho = [];
 end
end

function [rho] = phase_density(ph,Ph)
field  = fieldnames(ph); 

for is = 1:numel(field)
    A=ph.(field{is});
    if A(1)==Ph
        rho = A(2);
        break
    end
end
end



function [string] = save_topography(topo_M,Gr)
    
    
    dX = (max(Gr.x_g)-min(Gr.x_g))/1500;
    x_t = min(Gr.x_g):dX:max(Gr.x_g);
    topo_t=interp1(A.x,topo_M,x_t,'nearest','extrap');
    
    Topo = zeros(length(x_t),length(y_g));
    Topo(:,1)        = topo_t;
    Topo(:,2)        = topo_t;
    Easting     = x_t;
    Northing    = y_g;
    
    % For orientation
    % Topo(1,1):     SW
    % Topo(1,end):   SE
    % Topo(end,1):   NW
    % Topo(end,end): NE
    
    % compute grid spacing
    dx = (max(Easting) - min(Easting))/(length(Easting)-1);
    dy = (max(Northing) - min(Northing))/(length(Northing)-1);
    
    
    % transpose Topo to be read correctly by LaMEM
    %Topo    = Topo';
    
    % write binary to be read by LaMEM
    % (FILENAME,[# of points in x-dir, # of points in y-dir, x-coord of SW corner, y_coord of SW corner,
    % grid spacing in x-dir, grid spacing in y-dir, TopoMatrix in vector form])
    PetscBinaryWrite('topo.dat', [size(Topo,1); size(Topo,2); min(Easting);min(Northing); dx; dy; Topo(:)]);
    string = 'Isostatic balance finished, and the resulted topography has been saved in Topo.dat';

end 