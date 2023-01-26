
clear all;
close all; 

R        = 40;   %curvature radius 
theta_c  = 30;   %curvature radius ingested continental crust
theta_dc = 15;   % additional curvature to emulate passive margin (optional)
theta    = 90;   % curvature slab
Tk_WZ    = 20;   % thickness of the weak zone
L0       = 300;  % length of the slab from the bottom of the lithosphere
            
T_Age = 60e6;
Create_Setup(R,theta_c,T_Age,Tk_WZ,L0);
disp('The Setup is finished')



function Create_Setup(R,theta_c,T_Age,tk_WZ,L0); 
%Create setup:  Give R, theta, theta_c,tk_WZ,L0 create an initial setup.
%               The parameter defined locally in this function are considered constant
%               -> P1: Fill up the phase space with the background phase
%               and set its temperature equalt to the potential mantle
%               temperature. 
%              







    addpath matlab   %Here -> Folder2LaMEM/matlab => all the function that handle matlab files are there



    addpath(genpath('../geomio/src'));
    %


    %==========================================================================
    % OUTPUT OPTIONS
    %==========================================================================
    % See model setup in Paraview 1-YES; 0-NO
    Paraview_output        = 1;

    % Output parallel files for LaMEM, using a processor distribution file (msetup = parallel)
    LaMEM_Parallel_output  = 1;


    % Parallel partition file
    % for Debugging - BorisProcessorPartitioning_256cpu_8.8.4.bin

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

    % Scale : 
    secyear = 365.25*60*60*24; 

    % Phase list 

    Ph_Ar  = 0  ;% Air 
    Ph_UC  = 1  ;% Upper Crust 
    Ph_LC  = 2  ;% Lower Crust 
    Ph_Clt = 3  ;% Continental Lithosphere
    Ph_WZ  = 9  ;% Weak Zone 
    Ph_OLt = 6  ;% Ocean Lithosphere 
    Ph_UM  = 5  ;% Upper Mantle 
    Ph_OC  = 8   ;%place holder


    Thermal.T_P = 1300; 
    Thermal.T_S = 20; 
    Thermal.k   = 1e-6;
    Thermal.Age = T_Age; 
    % Basic geometry 
    D0     = 80;
    tk_UC  = 15; 
    tk_LC  = 15;    
    tk_OC  =  8; 
    % Slab_position 
    C  = [0.0-D0/2-R -D0-R]; 
    r  = [R R+D0]; 
    r_WZ = [r(2), r(2)+tk_WZ];

    % Create useful structure to avoid function input argument long and
    % tedius: 
    Continents.Phases_ST = [Ph_UC,Ph_LC,Ph_Clt];
    Continents.Depths    = [tk_UC,tk_LC,D0];
    Continents.Type      = Layer; 
    Continents.Thermal   = Thermal; 

    Slab.Phases          = [Ph_OC,Ph_OClt];
    Slab.Cont            = Continents;     %Better to have the information related to the continents as well 
    Slab.Depths          = [tk_OC,D0];
    Slab.Thermal         = Thermal   ; 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Set Mantle Phase
    Phase(:,:,:)  = Ph_UM;
    Temp(:,:,:)   = T_P;

    % Set up continental terranes: 
    % Temperature is computed only for the layer considered, while the rest
    % is assumed to be not affected by the half space cooling - just to
    % simplify a bit 
    % To Do => Convert this shit in Python
    % To Do => Create ad hoc function
    %=====================================================================%
    ind_z = Zvec>-D0 & Zvec<=0.0; 
    Phase(:,:,ind_z) = Ph_Clt; 

    %=====================================================================%
    % Set up temperature 
    erf_function = erf(A.Zpart(:,:,ind_z).*1000/2/(k*T_Age*secyear)^0.5);
    Temp(:,:,ind_z) = T_S - (T_P-T_S).*erf_function;
    ind_z = [] ;
    % Set up the phases pt 2 
    TOP    = 0.0;
    BOTTOM = TOP - tk_UC; 
    ind_z = Zvec < TOP & Zvec >=BOTTOM; 
    Phase(:,:,ind_z) = Ph_UC;
    ind_z = []; 
    % Set up the phases pt 3 
    TOP    = BOTTOM;
    BOTTOM = TOP - tk_LC; 
    ind_z = Zvec < TOP & Zvec >=BOTTOM; 
    Phase(:,:,ind_z) = Ph_LC;

    TOP    = BOTTOM;
    BOTTOM = TOP - (D0-tk_UC-tk_LC); 
    ind_z = Zvec < TOP & Zvec >=BOTTOM; 
    Phase(:,:,ind_z) = Ph_Clt;
    %% Set up the oceanic lithosphere, stalled geometry: 
    % 1st : construct the curvature slab, with the relative 
    ind_z = find(Zvec >= C(2)-L0 & Zvec<=0.0);
    ind_x = find(Xvec>=C(1) & Xvec <= C(1)+D0*4); 

    [Layout,arc_angleS] = find_slab_(A.Xpart,A.Zpart,ind_x,ind_z,C,r,L0,'Slab',R,D0,1);
    
    d_bot_abs = ones(size(Layout))*nan;

    
    TOP = 0.0; 
    BOTTOM = -D0; 
    ind = Layout < TOP & Layout >= BOTTOM ; 
    Phase(ind) = Ph_OLt;


    
    erf_function = erf(Layout(ind).*1000/2/(k*T_Age*secyear)^0.5);
    Temp(ind) = T_S - (T_P-T_S).*erf_function;
    ind = [] 
    ind = Layout < -D0; 
    Phase(ind) = Ph_UM; 
    Temp(ind) = T_P;
    ind = []
    %%%%%%%%%%%%%%%%
    % Fill up continental crust that is trapped in the subduction 
    TOP = 0.0; 
    BOTTOM = TOP-tk_UC; 
    ind = Layout< 0.0 & Layout >= BOTTOM & arc_angleS<theta_c; 
    Phase(ind) = Ph_UC; 
    
    TOP = BOTTOM; 
    BOTTOM = TOP-tk_LC; 
    ind = Layout< TOP & Layout >= BOTTOM & arc_angleS<theta_c; 
    Phase(ind) = Ph_LC; 
   
    
    %d_bot_abs(Layout<0.0) = -(tk_LC+tk_UC)+((tk_LC+tk_UC))/(30)*(arc_angleS(Layout<0.0)-theta_c); 
   
  %  d_bot_abs(arc_angleS>theta_c+30 & ~isnan(d_bot_abs) ) = - 0.0;
   % d_bot_abs(arc_angleS<=theta_c & ~isnan(d_bot_abs)) = -(tk_LC+tk_UC); 
    
    %Phase(Layout<d_bot_abs)=Ph_OLt;

    P1 = [C(1);-D0];
    P2 = [C(1)+(r(2)-tk_UC-tk_LC)*sin(theta_c*pi/180);C(2)+(r(2)-tk_UC-tk_LC)*cos(theta_c*pi/180)];
    m  = (P2(2)-P1(2))/(P2(1)-P1(1));
    z   = m*(A.Xpart-P1(1))+P1(2); 
    Phase(Layout<-tk_LC-tk_UC & (A.Zpart>z))= Ph_Clt; 



    Layout = [];

    




    % 
   [Layout,arc_angleS] = find_slab_(A.Xpart,A.Zpart,ind_x,ind_z,C,r_WZ,0.0,'Weak',R,tk_WZ,0);
   
    ind =(Layout<=0.0 & A.Zpart>=-D0);
    Phase(ind) = Ph_WZ; 



    

    %===========================================================================%
    % Set Air Phase
    ind_z         = find(Zvec>0);
    Phase(:,:,ind_z)  = Ph_Ar;
    Temp(:,:,ind_z)   = 20;

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

function [Layout,arc_angleS] = find_slab_(x,z,ind_x,ind_z,C,r,L0,type,R,D0,sl);
% 
%
% Really convoluted, but it works: 
% 1) create a polygon for the slab 
% 2) find the curvature 
% 3) save the relative distannce w.t.r the top
%
%
% find the area belonging to the curvature slab: 
ang_  = (90)*pi/180; 
ang_2 = (90+90)*pi/180;
L0    = L0 - R;

if strcmp(type,'Slab')
%=============
    p1    = [C(1)+r(1)*sin(ang_),C(2)+r(1)*cos(ang_)];
    p2    = [C(1)+r(2)*sin(ang_),C(2)+r(2)*cos(ang_)];
    p3    = [p1(1)+L0*sin(ang_2),p1(2)+L0*cos(ang_2)];
    p4    = [p2(1)+L0*sin(ang_2),p2(2)+L0*cos(ang_2)];
    % ============= 

else
    p1    = [C(1)+r(1)*sin(ang_),C(2)+r(1)*cos(ang_)];
    p2    = [C(1)+r(2)*sin(ang_),C(2)+r(2)*cos(ang_)];
    p3    = [p1(1)+L0*sin(ang_2),-D0];
    p4    = [p2(1)+L0*sin(ang_2),-D0];

end
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
                if(arc_angleS(i,lx,lz)>=0.0 & arc_angleS(i,lx,lz)<=90)
                  
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

function [Phase,Temperature] =  Set_Phase_Temperature(Phase,Temperature)
end