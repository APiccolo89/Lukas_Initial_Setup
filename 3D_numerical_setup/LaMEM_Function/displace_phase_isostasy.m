
function [A,surf] = displace_phase_isostasy(ph,A,Gr,TI)
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
[rho_ij] = Compute_density(A,ph);
[topo,~] = compute_topography_(A,rho_ij);
%==========================================================================
% Compute the length scale of the moving average using a slab of 90 km as
% reference. 
%==========================================================================
x  = squeeze(A.Xpart(:,:,1));
y  = squeeze(A.Ypart(:,:,1));
% My bad and ignorance: I apply moving average along x, and y, with a
% number of nodes that grossly is equivalent to 200 km (two times more or
% less a typical slab). 

dx = diff(x(:,1));
dy = diff(y(1,:));
dx_MV = floor(200./mean(dx));
dy_MV = floor(200./mean(dy));
topo_M = movmean(topo,dx_MV,2);
topo_M = movmean(topo_M,dy_MV,1);
% Substract the median of the topography: 
% Here I would like to open a philosophical discussion on the significance
% of the base level and how to conceive a continental lithosphere. Atm, the
% average continental topography is above the sea level. One cool thing is
% adding the water to LaMEM solving the main issue with this matter. 
topo_M = topo_M-median(topo_M,"all");
Ph = squeeze(A.Phase(:,:,:));
T  = squeeze(A.Temp(:,:,:));
Z  = squeeze(A.Zpart(:,:,:));
%x  = squeeze(A.Xpart(:,1,:));
%y  = squeeze(A.Ypart(:,1,:));
ilx = length(x(:,1));
ily = length(y(1,:));
Ph2 = 0.*Ph;
T2  = 0.*Ph;
% loop over the 2D nodes, interpolate along z 
for i = 1:ilx
    for j =1:ily
        z = squeeze(Z(j,i,:))+topo_M(j,i);
        Ph2(j,i,:)=interp1(z,squeeze(Ph(j,i,:)),squeeze(Z(i,j,:)),'nearest');
        T2(j,i,:)=interp1(z,squeeze(T(j,i,:)),squeeze(Z(j,i,:)),'nearest');
        z = [];
    end
end
    Ph2(isnan(Ph2)& Z >0.0)=ph.Ph_Ar(1);
    Ph2(isnan(Ph2)& Z <0.0)=ph.Ph_UM(1);
    T2(isnan(T2) & Z >0.0)=TI.TS;
    T2(isnan(T2)& Z <0.0)=TI.TP;
    % = plot several section of the setup: temperature, phase and
    % topography and the new update topography 
    
    % for iy = 1:ily
    %     A.Phase(:,iy,:) = Ph2;
    %     A.Temp(:,iy,:) = T2;
    % end
    %

    [s]=save_topography(A,topo_M,Gr);
    surf(1,:)=x(:,1);
    surf(2,:)=topo_M;
    %disp(s);

end

    function [topo,rho] = compute_topography_(A,rho_ij)

    x = (squeeze(A.Xpart(:,1,1)));
    y = (squeeze(A.Ypart(1,:,1)));
    rho_x = 0.*x;
    topo  = 0.*x;

    for i = 1:length(x)
        for j=1:length(y)
            rho= nanmean(rho_ij(i,j,:));
            if i == 1 && j==1
                rho_0 = rho;
            end
            topo(i,j) = 1000.*(rho_0-rho)./rho;
            rho_x(i,j)  = (rho);
            Ph = [];
            l_ph = [];
            rho = [];
        end
    end
    end



    function [string] = save_topography(A,topo_M,Gr)


    dX = (max(Gr.x_g)-min(Gr.x_g))/1500;
    dY = (max(Gr.y_g)-min(Gr.y_g))/1500;
    x_t = min(Gr.x_g):dX:max(Gr.x_g);
    y_t = min(Gr.y_g):dX:max(Gr.y_g);
    topo_t=interp2(squeeze(A.Ypart(1,:,1)),squeeze(A.Xpart(:,1,1)),topo_M,y_t,x_t,'nearest','extrap');

    Topo = zeros(length(x_t),length(y_t));
    Topo(:,1)        = topo_t;
    Topo(:,2)        = topo_t;
    Easting     = x_t;
    Northing    = y_t;

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


    function [rho_ij] = Compute_density(A,ph)
        rho_ij = A.Xpart.*0.0; 
        
        % loop over the phases contained in ph
        field_names = fields(ph);
        ip = numel(field_names); 
        for i = 1:ip 
            CPU_A = cputime;
            % Select Phase
            current_phase = ph.(field_names{i});
            % Extract unrelevant information such as the reference density 
            Density_ref = current_phase(2); 
            % Compute the density as a function of the temperature 
            if current_phase(1) ~=0
                rho_ij(A.Phase==current_phase(1)) = density_node(Density_ref,A.Temp(A.Phase==current_phase(1)));
            else
                rho_ij(A.Phase==current_phase(1))=Density_ref;
                disp('Air does not deserve to be computed, but, as a quick reminder:'); disp('do not set 0 friction angle or cohesion to air, you will not notice my help.');disp(' But believe me KSP solver will be grateful and gives candy')
            end
            % Just random information for losing time on Monday, when I am
            % supposed to take care of myself. 
            CPU_B = cputime; 
            
            disp(['Phase nr. ', num2str(current_phase(1)),' has been processed in ',num2str(CPU_B-CPU_A,2), 'seconds']); 
        end
        disp('Additional note: density is computed assuming a thermal expansion of:');
        disp('3e-5 [1/K], if you want to change, just go to function density node, in displace_isostasy.m');
       
    end

    function [rho_n] = density_node(Density_ref,T)
        %=================================================================%
        % Hard coded stuff, but I doubt that anyone would really enjoy
        % himself on the daunting task to create ad hoc alpha for each
        % phase. There are worst crime in our community. For example melt
        % extraction or other funny approximation (like Venus with a
        % surface temperature of 0 degree celsius) 
        %=================================================================%
        % Prior: selected Phases -> select Markers -> send T
        % informationhere and phase information here: compute the density
        %=================================================================%
        rho_n = Density_ref.*(1-3e-5.*(T-25.0));
    end
