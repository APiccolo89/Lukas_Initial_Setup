function [Temp]= compute_temperature_profile(obj,A,ind,Temp)
% This is a public function that is accessed by several other classes. 
% obj: -> Terrane/Boundary that needs to be filled with temperature 
% D  : -> Depth {Z coordinates or distance from the top surface in case of subduction}
% ind : where to fill the temperature field. 

% T_tk = [0.0, obj.Stratigraphy.Tk(end)]; % Stratigraphy of the object
% kappa  = obj.kappa; %Thermal conductivity of the object 
% T_Age = obj.Age;    % Age of the object in [s]
% T_P   = obj.Thermal_information.TP; % Mantle potential temperature in C
% T_S   = obj.Thermal_information.TS; % Surface temperature in C
% Tsl   = obj.Thermal_type.AverageT; % Average Temperature in C 
% 
Type = obj.Thermal_Type.Type; 
if isa(class(obj),'Terrane')
    D = A.Zpart;  % For Layered structure Z coordinate make its job
elseif isa(class(obj),'trench')
    D = obj.Layout; % Layout is the member of the trench that collects all the distance from the top surface of the slab
end

if strcmp(Type,'Continental') % Is continental


elseif strcmp(Type,'HalfSpace')
    [Temp] = HalfSpaceCooling(obj,D,ind,Temp);
elseif strcmp(Type,'McKenzie')
    % McKenzie type require to run the temperature field of Halfspace
    % cooling model before. As the quick hack to have a smooth temperature
    % field is to introduce a weighted average as a function of the
    % decoupling depth. 
    [Temp] = HalfSpaceCooling(obj,D,ind,Temp);
    [Temp] = compute_temperature_profile_McKenzie(obj,D,ind,Temp);
elseif strcmp(Type,'Average_T')
    [Temp] = Average_T(obj,D,ind,Temp);
else
    error('You do not provide a temperature type:{Continental},{HalfSpace},{McKenzie},{Average_T}, please correct your ways and repent yourself.')
end

if Type == 1
    ind_z = find(D<T_tk(1) & D>=T_tk(2) & ind == 1);
    ind_o = A.Zpart< T_tk(2);

else
    ind_z = D < T_tk(1) & D>=T_tk(2);
    ind_o = D < T_tk(2);
    
end

if isnan(T_Age)
    Temp(ind_z)=Tsl;
else

    erf_function = erf(D(ind_z).*1000/2/(k*T_Age)^0.5);
    Temp(ind_z) = T_S - (T_P-T_S).*erf_function;
    Temp(ind_o) = T_P;

end

end
