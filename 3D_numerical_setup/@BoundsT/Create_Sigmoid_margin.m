function [obj] = Create_Sigmoid_margin(obj,Boundary,x0,B,v,K1,K2)
% Create all the other data.
obj=obj.compute_coordinate_boundary_composite;
% This part of the code should be cleansed whenever someone has time. The
% goals justify the means, and I would like to finish the 3D project by the
% end of January with already part of the manuscript written
if strcmp(Boundary,'A')
elseif strcmp(Boundary,'B')
    obj.B{2} = 'Sigmoid'; %function handle boundary
    obj.B{3} = [x0,K1,K2,B,v]; % save the important information.  
elseif strcmp(Boundary,'C')
elseif strcmp(Boundary,'D')
    obj.D{2} = 'Sigmoid'; %function handle boundary
    obj.D{3} = [x0,K1,K2,B,v]  ; % save the important information.  
end
end
