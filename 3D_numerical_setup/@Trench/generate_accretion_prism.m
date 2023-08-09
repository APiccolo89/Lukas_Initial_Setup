function  [Phase] = generate_accretion_prism(obj,A,Phase,Boundary)
x = A.Xpart(:);
z = A.Zpart(:);
data_boundary = obj.Boundary.(Boundary);
%Select the boundary = limits of the boundary
[loc1, ~] = find(tril(data_boundary{1} == data_boundary{1}', -1)); % find the main boundary.
if rem(loc1,2)==0
    xc = data_boundary{1}(2);
else
    xc = data_boundary{1}(3);
end

C = [obj.Boundary(1),-obj.R];
d_p = [C(1)+obj.position_end_prism 0.0];
s   = (d_p(2)-C(2))./(d_p(1)-C(1));
ind2 = z(:)>s.*(x(:)-C(1))+C(2)  & (Phase(:) == ~isnan(Phase(:)) | Phase(:) ~= obj.Thermal_information.Ph_Ast) & isnan(obj.d_slab(:)) &z <0.0 & x>C(1);
Phase(ind2) = obj.phase_prism;

end
