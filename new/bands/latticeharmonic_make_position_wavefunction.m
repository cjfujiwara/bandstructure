function [psi] = latticeharmonic_make_position_wavefunction(eigen_vector,eigen_position,eigen_band,wannier_vector,wannier_position)

if max(eigen_position)>max(wannier_position) 
    error('need larger wannier mesh')
end
if min(eigen_position)<min(wannier_position) 
    error('need larger wannier mesh')
end

% Number of elements to shift perlattice site
N = round(1/mode(diff(wannier_position)));

psi = zeros(length(wannier_position),1);


for jj=1:length(eigen_vector)
    J = eigen_position(jj);
    n = eigen_band;
    psi = psi + eigen_vector(jj)*circshift(wannier_vector(:,n),N*J);
end

end

