function [npt] = calculateWannierMoments(npt)


% calculate dipole matrix elements of wannier
X = npt.X_extended; % in lattice sites
X = X(:);

UX= zeros(size(npt.Wannier_X,2));

UX_HO = zeros(size(npt.Wannier_X,2));

UX2= zeros(size(npt.Wannier_X,2));
for mm=1:size(npt.Wannier_X,2)
    psi_mm = npt.Wannier_X(:,mm);
    for nn=1:size(npt.Wannier_X,2)
        psi_nn = npt.Wannier_X(:,nn);
        UX(mm,nn) = ctranspose(psi_mm)*(X.*psi_nn);
        % UX2(mm,nn) = ctranspose(psi_mm)*((X.^2).*psi_nn); % appears to diverage?

        % ladder operator gets an extra one beacuse of zero indexing
        UX_HO(mm,nn) = sqrt(nn+1-1)*(mm==(nn+1)) + sqrt(nn-1)*(mm==(nn-1));
    end    
end
UX=real(UX);
npt.WannierDipole = UX;
% npt.WannierDipole2 = UX2;
 npt.WannierDipole2 = UX*UX;

 npt.WannierDipole_Harmonic = UX_HO*npt.Harmonic_Length/sqrt(2);

end

