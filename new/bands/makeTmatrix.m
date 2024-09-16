function Tmat = makeTmatrix(n,jmax)
% Kinetic energy tunneling operator with periodic boundary conditions;
%
%
    Tmat = zeros(n,n,jmax);
    for kk=1:jmax
        Tmat(:,:,kk) = diag(ones(n-kk,1),kk)+diag(ones(n-kk,1),-kk);
        for jj  = 1:kk
            Tmat(n-(jj-1),kk-jj+1,kk) = 1;
            Tmat(kk-jj+1,n-(jj-1),kk) = 1;
        end
    end

    Tmat = -Tmat;
end

