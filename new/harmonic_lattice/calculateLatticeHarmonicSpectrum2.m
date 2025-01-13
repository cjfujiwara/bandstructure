function [npt,output] = calculateLatticeHarmonicSpectrum2(npt,opts)
n = opts.NumSites;
jjMax = opts.MaxTunnelingOrder;
Tmat = makeTmatrix(n,jjMax);
Omega = opts.Omega;
nBands = opts.NumBands;
fr = npt.fr;

%% Construct Harmonic Position Operator
r = (n-1)/2;        % [-r,r] is the position in size
x = -r:1:r;         % position operator is symmetric
x2 = x.^2;          % position squared

% Change of basis operator
U = make_U_siteToOddEven(opts.NumSites);

%% Initialize ouput
output = struct;
output.Depth = npt.depth;
output.Tunneling = npt.Tunneling;
output.NumSites = n;
output.NumBands = nBands;
output.MaxTunnelingOrder = jjMax;
output.Omega = Omega;   
output.omega = opts.omega;
output.fr = fr;
output.PositionVector = repmat(x,[1 nBands]);

%% Interate over Lattice depths
output.EigenValues = zeros(output.NumSites*output.NumBands,length(npt.depth));
output.EigenVectors = zeros(output.NumSites*output.NumBands,output.NumSites*output.NumBands,length(npt.depth));

output.BandRanges = zeros(output.NumBands,2,length(npt.depth));
output.BandProjection = zeros(output.NumSites*output.NumBands,output.NumBands,length(npt.depth));

for uu = 1:length(npt.depth)
    tic
    fprintf(['computing lattice+harmonic (U=' num2str(npt.depth(uu)) 'Er)']);
    t   = npt.Tunneling(:,:,uu)*npt.fr; % Tunneling matrix elements   

    %% Kinetic Energy operator for each band
    T = zeros(n,n,nBands);
    for nn = 1:nBands % Iterate over all bands
        for jj = 1:jjMax % Iterate over all tunneling order
            T(:,:,nn) = T(:,:,nn) + t(nn,jj)*Tmat(:,:,jj);           
        end
    end

    % Change Basis fo Kinetic Energy
    for nn=1:nBands
        T(:,:,nn)=U*T(:,:,nn)*ctranspose(U);
    end
    T(abs(T)<1e-6)=0; % Reduce numerical noise
       
    %%
    % Potential Energy Operator for each band
    X1 = diag(x);
    X2 = diag(x2);

    % Dipole matrix operator band x band x delta sites
    mdip = ((npt.m1+npt.m2)/2)/pi; % normalized to lattice spacing

    % Construct potential energy operator w band couplings
    alpha = npt.WannierDipole; % couples different bands same site
    Vdiag = zeros(n*nBands);
    Vdip = zeros(n*nBands);
    Vdip2 = zeros(n*nBands);


    n_dipole_hopping = size(mdip,3);
    pl = -(n_dipole_hopping-1)/2:(n_dipole_hopping-1)/2;

    for rr = 1:nBands
        for cc = 1:nBands
            Vsub = zeros(length(x),length(x));
            for nI=1:n_dipole_hopping
                dI = pl(nI);
                dI = abs(dI);
                coupling_rc_dI = mdip(rr,cc,nI);                
                n1 = diag(x((1+dI):end),dI);
                n2 = diag(x(1:(end-dI)),-dI);
                N_hop = n1 + n2;
                mat1 = coupling_rc_dI*N_hop;
                Vsub = Vsub + mat1;
            end

            % Indeces in super matrix
            r1 = n*(rr-1)+1;r2 = n*rr;
            c1 = n*(cc-1)+1;c2 = n*cc;  
            Vdip(r1:r2,c1:c2) = Vsub;
            if rr==cc
                Vdiag(r1:r2,c1:c2) = X2;
            end
        end
    end
    V = (Vdiag + Vdip + Vdip2)*Omega;
%% asdf 
    % Energy offset Opeator for each band
    E0 = zeros(n,n,nBands);
    for nn = 1:nBands
        % Ebar = mean(npt.bandEigenValue(nn,:,uu))*fr;
        Ebar = npt.bandEigenValueAverage(nn,uu)*fr;
        output.BandRanges(nn,1,uu) = min(npt.bandEigenValue(nn,:,uu))*fr; 
        output.BandRanges(nn,2,uu) = max(npt.bandEigenValue(nn,:,uu))*fr; 
        E0(:,:,nn) = eye(n)*Ebar;
    end

    % Lowest Band
    H = T(:,:,1) + V(:,:,1) + E0(:,:,1);

    % Add on additional bands
    for nn = 2:nBands
        H = blkdiag(H,T(:,:,nn)+V(:,:,nn)+E0(:,:,nn));       
    end
%% asf

    % Lowest Band
    H = T(:,:,1)  + E0(:,:,1);
    % Add on additional bands
    for nn = 2:nBands
        H = blkdiag(H,T(:,:,nn)+E0(:,:,nn));       
    end
    H = H + V;



    % Solve Eigenvalue problem
    [C,D]=eig(H);    
    



    [Ed,inds] = sort(real(diag(D)));
    C = C(:,inds);

    output.EigenValues(:,uu) = Ed;
    output.EigenVectors(:,:,uu) = C;

    for nn = 1:nBands
        i1 = n*(nn-1)+1;
        i2 = n*nn;
        output.BandProjection(:,nn,uu) = sum(abs(C(i1:i2,:)).^2,1);
    end
    disp(['done( ' num2str(toc,3) 's)'])
end



end

% Create a change of basis matrix that changes from basis of odd+even combination of sites 
% to the basis of sites
function U = make_U_siteToOddEven(N)
% N : size of basis (assumed to be odd to include n=0);
    nc = (N+1)/2;               % Index of n=0;    
    % Change of basis matrix
    u1 = eye(N);
    u2 = flip(u1,1);
    u1(nc:end,nc:end)=-u1(nc:end,nc:end);
    U=u1+u2;
    U=U/sqrt(2);
    U(nc,nc)=1;
end

