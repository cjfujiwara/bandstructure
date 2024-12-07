function [npt,output] = calculateLatticeHarmonicSpectrumOdd(npt,opts)
n = opts.NumSites;
jjMax = opts.MaxTunnelingOrder;
Tmat = makeTmatrixNOPBC(n,jjMax);
Omega = opts.Omega;
nBands = opts.NumBands;
fr = npt.fr;


%% Construct Harmonic Position Operator
x = 0:1:(n-1);         % position operator is symmetric
x2 = x.^2;          % position squared

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

    % Kinetic Energy operator for each band
    T = zeros(n,n,nBands);
    for nn = 1:nBands % Iterate over all bands
        for jj = 1:jjMax % Iterate over all tunneling order
            T(:,:,nn) = T(:,:,nn) + t(nn,jj)*Tmat(:,:,jj);           
        end
    end
        
    %% Diagonal Potential Energy Operator       
    V = diag(x2)*Omega;
    
    % Add on infinite energy at origin (forces a node)
    V(1,1)=1e10;
%       V(1,1)=-1e9;

    %% Construct Hamiltonian
    H = T(:,:,1)+V;
    
    %% Solve Eigenvalue problem  
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

