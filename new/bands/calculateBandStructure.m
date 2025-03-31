function npt=calculateBandStructure(npt)

if nargin<1
    npt=constants;
    npt.depth=5;
    npt.numK=501;
    npt.numStates=25;
end

K           = npt.K;            % quasimomentum vector
numStates   = npt.numStates;    % plane wave basis size
eng         = zeros(numStates,length(K));           % Eigenvalues
vec         = zeros(numStates,numStates,length(K)); % Eigenvectors init

% Initialize Other Stuff
nfo             = struct;
nfo.depth       = 1;
nfo.numStates   = numStates;
nfo.k           = 0;
[~,pmat]=makeHmatrix(nfo); 
npt.PMatrix = pmat;


%% Calculate the band structure at each quasimomentum
for nn = 1:length(npt.depth)
    % Lattice depth
    depth=npt.depth(nn);
    
    fprintf(['computing bands ' ...
        '(U=' num2str(depth) 'Er,' ...
        'Nk = ' num2str(npt.numK) ',' ...
        'Nstates = ' num2str(npt.numStates) ') ...']);
    tic    
    for ii=1:length(K)    
        nfo.k       = K(ii);                % quasimomentum 
        nfo.depth   = depth;                % depth
        [H0,~]      = makeHmatrix(nfo);     % Hamiltonian
        [vS0,eng0]  = eig(H0);              % Solve
        eng(:,ii)   = diag(eng0);           % Assign energies           

        % Modify Eigenvectors for Sign
        for cc=1:size(vS0,2)  
            if mod(cc,2) % even parity band
                vS0(:,cc)=vS0(:,cc)*exp(-1i * angle(vS0(1,cc)));                
                % Forcing to be real, attempting
                vS0(:,cc) = abs(vS0(:,cc)).*sign(real(vS0(:,cc)));
            else % odd parity band
                vS0(:,cc)=vS0(:,cc)*exp(-1i * angle(vS0(cc,cc)));
                vS0(:,cc)=1i*vS0(:,cc);
                vS0(:,cc)=(-1)^(cc/2+1)*vS0(:,cc);
            end
        end   
        vec(:,:,ii)=vS0;             % Assign eigenvectors     
    end
    t2=toc;
    disp([' done (' num2str(round(t2,3)) ' s)']);
    % Add the computed band structure to the output
    npt.bandEigenValue(:,:,nn)=eng;
    npt.bandEigenVectors(:,:,:,nn)=vec;  
end

end

