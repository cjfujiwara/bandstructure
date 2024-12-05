function npt=calculateBandStructure(npt)

if nargin<1
    npt=constants;
    npt.depth=5;
    npt.numK=501;
    npt.numStates=25;
end

% Quasimomentum vector
K=npt.K;

% Basis and display size
numStates=npt.numStates;

% Data Vectors
bandsStatic0=zeros(numStates,length(K));       

vecStatic0=zeros(numStates,numStates,length(K));

nfo=struct;
nfo.depth=1;
nfo.numStates=numStates;
nfo.k = 0;

[~,pmat]=makeHmatrix(nfo); 

npt.PMatrix = pmat;


%% Calculate the band structure at each quasimomentum
clear myFigs
for nn = 1:length(npt.depth)
    % Lattice depth
    depth=npt.depth(nn);
    
    fprintf(['computing bands ' ...
        '(U=' num2str(depth) 'Er,' ...
        'Nk = ' num2str(npt.numK) ',' ...
        'Nstates = ' num2str(npt.numStates) ') ...']);
    t1=now;     
    
    for ii=1:length(K)    
        nfo.k=K(ii);
        nfo.depth=depth;
        [H0,~]=makeHmatrix(nfo);                % Hamiltonian
        [vS0,eng0]=eig(H0);                 % Solve
        bandsStatic0(:,ii)=diag(eng0);      % Assign energies           

      for cc=1:size(vS0,2)  
            if mod(cc,2) % even parity band
                vS0(:,cc)=vS0(:,cc)*exp(-1i * angle(vS0(1,cc)));
            else % odd parity band
                vS0(:,cc)=vS0(:,cc)*exp(-1i * angle(vS0(cc,cc)));

                vS0(:,cc)=1i*vS0(:,cc);
                vS0(:,cc)=(-1)^(cc/2+1)*vS0(:,cc);
            end
      end   
        vecStatic0(:,:,ii)=vS0;             % Assign eigenvectors     
    end
    t2=now;
    disp([' done (' num2str(round((t2-t1)*24*60*60,3)) ' s)']);

    % Add the computed band structure to the output
    npt.bandEigenValue(:,:,nn)=bandsStatic0;
    npt.bandEigenVectors(:,:,:,nn)=vecStatic0;  
end

end

