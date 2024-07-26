function [npt,hF]=getBandStructure(npt,opts)

if nargin<1
    npt=constants;
    npt.depth=5;
    npt.numK=501;
    npt.numStates=25;
end

if nargin<2
   opts=struct;
   opts.doPlot = 1;
end



% Quasimomentum vector
K=npt.K;
dK=K(2)-K(1);

% Basis and display size
numStates=npt.numStates;
numBands=4;

% Data Vectors
% bandsStatic0=zeros(length(K),numStates);   
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
    fprintf([' done (' num2str(round((t2-t1)*24*60*60,3)) ' s)']);

    % Add the computed band structure to the output
    npt.bandEigenValue(:,:,nn)=bandsStatic0;
    npt.bandEigenVectors(:,:,:,nn)=vecStatic0;

    % keyboard

    %% Plot Band structure
    
    if opts.doPlot
        t1=now;
        % Plot the static band structure, harmonic oscillator energies, and the
        % transitions
        fprintf(' plotting ...');

        % Make the figure
        hF=figure(222);
        set(hF,'Name','band_structure','color','w');
        clf;
        hF.Position=[0 50 250 400];

        % Initialize the axis
        ax1=axes;
        cla
        co=get(gca,'colororder');
        set(ax1,'fontsize',10,'box','on','linewidth',1,'fontname','times');
        xlabel('quasimomentum ($\hbar k_L$)','interpreter','latex');
        ylabel('energy ($E_R$)','interpreter','latex');
        xlim([min(K) max(K)]);
        ylim([-nfo.depth-1 ceil(max(max(bandsStatic0(1:numBands,:))))]);
        hold on

        % Plot harmonic oscilator energies
        engHO=-nfo.depth+2*sqrt(nfo.depth)*(0.5+(0:10));
        for ii=1:length(engHO)
            plot([-1 1],[1 1]*engHO(ii),'k:','linewidth',1); 
        end

        % Plot the bands
        for kk=1:numBands
           plot(K,bandsStatic0(kk,:),'-','linewidth',3,...
               'color',co(mod(kk-1,7)+1,:)); 
        end

        % Add depth label
        str=['$' num2str(depth) '~E_R$'];
        text(5,5,str,'interpreter','latex','units','pixels',...
            'verticalalignment','bottom','fontsize',16);
        t2=now;
        fprintf([' done (' num2str(round((t2-t1)*24*60*60,3)) ' s)']);

    else
        hF=[];
    end
    disp(' ');
end

end

