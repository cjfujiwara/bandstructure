function [npt,hF]=getBandStructure(npt)

if nargin==0
    npt=constants;
    npt.depth=5;
    npt.numK=501;
    npt.numStates=25;
end

% Lattice depth
depth=npt.depth;

% Quasimomentum vector
K=npt.K;
dK=K(2)-K(1);

% Basis and display size
numStates=npt.numStates;
numBands=3;

% Data Vectors
bandsStatic0=zeros(length(K),numStates);       
vecStatic0=zeros(numStates,numStates,length(K));

%% Calculate the band structure at each quasimomentum
fprintf('Computing static band structure...');
t1=now;
for ii=1:length(K)     
    nfo=struct;
    nfo.depth=depth;
    nfo.numStates=numStates;
    nfo.k=K(ii);
    
    H0=makeHmatrix(nfo);                % Hamiltonian
    [vS0,eng0]=eig(H0);                 % Solve
    bandsStatic0(ii,:)=diag(eng0);      % Assign energies   
    
    % Rotate each vectors to have real positive coefficient of exp(2 i phi)
    % or 1, (parity of band)
    for cc=1:size(vS0,2)   
        ind = 2-mod(cc,2);
        vS0(:,cc)=vS0(:,cc)*exp(-1i * angle(vS0(ind,cc)));
    end          
    vecStatic0(:,:,ii)=vS0;             % Assign eigenvectors        
end
t2=now;
disp(['done (' num2str(round((t2-t1)*24*60*60,3)) ' s)']);

% Add the computed band structure to the output
npt.bandEigenValue=bandsStatic0;
npt.bandEigenVectors=vecStatic0;


%% Plot Band structure

% Plot the static band structure, harmonic oscillator energies, and the
% transitions
fprintf('Plotting static bands with transitions...');

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
ylim([-nfo.depth-1 ceil(max(max(bandsStatic0(:,1:numBands))))]);
hold on

% Plot harmonic oscilator energies
engHO=-nfo.depth+2*sqrt(nfo.depth)*(0.5+(0:10));
for ii=1:length(engHO)
    plot([-1 1],[1 1]*engHO(ii),'k:','linewidth',1); 
end

% Plot the bands
for kk=1:numBands
   plot(K,bandsStatic0(:,kk),'-','linewidth',3,...
       'color',co(mod(kk-1,7)+1,:)); 
end

% Add depth label
str=['$' num2str(depth) '~E_R$'];
text(5,5,str,'interpreter','latex','units','pixels',...
    'verticalalignment','bottom','fontsize',16);

disp('done');
end

