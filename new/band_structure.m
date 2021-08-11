%% Static Lattice Calculations
K=linspace(-1,1,501);
numStates=21;
numBands=5;
depth=50;

bandsStatic0=zeros(length(K),numStates);       % static bands
bandsFold=zeros(length(K),numStates);

vecStatic0=zeros(numStates,numStates,length(K));


%%%%%%%%%%%% calculate static lattice bands %%%%%%%%%%%%
fprintf('computing static band structure...');
for ii=1:length(K)     
    nfo=struct;
    nfo.depth=depth;
    nfo.numStates=numStates;
    nfo.k=K(ii);
    
    H0=makeHmatrix(nfo);              % get Hamiltonian
    [vS0,eng0]=eig(H0);               % calculate eigenvalues
    bandsStatic0(ii,:)=diag(eng0);
    vecStatic0(:,:,ii)=vS0;
end
disp('done');

% Add the computed band structure to the output
output.bandsStatic=bandsStatic0;
output.vecStatic=vecStatic0;

%% Figure 1: Static Bands
% Plot the static band structure, harmonic oscillator energies, and the
% transitions
fprintf('Plotting static bands with transitions...');

hF1=figure('Name','static_bands','color','w',...
    'units','pixels');
clf;
hF1.Position=[10 50 350 600];

ax1=axes;
cla
co=get(gca,'colororder');
set(ax1,'fontsize',14,'box','on','linewidth',1,'fontname','times');
xlabel('quasimomentum ($\hbar k_L$)','interpreter','latex');
ylabel('energy ($E_R$)','interpreter','latex');
xlim([min(K) max(K)]);
ylim([-nfo.depth ceil(max(max(bandsStatic0(:,1:numBands))))]);
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


disp('done');

%% Wannier Function
n=1;
nSite=0;

figure(10)
clf

% position vector
xVec=linspace(-pi,pi,1000);

% wanier function
wData=zeros(1,length(xVec));
for kk=1:length(K)
    % Get basis vectors
    v=vecStatic0(:,:,kk);
    
    % Get weights for this band index
    vn=v(:,n);
        
    % Rotate vectors to have real positive coefficient of exp(2 i phi)
    vn = vn * exp(-1i * angle(vn(2)));
        
    u_nk   = makeU_nk(vn);
    psi_nk = u_nk(xVec).*exp(1i*K(kk)*xVec).*exp(-1i * K(kk) * pi*nSite); 
        
    % Normalize the wavefunction
    N=sum(psi_nk.*conj(psi_nk));    
    psi_nk=psi_nk/sqrt(N);
        
    wData=wData+psi_nk;    
end

wData=wData/sqrt(length(K));

N=sum(wData.*conj(wData));    
wData=wData/sqrt(N);

p=plot(xVec/pi,real(wData.*conj(wData)),'LineWidth',2); 

% 
HO=exp(-xVec.^2*sqrt(depth)/2);
N=sum(HO.*conj(HO));  
HO=HO/sqrt(N);

hold on
p=plot(xVec/pi,real(HO.*conj(HO)),'LineWidth',2); 

