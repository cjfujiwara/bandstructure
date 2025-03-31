function [npt,hF]=calculateGapTunneling(npt)

%% Tunneling, Gap

U=linspace(.1,300,501);
U=U';

numStates=npt.numStates;
numStates=25;

BW=zeros(length(U),numStates);
Eds_0=zeros(length(U),1);
Eds_pi=zeros(length(U),1);
Eps_0=zeros(length(U),1);
Eps_pi=zeros(length(U),1);
%%%%%%%%%%%% calculate static lattice bands %%%%%%%%%%%%
fprintf('computing gaps ...');
t1=now;
for ii=1:length(U)     
    nfo=struct;
    nfo.depth=U(ii);
    nfo.numStates=numStates;
    
    k0=0;
    kpi=1;
    
    nfo.k=k0;    
    H0=makeHmatrix(nfo);              % get Hamiltonian
    [vS0,eng0]=eig(H0);               % calculate eigenvalues
    
    nfo.k=kpi;    
    Hpi=makeHmatrix(nfo);             % get Hamiltonian
    [vSpi,engpi]=eig(Hpi);            % calculate eigenvalues     
    
    % Bandwidth
    BW(ii,:)=diag(engpi-eng0);
    
    Eds_0(ii) = eng0(3,3)-eng0(1,1);
    Eds_pi(ii) = engpi(3,3)-engpi(1,1);
    
    Eps_0(ii) = eng0(2,2)-eng0(1,1);
    Eps_pi(ii) = engpi(2,2)-engpi(1,1);
end
t2=now;
fprintf([' done (' num2str(round((t2-t1)*24*60*60,3)) 's)']);
 
%% Apend to output
bandwidth=struct;
bandwidth.depth=U;
bandwidth.bandwidth=BW;
bandwidth.tunneling=BW;npt.fr/4;

bandgap=struct;
bandgap.depth=U;
bandgap.Eps=[Eps_0 Eps_pi];
bandgap.Eds=[Eds_0 Eds_pi];

npt.bandwidth=bandwidth;
npt.bandgap=bandgap;

%% Plot Results
t1=now;
fprintf(' plotting ...');

hF=figure(224);
set(hF,'color','w','name','gap_tunneling');
hF.Position=[250 50 300 200];
clf

pS=plot(U,abs(BW(:,1))*npt.fr/4,'LineWidth',2); 
hold on
pP=plot(U,abs(BW(:,2))*npt.fr/4,'LineWidth',2); 

ylim([0 1000]);


xlabel('lattice depth (E_R)');
ylabel('tunneling rate (Hz)');
legend([pS,pP],{'$+\Delta_s/4$','$-\Delta_p/4$'},'interpreter','latex',...
    'fontsize',10);
set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on');

xlim([0 20]);
% Plot the probability density
hF4=figure(225);
set(hF4,'color','w','name','gap_ds');
hF4.Position=[550 50 300 200];
clf

pS=plot(U,Eds_0*npt.fr*1e-3,'k-','LineWidth',2); 
hold on
pP=plot(U,Eds_pi*npt.fr*1e-3,'k--','LineWidth',2); 

pHO=plot(U,2*sqrt(4*U),'k--','LineWidth',2); 

xlabel('lattice depth (E_R)');
ylabel('energy (kHz)');
legend([pS,pP],{'$E_{ds}(0)$','$E_{ds}(\pi)$'},'interpreter','latex',...
    'fontsize',10,'location','southeast');
set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on');
xlim([0 20]);

% Plot the probability density
hF4b=figure(230);
set(hF4b,'color','w','name','gap_ps');
hF4b.Position=[850 50 300 200];
clf

pS=plot(U,Eps_0*npt.fr*1e-3,'k-','LineWidth',2); 
hold on
pP=plot(U,Eps_pi*npt.fr*1e-3,'k--','LineWidth',2); 

pHO=plot(U,sqrt(4*U)*npt.fr*1e-3,'k:','LineWidth',2); 

% pHO=plot(U,sqrt(4*U)*npt.fr*1e-3,'r:','LineWidth',2); 

xlabel('lattice depth (E_R)');
ylabel('energy (kHz)');
legend([pS,pP,pHO],{'$E_{ps}(0)$','$E_{ps}(\pi)$','harmonic'},'interpreter','latex',...
    'fontsize',10,'location','southeast');
set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on');
t2=now;
fprintf([' done (' num2str(round((t2-t1)*24*60*60,3)) 's)']);
disp(' ');


% Plot the probability density
hF5b=figure(2311);
set(hF5b,'color','w','name','gap_ps');
hF5b.Position=[1150 50 300 200];
clf

yBand = 0.5*(Eps_pi+Eps_pi);
yHO   = sqrt(4*U);

pP=plot(U,yBand./yHO,'k-','LineWidth',1); 


xlabel('lattice depth (E_R)');
ylabel('energy (kHz)');
legend([pP],{'$\Delta_{ps}/\Delta_\mathrm{HO}$'},'interpreter','latex',...
    'fontsize',10,'location','southeast');
set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on');
t2=now;
fprintf([' done (' num2str(round((t2-t1)*24*60*60,3)) 's)']);
disp(' ');
end

