function [N,etot]=calculateThermodynamics(npt,x,y,z)
out = struct;

Ex = [x.EigenValues];
Ey = [y.EigenValues];
Ez = [z.EigenValues];

% Trap Frequency and effective band mass for effective frequency
meff = npt.BandCurvatureG(1);
fx = x.omega/(2*pi);
fy = y.omega/(2*pi);
fz = z.omega/(2*pi);
fbar = (fx*fy*fz)^(1/3)*sqrt(meff);

% Use fitted slope for effective trap frequency
fx = x.FirstBandLinearSlope;
fy = y.FirstBandLinearSlope;
fz = z.FirstBandLinearSlope;
fbar= (fx*fy*fz)^(1/3);

% energy of next band relative to bottom of first band
Ep = npt.fr*(min(npt.bandEigenValue(2,:))-min(npt.bandEigenValue(1,:)));

%Band width
BW = npt.fr*(max(npt.bandEigenValue(1,:))-min(npt.bandEigenValue(1,:)));

% Which states to look at
nmax = [300 300 150];

% Consturcust energy
[EX,EY,EZ]=meshgrid(Ex(1:nmax(1)),Ey(1:nmax(2)),Ez(1:nmax(3)));
Etot = EX+EY+EZ;
Etot=Etot(:);
[Etot,inds] = sort(Etot);

% State indeces
[sx,sy,sz]=meshgrid(1:nmax(1),1:nmax(2),1:nmax(3));
sx=sx(inds);
sy=sy(inds);
sz=sz(inds);


% Band indeces
bx = repmat(1:x.NumBands,[nmax(1) 1]);
nx = sum(x.BandProjection(1:nmax(1),:).*bx,2);
by = repmat(1:y.NumBands,[nmax(2) 1]);
ny = sum(x.BandProjection(1:nmax(2),:).*by,2);
bz = repmat(1:z.NumBands,[nmax(3) 1]);
nz = sum(x.BandProjection(1:nmax(3),:).*bz,2);

[NX,NY,NZ]=meshgrid(nx,ny,nz);
NX = NX(:);
NX = NX(inds);
NY = NY(:);
NY = NY(inds);
NZ = NZ(:);
NZ = NZ(inds);

Etot = Etot-Etot(1);
N=1:length(Etot);

out.Index = N;
out.Energy = Etot;


Nsum = NX+NY+NZ;
Nc = round(Nsum);
Nc = Nc(:);

out.BandSum = Nc;


% keyboard
t1= npt.Tunneling(1,1)*npt.fr;

%% Plotting
f=figure(20);
f.Color='w';
f.Position = [50 50 600 300];
co=get(gca,'colororder');
clf
isub = 1:2e5;


% Harmonic approximation
N2e = @(N) (6*N).^(1/3)*fbar;

% Using actual combinatorics
e2N = @(n) (n+1).*(n+2).*(n+3)/6;

Evec = linspace(0,max(Etot),1e4);
neff=Evec/fbar;
Neff = e2N(neff);

titstr = ['$' num2str(npt.depth) 'E_R, ' ...
    '(' num2str(round(x.omega/(2*pi),1)) ',' ...
    num2str(round(y.omega/(2*pi),1)) ',' ...
    num2str(round(z.omega/(2*pi),1)) ')\mathrm{Hz},' ...
    num2str(round(meff,2)) 'm_0' '$'];

pho = plot(N(isub),N2e(N(isub)-1),'k-','linewidth',1);
hold on

% pdata = plot(N(isub),Etot(isub),'.','linewidth',2,'color',co(1,:));
pdata = scatter(N(isub),Etot(isub),5,co(Nc(isub)-2,:),'linewidth',2);

xlabel('eigen index');
ylabel('energy (Hz)');
str1=['harmonic w/ $m_\mathrm{eff}$ $\epsilon = \bar{f}(6N)^{1/3}$' newline '$f=' num2str(round(fbar,1)) '\mathrm{Hz}$'];
str2=['numerical spectrum'];
set(gca,'box','on','fontsize',12,'fontname','times','linewidth',1);
ylim([0 3.5*BW]);
drawnow;
pP=plot(get(gca,'XLim'),[1 1]*Ep,'-','color',co(2,:));
p1=plot(get(gca,'XLim'),[1 1]*BW,'k-.');
p2=plot(get(gca,'XLim'),[1 1]*2*BW,'k--');
p3=plot(get(gca,'XLim'),[1 1]*3*BW,'k-');
legStr={str1, str2,'min$(E_p)$' , '1 bandwidth','2 bandwidth','3 bandwidth'};
legend([pho pdata pP p1 p2 p3],legStr,'interpreter','latex','location','southeast');
title(titstr,'interpreter','latex');
drawnow;

%% Chemical Potential, Total Energy Entropy
% Etot = Etot(1:1e6);

disp('Calcuating chemical potential');
N0 = 50e3;          % Particle number
E_fermi = Etot(N0); % Fermi Energy
Enorm = Etot/E_fermi;


% T/Tf vector
Tvec = linspace(0,.5,100);
Tvec = Tvec(:);


% Calculate the chemical potential
muGuess = 1;
muVec = zeros(length(Tvec),1);
muVec(1) = 1;

engTot = zeros(length(Tvec),1);
for kk = 2:length(Tvec)
    fprintf([num2str(kk) '/' num2str(length(Tvec)) ' ' num2str(Tvec(kk)) ' ...']);
    myfunc = @(mu) N0-mu2N(Enorm,Tvec(kk),mu);
    muThis = fzero(myfunc,muGuess);
    muVec(kk) = muThis;
    muGuess= muThis;    
    engTot(kk) = mu2E(Enorm,Tvec(kk),muThis);
    disp(['mu = ' num2str(round(muThis,3)) ' done']);
end
engTot(1) = mu2E(Enorm,1e-9,1);
out.T2mu = @(x) interp1(Tvec,muVec,x);


% Fugacity
z = exp(muVec./Tvec);
% Entropy
fprintf('calculating entropy ... ')
entropy = zeros(length(Tvec),1);
for kk=1:length(Tvec)
    mu = muVec(kk);
    s0 = (engTot(kk)-mu*N0)/Tvec(kk);
    s1 = sum(log(1+exp(-(Enorm-mu)/Tvec(kk))));
    entropy(kk) = (s0+s1)/N0;
end
disp('done')

%% Harmonic Oscillator Chemical Potential

fprintf('calculating fugacity harmonic ...')
logZ = linspace(0,1000,1000);
logZ=logZ(:);
zho = exp(logZ);
Tho = (-factorial(3)*polylog(3,-zho)).^(-1/3);
Tho=Tho(:);
disp('done');

%% Plot Chemical Potential

kBh = 20.84; % Hz/nK

f2 = figure(21);
f2.Color='w';
f2.Position=[5 100 600 300];
clf
pHO=plot(Tho,logZ.*Tho,'k-','linewidth',1);
hold on
pData=plot(Tvec,muVec,'.','color',co(1,:));
hold on

xlabel('temperature $T/T_\mathrm{F}$','interpreter','latex')
ylabel('chemical potential $\mu/T_\mathrm{F}$','interpreter','latex')
set(gca,'fontsize',12,'box','on','linewidth',1,'fontname','times');
title(titstr,'interpreter','latex');
str1=['$N=' num2str(N0) ', E_F= ' num2str(round(E_fermi/1e3,2)) ' \mathrm{kHz} (' num2str(round(E_fermi/kBh)) ' \mathrm{nK})$'];
text(.01,.01,str1,'fontsize',14,'units','normalized','interpreter','latex',...
    'horizontalalignment','left','verticalalignment','bottom');
legend({'harmonic','actual'})

%% Fugacity
f3 = figure(22);
f3.Color='w';
f3.Position=[610 100 600 300];
clf
plot(Tvec,z,'.','color',co(1,:));

xlabel('temperature $T/T_\mathrm{F}$','interpreter','latex')
ylabel('fugacity  $z$','interpreter','latex')
set(gca,'fontsize',12,'box','on','linewidth',1,'fontname','times',...
    'yscale','log');
title(titstr,'interpreter','latex');
str1=['$N=' num2str(N0) ', E_F= ' num2str(round(E_fermi/1e3,2)) ' \mathrm{kHz} (' num2str(round(E_fermi/kBh)) ' \mathrm{nK})$'];
text(.01,.01,str1,'fontsize',14,'units','normalized','interpreter','latex',...
    'horizontalalignment','left','verticalalignment','bottom');
ylim([0 1e4]);

%% Energy
f4 = figure(23);
f4.Color='w';
f4.Position=[1215 100 600 300];
clf


plot(Tvec,engTot/N0,'.','color',co(1,:));

xlabel('temperature $T/T_\mathrm{F}$','interpreter','latex')
ylabel('$\mathcal{E}/N~(\epsilon_\mathrm{F})$','interpreter','latex')
set(gca,'fontsize',12,'box','on','linewidth',1,'fontname','times',...
    'yscale','linear');
title(titstr,'interpreter','latex');
str1=['$N=' num2str(N0) ', E_F= ' num2str(round(E_fermi/1e3,2)) ' \mathrm{kHz} (' num2str(round(E_fermi/kBh)) ' \mathrm{nK})$'];
text(.01,.01,str1,'fontsize',14,'units','normalized','interpreter','latex',...
    'horizontalalignment','left','verticalalignment','bottom');
% ylim([0 1e4]);
%% Entropy
f5 = figure(24);
f5.Color='w';
f5.Position=[15 410 600 300];
clf

SHO = 4*polylog(4,-zho)./polylog(3,-zho)-log(zho);

pHO=plot(Tho,SHO,'k-','linewidth',1);
hold on

pHOapprox=plot(Tho,pi^2*Tho,'k--','linewidth',1);


pNumeric = plot(Tvec,entropy,'.','color',co(1,:));

xlabel('temperature $T/T_\mathrm{F}$','interpreter','latex')
ylabel('$S/N~(\mathrm{k_B})$','interpreter','latex');
set(gca,'fontsize',12,'box','on','linewidth',1,'fontname','times',...
    'yscale','linear');
title(titstr,'interpreter','latex');
str1=['$N=' num2str(N0) ...
    ', E_F= ' num2str(round(E_fermi/1e3,2)) ' \mathrm{kHz} (' ...
    num2str(round(E_fermi/kBh)) ' \mathrm{nK})$' newline  ...
    '$t=' num2str(round(t1,1)) '~\mathrm{Hz}, E_F/t=' num2str(round(E_fermi/t1,2)) '$'];
text(.99,.01,str1,'fontsize',14,'units','normalized','interpreter','latex',...
    'horizontalalignment','right','verticalalignment','bottom');

legend([pHO pHOapprox pNumeric],{'harmonic', '$\pi^2 T/T_\mathrm{F}$','numeric'},...
    'interpreter','latex','fontsize',12,'location','northwest');
xlim([0 .4]);

%% Example Occupation
myT = 0.15;
myMu = interp1(Tvec,muVec,myT);
myP = mu2P(Enorm,myT,myMu);
myN=sum(myP);

NNX = round(NX);
NNY = round(NY);
NNZ = round(NZ);

Nxp = sum(myP.*[NNX==2]);
Nyp = sum(myP.*[NNY==2]);
Nzp = sum(myP.*[NNZ==2]);
Np = round(Nxp + Nyp + Nzp);

N3dlocal = sum(myP.*[sy>=45].*[sx>=45].*[sz>=15])

Nylocal = sum(myP.*[sy>=45]);
Nxlocal = sum(myP.*[sx>=45]);
Nzlocal = sum(myP.*[sz>=15]);

Nactive = sum(myP.*[sy<45].*[sx<45].*[sz<15]);

iXY_plane = [sz == 1];
iXZ_plane = [sy == 1];

iX = logical(iXY_plane.*iXZ_plane);
X = 1:nmax(1);


%%
strS = ['$N= ' num2str(round(myN))  ', T = ' num2str(myT) 'T_\mathrm{F}$'];
strEF = ['$\epsilon_\mathrm{F} = ' num2str(round(E_fermi)) '~\mathrm{Hz}$'];
strMu = ['$\mu = ' num2str(round(myMu*E_fermi)) '~\mathrm{Hz}$'];
strBW = ['$3\mathrm{BW} = ' num2str(round(3*BW)) '~\mathrm{Hz}$'];

strP = ['$\mathrm{min}\epsilon_p = ' num2str(round(Ep)) '~\mathrm{Hz}$'];

f6 = figure(25);

clf
f6.Color='w';
f6.Position=[50 50 600 300];
pS=plot(Enorm*E_fermi,myP,'.');
hold on
xlim([0 2*E_fermi])
xlabel('energy (Hz)')
ylabel('occupation');
set(gca,'fontsize',12,'fontname','times','box','on',...
    'linewidth',1);
pEF = plot([1 1]*myMu*E_fermi,[0 1],'-','linewidth',1,'color',...
    co(1,:)*.6);
pMu = plot([1 1]*myMu*E_fermi,[0 1],'b:','linewidth',1);

pBW = plot([1 1]*3*BW,[0 1],'k-','linewidth',1);

pP = plot([1 1]*Ep,[0 1],'r-','linewidth',1);


legend([pS, pEF, pMu pBW pP],{strS strEF strMu strBW strP},...
    'interpreter','latex','location','southwest');

%%

keyboard



end



