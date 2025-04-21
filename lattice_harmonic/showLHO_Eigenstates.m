function ff = showLHO_Eigenstates(lattice,npt)

t1=npt.Tunneling(1,1)*npt.fr;
t2=npt.Tunneling(1,2)*npt.fr;

ustr = ['V_0=' num2str(npt.Depth) '~E_\mathrm{R}'];
omegastr=['\omega=2 \pi \cdot' num2str(npt.omega/(2*pi)) '~\mathrm{Hz}'];
title_str = ['$' ustr ';' omegastr '$'];
tstr = ['$t=(' num2str(round(t1)) ',' num2str(round(t2,1)) ',...)~\mathrm{Hz}$'];

%% Init Figure Visualize

ff=figure;
ff.Color='w';
ff.Position=[50 50 1200 700];
clf

%% Eigenvalues

ax2=subplot(2,3,1,'parent',ff);
co=get(gca,'colororder');
eng = npt.EigenValues;
e0=min(eng);
legStr={};
clear plist
plist=[];

% Band Patch - make a patch for each band
nn=1;
e1=npt.fr*min(lattice.bandEigenValue(nn,:));
e2=npt.fr*max(lattice.bandEigenValue(nn,:));
i1=find(eng>=e1,1);
i2=find(eng>=e2,1);        
if ~isempty(i1) && ~isempty(i2)
    x2 = [i1 i1 i2 i2];
    y2 = [e1 e2 e2 e1];
   pp(nn)=patch(x2,y2-e0,'r');
   set(pp(nn),'facecolor',co(mod(nn-1,7)+1,:),'FaceAlpha',.05,'linestyle',':','linewidth',.5);
   hold on
   legStr{end+1} =['$s$-band'];
   plist(end+1)=pp(nn);
end

% Data
plist(end+1)=plot(1:npt.NumSites,eng-e0,'o','markeredgecolor','k',...
    'markerfacecolor',co(1,:));
legStr{end+1}='data';
xlabel('eigenindex');
ylabel('energry (Hz)');

% Limits
nMax = x2(end) + 5;
xlim([1 nMax]);

% Effective Mass
% Energy spectrum assuming k=0 effective mass
mstar=lattice.BandMassGamma(1);
plist(end+1)=plot(i1:i2,((i1:i2)-1)*(npt.omega/(2*pi))/sqrt(mstar),'.-',...
    'color',co(2,:),'linewidth',1);
legStr{end+1}=['$n \hbar\omega\sqrt{m/m^*}\approx nh\cdot' ...
    num2str(round(npt.omega/(2*pi*sqrt(mstar)))) '~\mathrm{Hz}$'];

% Best Fit Line
% if isfield(npt,'FirstBandLinearSlope')
%     m=npt.FirstBandLinearSlope;
%     nmax = find(npt.EigenValues>npt.BandRanges(1,2),1);
%     pHOFit=plot([1 nmax],[0 (nmax-1)]*m,'-','linewidth',1,...
%         'color','k');
%     strHOFit = ['$\mathrm{fit:}~nh\times' num2str(round(m,1)) '\mathrm{Hz}$'];
%     legStr{end+1}=strHOFit;
%     plist(end+1)=pHOFit;
% end
    
yL = ax2.YLim;
set(ax2,'YLim',[0 yL(2)])

legend(plist,legStr,'interpreter','latex','location','northwest',...
        'fontsize',8)
    
set(gca,'box','on','linewidth',1);
title('eigenvalues','interpreter','latex')

%% Eigenvector Visualize

ax1=subplot(2,3,4);
imagesc(1:npt.NumSites,npt.PositionVector,abs(npt.EigenVectors).^2);
set(gca,'YDir','normal')
ylabel('site');
xlabel('eigenindex');
title(title_str,'interpreter','latex')
caxis([0 .05]);
ylim([-60 60])
xlim([1 nMax]);

%% Neighboring Eigenstates Energy
ax3=subplot(2,3,2);

Nvec = 1:npt.NumSites;
Y    = eng-e0;
dY   = diff(Y,1);


plot(Nvec(1:(end-1)),dY,'o','markeredgecolor','k',...
    'markerfacecolor',co(1,:));
xlim([1 nMax])
ylim([0 dY(1)*1.1])
hold on
xlabel('eigenindex $i$','interpreter','latex')
title('energy of $\Delta i = 1$','interpreter','latex')
ylabel('energy difference [Hz]','interpreter','latex')

%% Neighboring Eigentstes Energy
ax4=subplot(2,3,3);

Ysub = Y(4:end)-Y(1:end-3);
plot(Nvec(1:(end-3)),Ysub,'o','markeredgecolor','k',...
    'markerfacecolor',co(1,:));
xlim([1 nMax])
ylim([0 Ysub(1)*1.1])

hold on
xlabel('eigenindex $i$','interpreter','latex')
title('energy of $\Delta i = 3$','interpreter','latex')
ylabel('energy difference [Hz]','interpreter','latex')

%% Dipole Operator
ax5=subplot(2,3,5);

D = npt.DipoleOperator;
D1 = diag(D,1);
D3 = diag(D,3);


plot(Nvec(1:(end-1)),abs(D1),'o','markeredgecolor','k',...
    'markerfacecolor',co(1,:));
xlim([1 nMax])
xlabel('eigenindex $i$','interpreter','latex')
title('dipole of $\Delta i = 1$','interpreter','latex')
ylabel('$|D_{i+1,i}|~[a_L]$','interpreter','latex')

%% Dipole Operator
ax6=subplot(2,3,6);

plot(Nvec(1:(end-3)),abs(D3),'o','markeredgecolor','k',...
    'markerfacecolor',co(1,:));
ylabel('$|D_{i+3,i}|~[a_L]$','interpreter','latex')
xlabel('eigenindex $i$','interpreter','latex')
xlim([1 nMax])
title('dipole of $\Delta i = 3$','interpreter','latex')


%% FInish Up
linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'x');

end

