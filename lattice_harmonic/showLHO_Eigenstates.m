function f = showLHO_Eigenstates(lattice,npt)

t1=npt.Tunneling(1,1)*npt.fr;
t2=npt.Tunneling(1,2)*npt.fr;


ustr = ['V_0=' num2str(npt.Depth) '~E_\mathrm{R}'];

omegastr=['\omega=2 \pi \cdot' num2str(npt.omega/(2*pi)) '~\mathrm{Hz}'];

title_str = ['$' ustr ';' omegastr '$'];


tstr = ['t=(' num2str(round(t1)) ',' num2str(round(t2,1)) ',...)~\mathrm{Hz}'];


ff=figure;
ff.Color='w';
ff.Position=[50 50 600 600];
clf

ax1=subplot(2,1,1);
imagesc(1:npt.NumSites,npt.PositionVector,abs(npt.EigenVectors).^2);
set(gca,'YDir','normal')
ylabel('site');
xlabel('eigenindex');
title(title_str,'interpreter','latex')
caxis([0 .05]);
ylim([-60 60])
%%
legStr={};
clear plist
    plist=[];

ax2=subplot(2,1,2);
co=get(gca,'colororder');

eng = npt.EigenValues;
e0=min(eng);

% Band Patch
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
linkaxes([ax1 ax2],'x');
xlim([1 65]);

% Effective Mass
mstar=lattice.BandMassGamma(1);
plist(end+1)=plot(i1:i2,((i1:i2)-1)*(npt.omega/(2*pi))/sqrt(mstar),'.-',...
    'color','r');
legStr{end+1}=['$n \hbar\omega\sqrt{m/m^*}\approx nh\cdot' ...
    num2str(round(npt.omega/(2*pi*sqrt(mstar)))) '~\mathrm{Hz}$'];
% keyboard


if isfield(npt,'FirstBandLinearSlope')
    m=npt.FirstBandLinearSlope;
    nmax = find(npt.EigenValues>npt.BandRanges(1,2),1);
%     nmax = nmax-1;
    pHOFit=plot([1 nmax],[0 (nmax-1)]*m,'-','linewidth',2,...
        'color','k');
    strHOFit = ['$nh\times' num2str(round(m,1)) '\mathrm{Hz}$'];
    legStr{end+1}=strHOFit;
    plist(end+1)=pHOFit;
end
    

 legend(plist,legStr,'interpreter','latex','location','southeast',...
        'fontsize',10)
    
set(gca,'box','on','linewidth',1);
end

