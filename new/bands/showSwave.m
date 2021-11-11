function [hF] = showSwave(npt)

%% Interaction Saturation
[funcs] = trapped_interaction;
aVec=linspace(-2,4,1000);

hF=figure;
clf
% plot(aVec,funcs{1}(aVec))
set(gcf,'color','w');
hF.Position=[1150 50 600 400];
hF.Name='saturation_harmonic_interaction';

U_Wannier = npt.OverlapS.U_Wannier_Er(1,:);
co=get(gca,'colororder');




myc=parula(length(npt.depth));
for nn=1:length(U_Wannier)
    cc = myc(nn,:);
    U_HO = U_Wannier(nn)/(sqrt(4*npt.OverlapS.depth(nn)));
    pLinearWannier(nn)=plot(aVec,aVec*U_HO,'-','linewidth',1,...
        'color',cc);
    hold on
end

pLinearHO=plot(aVec,aVec*sqrt(2/pi),'k-','linewidth',2);
hold on

p=plot(aVec,funcs{2}(aVec),'k:','linewidth',2);
hold on
% legStr={'bound','first unbound','linear'};
% legend([ps(1:2) pFree],legStr,'location','southwest');

legend([pLinearHO p pLinearWannier(1) pLinearWannier(end)],...
    {'$\sqrt{2/\pi}$ HO', 'sat. HO',['lin. w ' num2str(npt.depth(1)) '$E_R$' ],['lin. w ' num2str(npt.depth(end)) '$E_R$' ]},'location','northwest','interpreter','latex',...
    'fontsize',10);

ylim([0 1]);
xlim([0 1.5]);

set(gca,'xgrid','on','ygrid','on','fontsize',10,'box','on',...
    'linewidth',1);

xlabel('scattering length $(a_\mathrm{HO})$','interpreter','latex');
ylabel('binding energy ($\hbar \omega$)','interpreter','latex');

cb=colorbar;
cb.Label.String = 'depth ($E_R$)';
cb.Label.Interpreter='latex';
colormap(parula)
caxis([min(npt.depth) max(npt.depth)]);
end

