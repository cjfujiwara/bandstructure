function hFw = showTrappedInteraction
%% Interaction Saturation
[funcs] = trapped_interaction;
aVec=linspace(-2,4,1000);

hFw=figure;
clf
% plot(aVec,funcs{1}(aVec))
set(gcf,'color','w');
hFw.Position=[1150 50 250 250];
hFw.Name='saturation_harmonic_interaction';

pFree=plot(aVec,aVec*sqrt(2/pi),'k--','linewidth',2);
hold on
co=get(gca,'colororder');
for kk=1:4
    ps(kk)=plot(aVec,funcs{kk}(aVec),'-','linewidth',2,...
        'color',co(kk,:));
    hold on
end

% legStr={'bound','first unbound','linear'};
% legend([ps(1:2) pFree],legStr,'location','southwest');

legend([pFree],'$\sqrt{2/\pi}$','location','southwest','interpreter','latex');

ylim([-3 3]);
xlim([-2 2]);

set(gca,'xgrid','on','ygrid','on','fontsize',10,'box','on',...
    'linewidth',1);

xlabel('scattering length $(a_\mathrm{HO})$','interpreter','latex');
ylabel('binding energy ($\hbar \omega$)','interpreter','latex');

end

