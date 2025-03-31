function hFw = showTrappedInteraction
%% Interaction Saturation
[funcs] = trapped_interaction;
aVec=linspace(-10,10,1000);

hFw=figure;
clf
% plot(aVec,funcs{1}(aVec))
set(gcf,'color','w');
hFw.Position=[1150 50 400 400];
hFw.Name='saturation_harmonic_interaction';

pFree=plot(aVec,aVec*sqrt(2/pi),'k--','linewidth',2);
hold on
co=get(gca,'colororder');
for kk=1:5
    ps(kk)=plot(aVec,funcs{kk}(aVec),'-','linewidth',2,...
        'color',co(kk,:));
    ps(kk).Color='k';
    hold on
    
    plot([-10 10],2*[1 1]*(kk-1.5),'k:');
end

% legStr={'bound','first unbound','linear'};
% legend([ps(1:2) pFree],legStr,'location','southwest');

% legend([pFree],'$\sqrt{2/\pi}$','location','southwest','interpreter','latex');
legend([pFree],'linear','location','southwest','interpreter','latex');

ylim([-3 3]);
xlim([-8 8]);

set(gca,'xgrid','on','ygrid','on','fontsize',14,'box','on',...
    'linewidth',1);

xlabel('scattering length $(a_\mathrm{HO})$','interpreter','latex');
ylabel('energy ($\hbar \omega$)','interpreter','latex');
title('3d hamornic + contact s-wave')

end

