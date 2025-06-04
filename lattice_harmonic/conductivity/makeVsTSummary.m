
%% Vesus temperature figures

hF_vs_T=figure(3);
hF_vs_T.Color='w';
clf
hF_vs_T.Name='versus_temp_bootstrap';
co=get(gca,'colororder');

subplot(121);
errorbar(temp(:,1),gamma(:,1),gamma(:,2),gamma(:,2),temp(:,2),temp(:,2),...
    'o','color',co(1,:),'markerfacecolor',co(1,:))

% errorbar(Tx(:,1),gamma(:,1),gamma(:,2),gamma(:,2),Tx(:,2),Tx(:,2),...
    % 'o','color',co(1,:),'markerfacecolor',co(1,:))


xlabel('spectral T/t [t]')
ylabel('\Gamma [1/s]')
set(gca,'fontsize',14,'fontname','times')
title('current dissipation');
xlim([0 4])


subplot(122);
p0=errorbar(temp(:,1),rho0(:,1),rho0(:,2),rho0(:,2),temp(:,2),temp(:,2),...
    'o','color',.5*co(3,:),'markerfacecolor',co(3,:));
% errorbar(Tx(:,1),rho0(:,1),rho0(:,2),rho0(:,2),Tx(:,2),Tx(:,2),...
    % 'o','color',co(1,:),'markerfacecolor',co(1,:))
hold on
% pinf=errorbar(temp(:,1),rhoinf(:,1),rhoinf(:,2),rhoinf(:,2),temp(:,2),temp(:,2),...
    % 'o','color',.5*co(4,:),'markerfacecolor',co(4,:));
xlabel('spectral T/t [t]')
ylabel('\rho [1/\sigma_0]')
% legend([p0 pinf],{'$\rho(\mathrm{Im}[\sigma]=0)$','$\rho(\omega\rightarrow \infty)$'},'interpreter','latex',...
%     'location','northwest')

set(gca,'fontsize',14,'fontname','times')
title('resistivity');
xlim([0 4])