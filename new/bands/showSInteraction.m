function [hF] = showSInteraction(npt,opts)

%% Interaction Saturation
aVec=linspace(-2,4,1e4);% ion unites of lattice site


hF=figure;
clf
% plot(aVec,funcs{1}(aVec))
set(gcf,'color','w');
hF.Position=[50 50 800 400];
hF.Name='saturation_harmonic_interaction';

U_Wannier = npt.OverlapS.U_Wannier_Er(1,:);
co=get(gca,'colororder');


myc=parula(length(npt.depth));

clear ps
legStr={};

for nn=1:length(U_Wannier)
    cc = myc(nn,:);
    
    aHO = aVec./npt.Harmonic_Length(nn);
    ps(nn)=plot(aVec,U_Wannier(nn).*aHO,'-','linewidth',1,...
        'color',co(nn,:));
    legStr{nn}=['wannier ' num2str(npt.OverlapS.depth(nn)) ' Er'];
    hold on
end

% pLinearHO=plot(aVec,aVec*sqrt(2/pi),'k--','linewidth',1);
hold on

% legStr={'bound','first unbound','linear'};
% legend([ps(1:2) pFree],legStr,'location','southwest');

 [outHO,outEr]=buchler_interaction;

ps(end+1)=plot(aVec,outEr.U_04(aVec),'--','linewidth',1,'color',co(1,:));legStr{end+1}='buchler 4Er';
ps(end+1)=plot(aVec,outEr.U_08(aVec),'--','linewidth',1,'color',co(2,:));legStr{end+1}='buchler 8Er';
ps(end+1)=plot(aVec,outEr.U_12(aVec),'--','linewidth',1,'color',co(3,:));legStr{end+1}='buchler 12Er';
ps(end+1)=plot(aVec,outEr.U_16(aVec),'--','linewidth',1,'color',co(4,:));legStr{end+1}='buchler 16Er';
ps(end+1)=plot(aVec,outEr.U_20(aVec),'--','linewidth',1,'color',co(5,:));legStr{end+1}='buchler 20Er';


ylim([0 3]);
xlim([0 .15]);

set(gca,'xgrid','on','ygrid','on','fontsize',10,'box','on',...
    'linewidth',1);

xlabel('scattering length $(a_L)$','interpreter','latex');
ylabel('energy ($E_R$)','interpreter','latex');

legend([ps],...
    [legStr],'location','eastoutside','interpreter','latex',...
    'fontsize',10);

end