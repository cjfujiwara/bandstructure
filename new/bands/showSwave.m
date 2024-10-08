function [hF] = showSwave(npt,opts)

%% Interaction Saturation
[busch_harmonic_energy] = trapped_interaction;
aVec=linspace(-2,4,1000);


hF=figure;
clf
% plot(aVec,funcs{1}(aVec))
set(gcf,'color','w');
hF.Position=[50 50 800 400];
hF.Name='saturation_harmonic_interaction';

U_Wannier = npt.OverlapS.U_Wannier_Er(1,:);
co=get(gca,'colororder');


% U_buchler = @(W,w02,G0,X0) W*w02./(

sWannier = @(U) ['$g\int |w_0|^4\mathrm{d}^3x$ ' num2str(U) 'Er'];

myc=parula(length(npt.depth));

clear ps
legStr={};

for nn=1:length(U_Wannier)
    cc = myc(nn,:);
    U_HO = U_Wannier(nn)/(sqrt(4*npt.OverlapS.depth(nn)));
%     pLinearWannier(nn)=plot(aVec,aVec*U_HO,'-','linewidth',1,...
%         'color',cc);
    ps(nn)=plot(aVec,aVec*U_HO,'-','linewidth',1,...
        'color',co(nn,:));
    legStr{nn}=['wannier ' num2str(npt.OverlapS.depth(nn)) ' Er'];
    hold on
end


pLinearHO=plot(aVec,aVec*sqrt(2/pi),'k--','linewidth',1);
hold on

pBusch=plot(aVec,busch_harmonic_energy{1}(aVec),'k-','linewidth',1);
pBusch=plot(aVec,busch_harmonic_energy{3}(aVec),'k-','linewidth',1);
pBusch=plot(aVec,busch_harmonic_energy{2}(aVec),'k-','linewidth',1);

% for nn=1:length(U_Wannier)
%     da = 0.01;
%     E_ho_small = busch_harmonic_energy{2}(da);
%     
%     U_wannier_HO = U_Wannier(nn)/(sqrt(4*npt.OverlapS.depth(nn)));
%     E_wannier_small = U_wannier_HO*da;
%     
%     R = E_wannier_small/E_ho_small;
%  
%     pBuschWannier=plot(aVec,R*busch_harmonic_energy{2}(aVec),'-','linewidth',1,'color',co(nn,:));
% 
% end

hold on
% legStr={'bound','first unbound','linear'};
% legend([ps(1:2) pFree],legStr,'location','southwest');



 out=buchler_interaction;

ps(end+1)=plot(aVec,out.U_04(aVec),'--','linewidth',1,'color',co(1,:));legStr{end+1}='buchler 4Er';
ps(end+1)=plot(aVec,out.U_08(aVec),'--','linewidth',1,'color',co(2,:));legStr{end+1}='buchler 8Er';
ps(end+1)=plot(aVec,out.U_12(aVec),'--','linewidth',1,'color',co(3,:));legStr{end+1}='buchler 12Er';
ps(end+1)=plot(aVec,out.U_16(aVec),'--','linewidth',1,'color',co(4,:));legStr{end+1}='buchler 16Er';
ps(end+1)=plot(aVec,out.U_20(aVec),'--','linewidth',1,'color',co(5,:));legStr{end+1}='buchler 20Er';


ylim([0 1]);
xlim([0 1.5]);

set(gca,'xgrid','on','ygrid','on','fontsize',10,'box','on',...
    'linewidth',1);

xlabel('scattering length $(a_\mathrm{HO})$','interpreter','latex');
ylabel('energy ($\hbar \omega_{\mathrm{ho}}$)','interpreter','latex');

legend([pLinearHO pBusch ps],...
    [{'$\sqrt{2/\pi}a_s/a_\mathrm{ho}$', 'Busch'},legStr],'location','eastoutside','interpreter','latex',...
    'fontsize',10);

end

