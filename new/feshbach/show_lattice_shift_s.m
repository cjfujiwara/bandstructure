function hF = show_lattice_shift_s

U=linspace(1,500,1000);

hF=figure(512);
set(hF,'color','w','name','lattice_shift_s');
hF.Position=[0 710 575 250];
clf


[dB1,~]=lattice_shift_79(U,U,U);
[dB2,~]=lattice_shift_95(U,U,U);



% subplot(121);
p1=plot(U,1e3*dB1,'-','LineWidth',2); 
hold on
p2=plot(U,1e3*dB2,'-','LineWidth',2); 

set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on');
xlabel('lattice depth (E_R)');
ylabel('resonance shift (mG)');

str='$\Delta B = \Delta/\left(1-\frac{2\hbar^2}{2\mu R \hbar \bar{\omega} a_\mathrm{bg}}\right)$';
text(.98,.02,str,'units','normalized','fontsize',10,'interpreter','latex',...
    'verticalalignment','bottom','horizontalalignment','right');

legend([p1 p2],{'79 ','95'},'location','northwest');

% 
% subplot(122);
% plot(U,B,'-','LineWidth',2); 
% set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on');
% xlabel('lattice depth (E_R)');
% ylabel('resonance (G)');

drawnow;



set(gca,'fontsize',10,'xgrid','on','ygrid','on','box','on','linewidth',1);
end

