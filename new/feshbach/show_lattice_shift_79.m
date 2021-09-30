function hF = show_lattice_shift_79

U=linspace(1,500,1000);

hF=figure(226);
set(hF,'color','w','name','lattice_shift_79');
hF.Position=[10 540 800 250];
clf


[dB,B]=lattice_shift_79(U,U,U);

subplot(121);
plot(U,1e3*dB,'-','LineWidth',2); 
set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on');
xlabel('lattice depth (E_R)');
ylabel('resonance shift (mG)');

str='$\Delta B = \Delta/\left(1-\frac{2\hbar^2}{2\mu R \hbar \bar{\omega} a_\mathrm{bg}}\right)$';
text(.98,.02,str,'units','normalized','fontsize',16,'interpreter','latex',...
    'verticalalignment','bottom','horizontalalignment','right');

subplot(122);
plot(U,B,'-','LineWidth',2); 
set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on');
xlabel('lattice depth (E_R)');
ylabel('resonance (G)');

drawnow;



set(gca,'fontsize',10,'xgrid','on','ygrid','on','box','on','linewidth',1);
end

