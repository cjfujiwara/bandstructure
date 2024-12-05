function hF = showWannierFunction(npt)
  
fprintf(' plotting ... ');
tstart=now;

hF_W_probability=figure(10000+ii);
hF_W_probability.Color='w';
hF_W_probability.Position=[250 50 800 400];
hF_W_probability.Name = 'wannier_probability';
clf

t=uicontrol('style','text','string',[num2str(depth) 'Er'],...
    'fontname','times','fontsize',12,'backgroundcolor','w');
t.Position(3:4)=t.Extent(3:4);
t.Position(1:2)=[5 hF_W_probability.Position(4)-t.Extent(4)];
for nn=1:length(bands)
    subplot(2,length(bands),nn)
    psi_wannier_k = wannier_k(:,nn);
    rho_wannier_k = psi_wannier_k.*conj(psi_wannier_k);
    rho_wannier_k = rho_wannier_k/max(rho_wannier_k);

    plot(K_eff,rho_wannier_k,'-','linewidth',1);
    xlabel('momentum ($\hbar k_L$)','interpreter','latex');
    str=['$|w_' num2str(bands(nn)) '(k)|^2$'];
    text(.01,.98,str,'units','normalized','interpreter','latex',...
        'verticalalignment','top','fontsize',10);
    set(gca,'xgrid','on','ygrid','on','fontname','times',...     
        'fontsize',8,'yticklabel',{});
    xlim([-20 20]);

    subplot(2,length(bands),length(bands)+nn)
    psi_wannier = wannier_x(:,nn);
    rho_wannier = psi_wannier.*conj(psi_wannier);
    rho_wannier = rho_wannier/max(rho_wannier);

    plot(X/pi,rho_wannier,'-','linewidth',1);
    xlabel('position (site)','interpreter','latex');
    str=['$|w_' num2str(bands(nn)) '(x)|^2$'];

    text(.01,.98,str,'units','normalized','interpreter','latex',...
        'verticalalignment','top','fontsize',10);
    set(gca,'xgrid','on','ygrid','on','fontname','times',...     
        'fontsize',8,'yticklabel',{});
    xlim([-.5 .5]);
    hold on    

    psi_harmonic_x = wannier_x_harmonic(:,nn);
    rho_ho = psi_harmonic_x.*conj(psi_harmonic_x);
    rho_ho = rho_ho/max(rho_ho);
%     plot(X/pi,rho_ho,'k--','linewidth',1);    
    plot(X/pi,rho_ho,'-','linewidth',2,'color',[.3 .3 .3 .7]);  

end

%% Wavefunction Plot
hF_W_psi=figure(10001);
hF_W_psi.Color='w';
hF_W_psi.Position=[1055 50 800 400];
hF_W_psi.Name = 'wannier_wavefunction';
clf
t=uicontrol('style','text','string',[num2str(depth) 'Er'],...
    'fontname','times','fontsize',12,'backgroundcolor','w');
t.Position(3:4)=t.Extent(3:4);
t.Position(1:2)=[5 hF_W_psi.Position(4)-t.Extent(4)];
for nn=1:length(bands)
    subplot(2,length(bands),nn)
    psi_wannier_k = wannier_k(:,nn);

    plot(K_eff,real(psi_wannier_k),'-','linewidth',1);
    hold on
    plot(K_eff,imag(psi_wannier_k),'-','linewidth',1);

    xlabel('momentum ($\hbar k_L$)','interpreter','latex');
    str=['$w_' num2str(bands(nn)) '(k)$'];
    text(.01,.98,str,'units','normalized','interpreter','latex',...
        'verticalalignment','top','fontsize',10);
    set(gca,'xgrid','on','ygrid','on','fontname','times',...     
        'fontsize',8,'yticklabel',{});
    xlim([-20 20]);

    subplot(2,length(bands),length(bands)+nn)
    psi_wannier_x = wannier_x(:,nn);


    plot(X/pi,real(psi_wannier_x),'-','linewidth',1);
    hold on
    plot(X/pi,imag(psi_wannier_x),'-','linewidth',1);

    xlabel('position (site)','interpreter','latex');
    str=['$w_' num2str(bands(nn)) '(x)$'];

    text(.01,.98,str,'units','normalized','interpreter','latex',...
        'verticalalignment','top','fontsize',10);
    set(gca,'xgrid','on','ygrid','on','fontname','times',...     
        'fontsize',8,'yticklabel',{});
    xlim([-.5 .5]);
    hold on    

    psi_harmonic_x = wannier_x_harmonic(:,nn);
    plot(X/pi,psi_harmonic_x,'-','linewidth',2,'color',[.3 .3 .3 .7]);  
end

tend = now;
fprintf([' done (' num2str(round((tend-tstart)*24*60*60,3)) ' s)']);


end

