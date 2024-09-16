function [hF_W_psi] = showWannier(npt,opts)

    bands = opts.bands;

for uu=1:length(npt.depth)
    tstart=now;
    fprintf(['plotting wannier ' num2str(uu) '/' num2str(length(npt.depth)) ' ...']);

    U0 = npt.depth(uu);
    %% Wavefunction Plot
    hF_W_psi(uu)=figure(10000+uu);
    hF_W_psi(uu).Color='w';
    hF_W_psi(uu).Position=[260 50 1200 400];
    hF_W_psi(uu).Name = 'wannier_wavefunction';
    clf
    % t=uicontrol('style','text','string',[num2str(depth) 'Er'],...
    %     'fontname','times','fontsize',12,'backgroundcolor','w');
    % t.Position(3:4)=t.Extent(3:4);
    % t.Position(1:2)=[5 hF_W_psi.Position(4)-t.Extent(4)];
    for nn=1:length(bands)
        subplot(2,length(bands),nn)    
        p1=plot(npt.K_extended,real(npt.Wannier_K(:,nn,uu)),'-','linewidth',1);
        hold on
        p2=plot(npt.K_extended,imag(npt.Wannier_K(:,nn,uu)),'-','linewidth',1);    
        xlabel('momentum ($\hbar k_L$)','interpreter','latex');
        str=['$' num2str(U0) 'E_R$' newline '$w_' num2str(bands(nn)) '(k)$'];
        text(.01,.98,str,'units','normalized','interpreter','latex',...
            'verticalalignment','top','fontsize',12);
        set(gca,'xgrid','on','ygrid','on','fontname','times',...     
            'fontsize',10,'yticklabel',{});
        xlim([-5 5]);

        % yyaxis right
        % p3=plot(npt.K_extended,abs(npt.Wannier_K(:,nn,uu)).^2,'k-','linewidth',1);    
        % set(gca,'yColor','k');
        % legend([p1 p2 p3],{'Re','Im','norm^2'});

            % title(['$U_0:' num2str(U0) 'E_R,~w_' num2str('n') '(k)$'],'interpreter','latex');

        subplot(2,length(bands),length(bands)+nn)   
        plot(npt.X_extended,real(npt.Wannier_X(:,nn,uu)),'-','linewidth',1);
        hold on
        plot(npt.X_extended,imag(npt.Wannier_X(:,nn,uu)),'-','linewidth',1);    
        xlabel('position (site)','interpreter','latex');
        str=['$' num2str(U0) 'E_R$' newline '$w_' num2str(bands(nn)) '(x)$'];
        text(.01,.98,str,'units','normalized','interpreter','latex',...
            'verticalalignment','top','fontsize',10);
        set(gca,'xgrid','on','ygrid','on','fontname','times',...     
            'fontsize',10,'yticklabel',{});
        % xlim([-1 1]*10*npt.Harmonic_Length(uu))
        xlim([-5 5]);
        hold on    
        plot(npt.X_extended,npt.Wannier_X_Harmonic(:,nn,uu),'-','linewidth',2,'color',[.3 .3 .3 .7]);
    end
    tend = toc;
    disp([' done (' num2str(round(tend,2)) ' s)']);
    drawnow;
end

end

