function hF = showBandStructure(npt,opts)

for kk=1:length(npt.depth)   
    U0 = npt.depth(kk);
    eng = npt.bandEigenValue(:,:,kk);

    % Make the figure
    hF(kk)=figure(200+kk);
    set(hF(kk),'Name',['band_structure_' num2str(U0) 'Er'],'color','w');
    clf;
    hF(kk).Position=[0 50 250 400];

    % Initialize the axis
    ax1=axes;
    cla
    co=get(gca,'colororder');
    set(ax1,'fontsize',10,'box','on','linewidth',1,'fontname','times');
    xlabel('quasimomentum ($\hbar k_L$)','interpreter','latex');
    ylabel('energy ($E_R$)','interpreter','latex');
    xlim([min(npt.K) max(npt.K)]);

    U_low = floor(min(eng(min(opts.Bands),:))); 
    U_high = ceil(max(eng(max(opts.Bands),:))); 
    ylim([U_low U_high]+[-.5 .5]);


    hold on

    % Plot harmonic oscilator energies
    engHO=-U0+2*sqrt(U0)*(0.5+(0:10));
    for ii=1:length(engHO)
        plot([-1 1],[1 1]*engHO(ii),'k:','linewidth',1); 
    end

    % Plot the bands
    for nn=1:length(opts.Bands)
        ind = opts.Bands(nn);
       plot(npt.K,eng(ind,:),'-','linewidth',3,...
           'color',co(mod(ind-1,7)+1,:)); 
    end

    % Add depth label
    str=['$' num2str(U0) '~E_R$'];
    text(5,5,str,'interpreter','latex','units','pixels',...
        'verticalalignment','bottom','fontsize',16);

end
 
end

