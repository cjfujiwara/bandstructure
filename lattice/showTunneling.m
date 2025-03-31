function showTunneling(npt,opts)

hF=figure;
hF.Color='w';
clf

if ~isfield(opts,'bands')
opts.bands = 1;
end
co=get(gca,'colororder');

U = [npt.depth];
U=U(:);
line_spec={'-','--','-.','-.','-.','-.','-.'};
for band_index=1:length(opts.bands)
        subplot(1,length(opts.bands),band_index);
        % yyaxis left
        % set(gca,'YColor','k');
        % yLeft=gca;

    for site_index=1:5
        band_number = opts.bands(band_index);
        site_number = site_index;
        T = npt.Tunneling(band_number,site_number,:);
        T=T(:);


        plot(U,T*npt.fr*1e-3,line_spec{site_index},'linewidth',1,...
            'color',co(band_index,:));
        hold on;

    end
        % ylabel('tunneling (Er)')

    % yL=get(gca,'YLim');
    % yyaxis right
    % % yRight=gca;
    % set(gca,'YColor','k');
    % ylim(yL*npt.fr*1e-3)
    ylabel('tunneling (kHz)')
    xlabel('lattice depth (Er)');
 
end

end

