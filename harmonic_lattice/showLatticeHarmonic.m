function f=showLatticeHarmonic(input,npt)

for uu=1:length(input.Depth)
    f = figure;
    f.Color='w';
    f.Position=[10 550 500 400];
    co = get(gca,'colororder');
    inds = 1:size(input.EigenValues(:,uu));
    E_min = min(input.EigenValues(:,uu));

    nstates = input.NumSites*input.NumBands;
    eng = input.EigenValues(:,uu);

    legStr={};
    clear plist
    plist=[];

    strLabel = ['$U_0:' num2str(input.Depth) 'E_R' ...
        ',\Omega:' num2str(input.Omega,2) '~\mathrm{Hz}' ...
        ',\omega:2\pi\cdot' num2str(input.omega/(2*pi),2) '~\mathrm{Hz}' ...
        ',E_R:h\cdot' num2str(input.fr*1e-3,4) '~\mathrm{kHz}$'];
    strLabel = [strLabel newline ...
        num2str(input.NumSites) ' sites' ...
        ', up to ' num2str(input.MaxTunnelingOrder) ' sites tunneling' ...
        ', ' num2str(input.NumBands) ' bands'];

    text(.01,.98,strLabel,'units','normalized','verticalalignment','top',...
        'horizontalalignment','left','units','normalized','interpreter','latex');


    % Band Patches
    for nn = 1:input.NumBands
        % plot([1 nstates],[1 1]*min(npt.bandEigenValue(nn,:))*npt.fr,'-','color',co(mod(nn-1,7)+1,:));
        % plot([1 nstates],[1 1]*max(npt.bandEigenValue(nn,:))*npt.fr,'-','color',co(mod(nn-1,7)+1,:));
        e1=npt.fr*min(npt.bandEigenValue(nn,:));
        e2=npt.fr*max(npt.bandEigenValue(nn,:));
        i1=find(eng>=e1,1);
        i2=find(eng>=e2,1);        
        if ~isempty(i1) && ~isempty(i2)
            x2 = [i1 i1 i2 i2];
            y2 = [e1 e2 e2 e1];
           pp(nn)=patch(x2,y2,'r');
           set(pp(nn),'facecolor',co(mod(nn-1,7)+1,:),'FaceAlpha',.05,'linestyle',':','linewidth',.5);
           hold on
           legStr{end+1} =['band ' num2str(nn) ' bare'];
           plist(end+1)=pp(nn);
        end
    end


    % Quantum Harmonic Oscillator Approximation for each Band
    % uses effective mass (not sure what to do for negative masses)

    % Harmonic energy approximation
    i0 = find(eng>=0,1);
    foo_ho = @(ind) (ind-i0)*input.omega/(2*pi);
    pHO=plot([1 nstates],foo_ho([1 nstates]),'-','linewidth',2,...
        'color',[.3 .3 .3]);
    hold on
    strHO = ['$nh\times' num2str(input.omega/(2*pi)) '\mathrm{Hz}$'];
    legStr{end+1}=strHO;
    plist(end+1)=pHO;

    if isfield(input,'FirstBandLinearSlope')
        m=input.FirstBandLinearSlope;
        nmax = find(input.EigenValues>input.BandRanges(1,2),1);
        nmax = nmax-1;
        pHOFit=plot([1 nmax],[0 (nmax-1)]*m+min(input.EigenValues),':','linewidth',2,...
            'color',[.3 .3 .3]);
        strHOFit = ['$nh\times' num2str(round(m,1)) '\mathrm{Hz}$'];
        legStr{end+1}=strHOFit;
        plist(end+1)=pHOFit;
    end
    
    % Actual Energy Data
    [~,dominateBandIndex] = max(input.BandProjection(:,:,uu),[],2);
    colors = co(mod(dominateBandIndex-1,7)+1,:);
    pData=scatter(inds,eng,2,colors,'linewidth',2);

    % Kinetic Energy
    pKE = plot(ones(input.NumBands,1),input.fr*(1:input.NumBands).^2+E_min,'ko',...
        'markerfacecolor','k','markersize',5);
    strKE = '$\epsilon_\mathrm{KE} = E_R s^2+\mathrm{min}(\epsilon)$';
    legStr{end+1}=strKE;
    plist(end+1)=pKE;

    %
    e0 = npt.fr*min(npt.bandEigenValue(1,:));
    e1 = npt.fr*max(npt.bandEigenValue(1,:));
    i1 = find(eng>=e0,1);
    i1=1;
    i2 = find(eng>=e1,1);  

    if ~isempty(i1) && ~isempty(i2)  && ~isequal(i1,i2)
        omega_band = sqrt(npt.BandCurvatureG(1))*input.omega;  
        f_band = omega_band/(2*pi);
        
%         foo_ho_1_band = @(ind) (ind-1)*f_band+0.5*input.omega/(2*pi); 
%         pBandHO=plot([i1:1:i2],foo_ho_1_band([i1:1:i2])+e0,'.-','linewidth',1,...
%             'color',[.3 .3 .3]);
        
        foo_ho_1_band = @(ind) (ind-1)*f_band+min(input.EigenValues);
        pBandHO=plot([i1:1:i2],foo_ho_1_band([i1:1:i2]),'.-','linewidth',1,...
            'color',[.3 .3 .3]);
        hold on
        strHOLattice = ['$n \hbar \omega\sqrt{m_\mathrm{eff}/m}=n h \times ' ...
            num2str(round(f_band,1)) '~\mathrm{Hz}$'];
        strHOLattice = ['$n h \times ' num2str(round(f_band,1)) '~\mathrm{Hz}$'];
        legStr{end+1} = strHOLattice;
        plist(end+1)=pBandHO;
        
        
    end

    % Labels and Limits
    xlabel('eigen index');
    ylabel('energy (Hz)');
    set(gca,'box','on','linewidth',1,'fontsize',12,'fontname','times');
    Emax = input.fr*input.NumBands^2+E_min+2*input.fr;
    iL = find(eng>=Emax,1);
    if ~isempty(iL)
        xlim([1 iL]);
    else
        xlim([1 nstates]);
    end
    ylim([E_min Emax]);

 

    % legend([pHO pKE pBandHO pp],legStr,'interpreter','latex','location','southeast',...
        % 'fontsize',10)

    legend(plist,legStr,'interpreter','latex','location','southeast',...
        'fontsize',10)

end


end

