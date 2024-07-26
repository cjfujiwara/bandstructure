function showLatticeHarmonic(input,npt)

for uu=1:length(input.Depth)
    f = figure;
    f.Color='w';
    f.Position=[100 100 600 400];
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
    strHO = ['$\epsilon_\mathrm{HO} =nh\times' num2str(input.omega/(2*pi)) '\mathrm{Hz}$'];
    legStr{end+1}=strHO;
    plist(end+1)=pHO;

    if isfield(input,'FirstBandLinearSlope')
        m=input.FirstBandLinearSlope;
        nmax = find(input.EigenValues>input.BandRanges(1,2),1);
        nmax = nmax-1;
        pHOFit=plot([1 nmax],[0 (nmax-1)]*m+min(input.EigenValues),':','linewidth',2,...
            'color',[.3 .3 .3]);
        strHOFit = ['$\epsilon_\mathrm{HO} =nh\times' num2str(round(m,1)) '\mathrm{Hz}$'];
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
    i2 = find(eng>=e1,1);  

    if ~isempty(i1) && ~isempty(i2)  && ~isequal(i1,i2)
        omega_band = sqrt(npt.BandCurvatureG(1))*input.omega;
        foo_ho_1_band = @(ind) (ind-1)*omega_band/(2*pi)+0.5*input.omega/(2*pi);
        pBandHO=plot([i1 i2],foo_ho_1_band([i1 i2])+e0,'--','linewidth',2,...
            'color',[.3 .3 .3]);
        hold on
        strHOLattice = '$\epsilon_s \approx n \hbar \omega\sqrt{m_\mathrm{eff}/m} +\mathrm{min}(\epsilon)$';
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



% 
% p4 = plot((eng_4_bands),'-','linewidth',2,'color',co(3,:));
% hold on
% p3 = plot((eng_3_bands),'-','linewidth',2,'color',co(2,:));
% p2 = plot((eng_2_bands),'-','linewidth',2,'color',co(1,:));
% p1 = plot(eng,'k-','linewidth',2);
% 
% 
% E_min = min(eng_4_bands);
% 
% 
% xlim([1 800]);
% % xlim([1 1600]);
% 
% % ylim([0 5e4]);
% xlabel('eigenindex')
% ylabel('energy (Er)')
% ylim([E_min 12*fr]);
% 
% i1 = find(eng_4_bands>=(1*fr+0*E_min),1);
% i2 = find(eng_4_bands>=(2^2*fr+0*E_min),1);
% i3 = find(eng_4_bands>=(3^2*fr+0*E_min),1);
% 
% 
% plot([1 i1],[1 1]*(1+0*E_min)*fr,'--','color',co(1,:),'linewidth',2);
% plot([1 i2],[1 1]*(2^2+0*E_min)*fr,'--','color',co(2,:),'linewidth',2);
% plot([1 i3],[1 1]*(3^2+0*E_min)*fr,'--','color',co(3,:),'linewidth',2);
% 
% set(gca,'box','on','linewidth',1,'fontname','times','fontsize',14);
% legend([p1 p2 p3 p4],{'1 band','2 bands','3 bands','4 bands'}, ...
%    'location','southeast')
% 
% %%
% 
% Hstr = ['$\hat{H} = -\sum_{j} t^{(s)}_{j}\sum_i \hat{c}_{i+j}\hat{c}_i+\Omega\sum_i i^2$'];
% 
% tunneling_str_1 = ['$t^{(s)}=(' num2str(round(t(1,1),1)) ',' ...
%     num2str(round(t(1,2),1)) ',' ...
%     num2str(round(t(1,3),1)) ')~\mathrm{Hz}$'];
% tunneling_str_2 = ['$t^{(p)}=(' num2str(round(t(2,1),1)) ',' ...
%     num2str(round(t(2,2),1)) ',' ...
%     num2str(round(t(2,3),1)) ')~\mathrm{Hz}$'];
% omega_str = ['$\Omega=' num2str(round(Omega,2)) '~\mathrm{Hz}$'];
% 
% depth_str = ['$U_0=' num2str(U0) 'E_R$'];
% str_Ep = ['$\mathrm{min}(\epsilon_p) = ' num2str(round(E2+min(eng),1)) '~\mathrm{Hz}$'];
% 
% tunneling_str = [omega_str newline depth_str newline tunneling_str_1 newline tunneling_str_2];
% 
% 
% hF = figure(20);
% hF.Color='w';
% clf
% co=get(gca,'colororder');
% subplot(121);
% cla
% pD = plot(eng,'k-','linewidth',2);
% xlabel('eigenindex','interpreter','latex');
% ylabel('energy $\epsilon$ (Hz)','interpreter','latex')
% hold on
% 
% i1 = find(eng>=2*t(1,1),1);
% 
% if isempty(i1)
%     i1 = n;
% end
% 
% i2 = find(eng>=(E2+min(eng)),1);
% pBW = plot([1 i1],[1 1]*2*t(1,1),'b-','linewidth',.5);
% 
% pR = plot([1 1]*i1,2*t(1,1)*[.8 1.2],'b-','linewidth',.5);
% 
% text(i1,2*t1*1.2,num2str(i1),'interpreter','latex','units','data',...
%     'horizontalalignment','center','verticalalignment','bottom','fontsize',8);
% 
% if isempty(i2)
%     i2 = n;
% end
% 
% pNextBand = plot([1 i2],[1 1]*(E2+min(eng)),'-','color',co(3,:));
% 
% yH = min([15*t(1,1)+max(eng) (E2+min(eng)+2*t(1,1)) ]);
% % ylim([-2*t(1,1) yH])
% xlim([1 i2*1.1])
% 
% 
% text(.02,.98,Hstr,'interpreter','latex','units','normalized',...
%     'horizontalalignment','left','verticalalignment','top','fontsize',12);
% 
% text(.98,.02,tunneling_str,'interpreter','latex','units','normalized',...
%     'horizontalalignment','right','verticalalignment','bottom','fontsize',12);
% text(2,E2+min(eng),str_Ep,'interpreter','latex','units','data',...
%     'horizontalalignment','left','verticalalignment','cap','fontsize',12);
% text(2,2*t1,'$\mathrm{min}(\epsilon_s)+4t_{s1}$','interpreter','latex','units','data',...
%     'horizontalalignment','left','verticalalignment','bottom','fontsize',12);
% set(gca,'fontsize',14,'fontname','times','box','on');
% 
% 
%     subplot(122);
%     dN = 5;
% 
% for kk=1:dN:55
%     y = Csymm(:,kk);
%     y = y/max(abs(y));
%     y = y*dN*.4;    
%     % c = co(1+round(P(kk)),:);
%     c = co(1+mod(kk,2),:);
% 
%     plot(x,y+kk,'.','color',c,'linewidth',1);
%     % hold on
%         % plot(xsample,y2+kk,'color',c,'linewidth',1);
% 
%     str = ['$' num2str(kk) ':' num2str(round(eng(kk))) '~\mathrm{Hz}$'];
%     text(-58,kk,str,'verticalalignment','bottom','horizontalalignment','left',...
%         'fontsize',8,'interpreter','latex');
%     hold on    
% end
% xlim([-60 60]);
% set(gca,'fontsize',14,'fontname','times','box','on');
% xlabel('position (sites)');
% ylabel('wavefunction (arb.)');
% set(gca,'YTickLabel',{});
% 
% %% Lowest Lying Energies
% % WARNING DO NOT CHOOSE TOO LARGE, YOU WILL RUN OUT OF COMPUTER MEMORY
% %
% % n = 1000 --> n^3 combinations --> n^3*8bytes
% %
% % Also note that our imaging only allows for at most 200 sites
% % Also note that we only have around 50e3 atoms of single spin at most
% nmax = 51; % DO NOT MAKE LARGER THAN 301
% % nmax = 201; % DO NOT MAKE LARGER THAN 301
% 
% % nmax = n;e
% % 
% % eng = eng(1:2:end);
% % nmax = length(eng);
% 
% % eng_sub = eng(1:nmax);
% % [a b] = ndgrid([1 2 3 4],[5 6 7 8])
% 
% % eng = eng(1:2:end);
% % nmax = 101;
% 
% 
% % All combinations (index
% [i1, i2, i3]= ndgrid(1:nmax,1:nmax,1:nmax);
% eng1 = eng(i1);
% eng2 = eng(i2);
% eng3 = eng(i3);
% 
% % Calculate total energy
% eng_tot = eng1 + eng2 + eng3;
% 
% % Compress into a single list
% eng_tot = eng_tot(:);
% i1 = i1(:);
% i2 = i2(:);
% i3 = i3(:);
% 
% % Sort energies in ascending order
% [eng_tot,inds] = sort(eng_tot,'ascend');
% 
% % Sort all indeces by ascending order
% i1 = i1(inds);
% i2 = i2(inds);
% i3 = i3(inds);
% 
% f21=figure(21);
% clf
% f21.Color='w';
% 
% plot(eng_tot,'k-','linewidth',1);
% xlabel('index');
% ylabel('energy (kHz)')
% set(gca,'fontsize',14,'fontname','times','box','on');
% 
% % 
% % % density of states
% % f22= figure(22);
% % clf
% % f22.Color='w';
% % 
% % nsmooth = round(length(eng_tot)/100);
% % eng_smooth = smooth(eng_tot,nsmooth);
% % de = diff(eng_smooth);
% % evec = eng_smooth(1:end-1);
% % 
% % plot(evec,1./de,'-');
% 
end

