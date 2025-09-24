%% makeBootstrapSummary.m
% Author : CJ Fujiwara
%
% I made this a script since it is still in testing mode. It also makes
% things more flexible since "final" analyses can be complicated.

%% Structure and Instruction
%
% You must load a variable called composite_data. It should at mininum have
% the fields of composite data

%% Create Bootstrap Summary Figures

hF=figure(667);
hF.Name='Bootstrap Summary';
clf
tg = uitabgroup(hF);
t=563;

tempGRS=[];
gammaGRS=[];
trapGRS=[];
rho0GRS=[];
rhoInfGRS=[];
for gg=1:length(GibbsRescaledOut)
    
    tempGRS_me = GibbsRescaledOut(gg).bootstat(:,1)/t;
    binds=[tempGRS_me<=0.6];
    
    GibbsRescaledOut(gg).bootstat(binds,:)=[];
    
    str=[composite_data(gg).Name(1:10) ' ' num2str(GibbsRescaledOut(gg).SpectralFit.fout(1)/t,'%.1f') 't'];
    tb(gg)=uitab(tg,'Title',str,'backgroundcolor','w');
    s0=composite_data(gg).Name;
    t0=uicontrol('style','text','string',s0,'horizontalalignment','left',...
        'backgroundcolor','w','parent',tb(gg));
    t0.Position(3:4)=t0.Extent(3:4);
    t0.Position(1:2)=[1 1];

    f       = [bs_moments(gg).Frequency_Hz];
    s       = [bs_moments(gg).sigma];
    s_err   = [bs_moments(gg).sigmaErr];

    subplot(2,2,1,'parent',tb(gg));
    co=get(gca,'colororder');
    FREQ_THEORY = linspace(0,200,200);
    yF=conductivity_eval2(FREQ_THEORY, [GibbsRescaledOut(gg).SpectralFit.fout],2.5);
    plot(FREQ_THEORY,real(yF),'-','color',co(1,:));
    hold on
    plot(FREQ_THEORY,imag(yF),'-','color',co(2,:));
    errorbar(f,real(s),real(s_err),'o','color',co(1,:),'markerfacecolor',co(1,:));
    errorbar(f,imag(s),imag(s_err),'o','color',co(2,:),'markerfacecolor',co(2,:));
    xlim([0 150])
    xlabel('drive frequency (Hz)');
    ylabel('conductivity (\sigma_0)')
    title('spectrum');
    fitStr=['lsq fit : $T=' num2str(GibbsRescaledOut(gg).SpectralFit.fout(1)/t,'%.2f') 't,~\gammaGRS=' num2str(GibbsRescaledOut(gg).SpectralFit.fout(2),'%.1f') '/s,f_0=' num2str(GibbsRescaledOut(gg).SpectralFit.fout(3),'%.1f') '\mathrm{Hz}$'];
    text(.01,.01,fitStr,'units','normalized','interpreter','latex','horizontalalignment','left','verticalalignment','bottom')

    subplot(4,4,3,'parent',tb(gg));
    histfit(GibbsRescaledOut(gg).bootstat(:,1)/t);
    xlabel('T/t')
    title('temperature');
    pdT = fitdist(GibbsRescaledOut(gg).bootstat(:,1)/t,'normal');
    tempGRS(gg,1)=pdT.mu; 
    tempGRS(gg,2)=pdT.sigma;
    % xlim([0.5 5])
    text(.01,.99,[num2str(round(pdT.mu,2)) '\pm' num2str(round(pdT.sigma,2))],...
        'units','normalized','verticalalignment','top')
    ylabel('occurences')

    subplot(4,4,4,'parent',tb(gg));
    histfit(GibbsRescaledOut(gg).bootstat(:,2));
    xlabel('\Gamma (1/s)')
    title('Gamma');
    pdG = fitdist(GibbsRescaledOut(gg).bootstat(:,2),'normal');
    gammaGRS(gg,1)=pdG.mu; 
    gammaGRS(gg,2)=pdG.sigma;
    text(.01,.99,[num2str(round(pdG.mu,2)) '\pm' num2str(round(pdG.sigma,2))],...
        'units','normalized','verticalalignment','top')
    ylabel('occurences')


    subplot(4,4,7,'parent',tb(gg));
    histfit(GibbsRescaledOut(gg).bootstat(:,3));
    xlabel('trap frequency (Hz)')
    title('trap frequency');
    pdf = fitdist(GibbsRescaledOut(gg).bootstat(:,3),'normal');
    trapGRS(gg,1)=pdf.mu; 
    trapGRS(gg,2)=pdf.sigma;
    text(.01,.99,[num2str(round(pdf.mu,2)) '\pm' num2str(round(pdf.sigma,2))],...
        'units','normalized','verticalalignment','top')
    ylabel('occurences')


    subplot(4,4,8,'parent',tb(gg));
    histfit(GibbsRescaledOut(gg).bootstat(:,4));
    xlabel('\rho (\sigma_i=0) [1/\sigma_0]')
    title('resonant resitivity');
    pdrho0 = fitdist(GibbsRescaledOut(gg).bootstat(:,4),'normal');
    rho0GRS(gg,1)=pdrho0.mu; 
    rho0GRS(gg,2)=pdrho0.sigma;
    text(.01,.99,[num2str(round(pdrho0.mu,4)) '\pm' num2str(round(pdrho0.sigma,4))],...
        'units','normalized','verticalalignment','top')
    ylabel('occurences')


    subplot(4,4,9,'parent',tb(gg));
    histfit(GibbsRescaledOut(gg).bootstat(:,5));
    xlabel('\rho (\omega\rightarrow \infty) [1/\sigma_0]')
    title('high frequency limit resitivity');
    pdrhoinf = fitdist(GibbsRescaledOut(gg).bootstat(:,5),'normal');
    rhoInfGRS(gg,1)=pdrhoinf.mu; 
    rhoInfGRS(gg,2)=pdrhoinf.sigma;
    text(.01,.99,[num2str(round(pdrhoinf.mu,4)) '\pm' num2str(round(pdrhoinf.sigma,4))],...
        'units','normalized','verticalalignment','top')
    ylabel('occurences')

        subplot(4,4,10,'parent',tb(gg));
    histfit(GibbsRescaledOut(gg).bootstat(:,6)/t);
    xlabel('Tx size/t')
    title('temperature size x');
    pdTx = fitdist(GibbsRescaledOut(gg).bootstat(:,6)/t,'normal');
    TxGRS(gg,1)=pdTx.mu; 
    TxGRS(gg,2)=pdTx.sigma;
    text(.01,.99,[num2str(round(pdTx.mu,2)) '\pm' num2str(round(pdTx.sigma,2))],...
        'units','normalized','verticalalignment','top')
    ylabel('occurences')

        subplot(4,4,11,'parent',tb(gg));
    histfit(GibbsRescaledOut(gg).bootstat(:,7)/t);
    xlabel('Ty size/t')
    title('temperature size y');
    pdTy = fitdist(GibbsRescaledOut(gg).bootstat(:,7)/t,'normal');
    TyGRS(gg,1)=pdTy.mu; 
    TyGRS(gg,2)=pdTy.sigma;
    text(.01,.99,[num2str(round(pdTy.mu,2)) '\pm' num2str(round(pdTy.sigma,2))],...
        'units','normalized','verticalalignment','top')
    ylabel('occurences')

    % Calculate geometric mean harmonic temperature
    TxyGRS(gg,1) = sqrt(TxGRS(gg,1).*TyGRS(gg,1));
    TxyGRS(gg,2) = (TyGRS(gg,1).*TxGRS(gg,2)+TxGRS(gg,1).*TyGRS(gg,2))./(2*sqrt(TxGRS(gg,1).*TyGRS(gg,1)));
    
    % Gibbs
    % TGGRS(gg,1) = TGibbs(gg);
    % TGGRS(gg,2) = TGibbsErr(gg);

end

