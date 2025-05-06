%% makeBootstrapSummary.m
% Author : CJ Fujiwara
%
% I made this a script since it is still in testing mode. It also makes
% things more flexible since "final" analyses can be complicated.

%% Structure and Instruction
%
% You must load a variable called composite_data. It should at mininum have
% the fields of composite data

%% Run the bootstrap on the moments 
bs_moments=bootstrap_com(composite_data);

%% Rescale 

%% Run the Bootstrap on the Spectrum
% Only do this if you really mean to, since it will take your computer a
% few hours to run
doRunBootstrap = true;
if doRunBootstrap
    for nn=1:length(composite_data)
        sr=[composite_data(nn).conductivity.cond_real];
        sr_err=[composite_data(nn).conductivity.cond_real_err];        
        si=[composite_data(nn).conductivity.cond_imag];
        si_err=[composite_data(nn).conductivity.cond_imag_err];
        f=[composite_data(nn).conductivity.freq];
        s = sr+1i*si;
        s_err = sr_err+1i*si_err;        
        figure(20);
        clf
        errorbar(f,real(s),real(s_err),'o');
        hold on
        errorbar(f,imag(s),imag(s_err),'o');
        drawnow;
        out(nn)=conductivity_fit_bootstrap(f,s,s_err);
    end

end

%% Create Bootstrap Summary Figures

hF=figure(2);
hF.Name='Bootstrap Summary';
clf
tg = uitabgroup(hF);
t=563;

temp=[];
gamma=[];
trap=[];
rho0=[];
rhoinf=[];
for gg=1:length(out)

    s0=composite_data(gg).Name;
    t0=uicontrol('style','text','string',s0,'horizontalalignment','left',...
        'backgroundcolor','w');
    t0.Position(3:4)=t0.Extent(3:4);
    t0.Position(1:2)=[1 1];

    str=[composite_data(gg).Name(1:10) ' ' num2str(out(gg).SpectralFit.fout(1)/t,'%.1f') 't'];
    tb(gg)=uitab(tg,'Title',str,'backgroundcolor','w');

    subplot(2,3,1,'parent',tb(gg));
    histfit(out(gg).bootstat(:,1)/t);
    xlabel('T/t')
    title('temperature');
    pdT = fitdist(out(gg).bootstat(:,1)/t,'normal');
    temp(gg,1)=pdT.mu; 
    temp(gg,2)=pdT.sigma;
    % xlim([0.5 5])
    text(.01,.99,[num2str(round(pdT.mu,2)) '\pm' num2str(round(pdT.sigma,2))],...
        'units','normalized','verticalalignment','top')
    ylabel('occurences')

    subplot(2,3,2,'parent',tb(gg));
    histfit(out(gg).bootstat(:,2));
    xlabel('\Gamma (1/s)')
    title('Gamma');
    pdG = fitdist(out(gg).bootstat(:,2),'normal');
    gamma(gg,1)=pdG.mu; 
    gamma(gg,2)=pdG.sigma;
    text(.01,.99,[num2str(round(pdG.mu,2)) '\pm' num2str(round(pdG.sigma,2))],...
        'units','normalized','verticalalignment','top')
    ylabel('occurences')


    subplot(2,3,3,'parent',tb(gg));
    histfit(out(gg).bootstat(:,3));
    xlabel('trap frequency (Hz)')
    title('trap frequency');
    pdf = fitdist(out(gg).bootstat(:,3),'normal');
    trap(gg,1)=pdf.mu; 
    trap(gg,2)=pdf.sigma;
    text(.01,.99,[num2str(round(pdf.mu,2)) '\pm' num2str(round(pdf.sigma,2))],...
        'units','normalized','verticalalignment','top')
    ylabel('occurences')


    subplot(2,3,4,'parent',tb(gg));
    histfit(out(gg).bootstat(:,4));
    xlabel('\rho (\sigma_i=0) [1/\sigma_0]')
    title('resonant resitivity');
    pdrho0 = fitdist(out(gg).bootstat(:,4),'normal');
    rho0(gg,1)=pdrho0.mu; 
    rho0(gg,2)=pdrho0.sigma;
    text(.01,.99,[num2str(round(pdrho0.mu,4)) '\pm' num2str(round(pdrho0.sigma,4))],...
        'units','normalized','verticalalignment','top')
    ylabel('occurences')


    subplot(2,3,5,'parent',tb(gg));
    histfit(out(gg).bootstat(:,5));
    xlabel('\rho (\omega\rightarrow \infty) [1/\sigma_0]')
    title('high frequency limit resitivity');
    pdrhoinf = fitdist(out(gg).bootstat(:,5),'normal');
    rhoinf(gg,1)=pdrhoinf.mu; 
    rhoinf(gg,2)=pdrhoinf.sigma;
    text(.01,.99,[num2str(round(pdrhoinf.mu,4)) '\pm' num2str(round(pdrhoinf.sigma,4))],...
        'units','normalized','verticalalignment','top')
    ylabel('occurences')

    
    sr=[composite_data(gg).conductivity.cond_real];
    sr_err=[composite_data(gg).conductivity.cond_real_err];
    
    si=[composite_data(gg).conductivity.cond_imag];
    si_err=[composite_data(gg).conductivity.cond_imag_err];
    f=[composite_data(gg).conductivity.freq];
    
    s = sr+1i*si;
    s_err = sr_err+1i*si_err;


    subplot(2,3,6,'parent',tb(gg));
    co=get(gca,'colororder');

    FREQ_THEORY = linspace(0,200,200);
    yF=conductivity_eval2(FREQ_THEORY, [out(gg).SpectralFit.fout],2.5);
    plot(FREQ_THEORY,real(yF),'-','color',co(1,:));
    hold on
    plot(FREQ_THEORY,imag(yF),'-','color',co(2,:));
    errorbar(f,real(s),real(s_err),'o','color',co(1,:),'markerfacecolor',co(1,:));
    errorbar(f,imag(s),imag(s_err),'o','color',co(2,:),'markerfacecolor',co(2,:));
    xlim([0 150])
    xlabel('drive frequency (Hz)');
    ylabel('conductivity (\sigma_0)')
    title('spectrum');

    fitStr=['lsq fit : $T=' num2str(out(gg).SpectralFit.fout(1)/t,'%.2f') 't,~\Gamma=' num2str(out(gg).SpectralFit.fout(2),'%.1f') '/s,f_0=' num2str(out(gg).SpectralFit.fout(3),'%.1f') '\mathrm{Hz}$'];
    text(.01,.01,fitStr,'units','normalized','interpreter','latex','horizontalalignment','left','verticalalignment','bottom')
end

%% Vesus temperature figures

hF_vs_T=figure(3);
hF_vs_T.Color='w';
clf
hF_vs_T.Name='versus_temp_bootstrap';
co=get(gca,'colororder');

subplot(121);
errorbar(temp(:,1),gamma(:,1),gamma(:,2),gamma(:,2),temp(:,2),temp(:,2),...
    'o','color',co(1,:),'markerfacecolor',co(1,:))
xlabel('spectral T/t [t]')
ylabel('\Gamma [1/s]')
set(gca,'fontsize',14,'fontname','times')
title('current dissipation');


subplot(122);
p0=errorbar(temp(:,1),rho0(:,1),rho0(:,2),rho0(:,2),temp(:,2),temp(:,2),...
    'o','color',.5*co(3,:),'markerfacecolor',co(3,:));
hold on
pinf=errorbar(temp(:,1),rhoinf(:,1),rhoinf(:,2),rhoinf(:,2),temp(:,2),temp(:,2),...
    'o','color',.5*co(4,:),'markerfacecolor',co(4,:));
xlabel('spectral T/t [t]')
ylabel('\rho [1/\sigma_0]')
legend([p0 pinf],{'$\rho(\mathrm{Im}[\sigma]=0)$','$\rho(\omega\rightarrow \infty)$'},'interpreter','latex',...
    'location','northwest')
set(gca,'fontsize',14,'fontname','times')
title('resistivity');
