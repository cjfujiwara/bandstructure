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
% Analyze the 1st, 2nd, and 3rd moments. 
% Fit the 1st moment to sinuisoidal oscillation
% Fit the 2nd moment to a linear increase (for fitting)

[bs_moments,hF]=bootstrap_com(composite_data);

%% Get Conductivity
omega_xdt   = 2*pi*42;          % [1/s] XDT Trap Frequency
amu         = 1.66054e-27;      % [kg] atomic mass unit
m           = 40*amu;           % [kg] potassium-40 mass
aL          = 532e-9;           % [m] lattice spacing
umperv      = 3.6;              % [um/V] piezo conversion factor
hbar        = 1.05457182e-34;   % [Js] reduce planck's constant

for nn=1:length(bs_moments)
    v2=zeros(length(bs_moments(nn).Params),1);
    for jj=1:length([bs_moments(nn).Params])
        v2(jj) = unique([bs_moments(nn).Params{jj}.conductivity_ODT2_mod_amp]);
    end
    omega           = 2*pi*[bs_moments(nn).Frequency_Hz];omega=omega(:);
    x0              = (v2*umperv*1e-6);
    force_invsec    = (m*omega_xdt^2*x0*aL)/hbar;
    bs_moments(nn).force_invsec ...
                    = force_invsec;
    C_site          =  [bs_moments(nn).C_um]*1e-6/aL;C_site=C_site(:);
    S_site          =  [bs_moments(nn).S_um]*1e-6/aL;S_site=S_site(:);
    CErr_site       =  abs([bs_moments(nn).CErr_um]*1e-6/aL);CErr_site=CErr_site(:);
    SErr_site       =  abs([bs_moments(nn).SErr_um]*1e-6/aL);SErr_site=SErr_site(:);
    bs_moments(nn).sigma ...
                    = -1i*(omega./force_invsec).*(-1i*C_site+S_site);
    c(nn).sigmaErr ...
                    = 1i*(omega./force_invsec).*(-1i*CErr_site+SErr_site);
end

%% Run the Bootstrap on the Spectrum
% Only do this if you really mean to, since it will take your computer a
% few hours to run
doRunBootstrap = true;
if doRunBootstrap
    for nn=1:length(bs_moments)
        % nn=3;
        % Old data to gather
        % sr      = [composite_data(nn).conductivity.cond_real];
        % sr_err  = [composite_data(nn).conductivity.cond_real_err];        
        % si      = [composite_data(nn).conductivity.cond_imag];
        % si_err  = [composite_data(nn).conductivity.cond_imag_err];
        % f       = [composite_data(nn).conductivity.freq];
        % s       = sr+1i*si;
        % s_err   = sr_err+1i*si_err;       

        % f       = [bs_moments(nn).Frequency_Hz];
        % s       = [bs_moments(nn).sigma];
        % s_err   = [bs_moments(nn).sigmaErr];
        % 
        % % Interpreter R2 as an "error" for weighting purposes
        % R2      = [bs_moments(nn).rsquare];
        % R2_err  = sqrt(1./R2);
        % R2_err  = (1+1i)*R2_err.*abs(mean(real(s_err)));        
        % 
        % out(nn)=conductivity_fit_bootstrap(f,s,s_err);

        out(nn)=conductivity_fit_bootstrap(bs_moments(nn));

        % keyboard
    end
end

%% Create Bootstrap Summary Figures

hF=figure(2);
hF.Name='Bootstrap Summary';
clf
tg = uitabgroup(hF);
t = 563;

temp=[];
gamma=[];
trap=[];
rho0=[];
rhoinf=[];
for gg=1:length(composite_data)
    str=[composite_data(gg).Name ' ' num2str(out(gg).SpectralFit.fout(1)/t,'%.1f') 't'];
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




    subplot(4,4,3,'parent',tb(gg));
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

    subplot(4,4,4,'parent',tb(gg));
    histfit(out(gg).bootstat(:,2));
    xlabel('\Gamma (1/s)')
    title('Gamma');
    pdG = fitdist(out(gg).bootstat(:,2),'normal');
    gamma(gg,1)=pdG.mu; 
    gamma(gg,2)=pdG.sigma;
    text(.01,.99,[num2str(round(pdG.mu,2)) '\pm' num2str(round(pdG.sigma,2))],...
        'units','normalized','verticalalignment','top')
    ylabel('occurences')


    subplot(4,4,7,'parent',tb(gg));
    histfit(out(gg).bootstat(:,3));
    xlabel('trap frequency (Hz)')
    title('trap frequency');
    pdf = fitdist(out(gg).bootstat(:,3),'normal');
    trap(gg,1)=pdf.mu; 
    trap(gg,2)=pdf.sigma;
    text(.01,.99,[num2str(round(pdf.mu,2)) '\pm' num2str(round(pdf.sigma,2))],...
        'units','normalized','verticalalignment','top')
    ylabel('occurences')


    subplot(4,4,8,'parent',tb(gg));
    histfit(out(gg).bootstat(:,4));
    xlabel('\rho (\sigma_i=0) [1/\sigma_0]')
    title('resonant resitivity');
    pdrho0 = fitdist(out(gg).bootstat(:,4),'normal');
    rho0(gg,1)=pdrho0.mu; 
    rho0(gg,2)=pdrho0.sigma;
    text(.01,.99,[num2str(round(pdrho0.mu,4)) '\pm' num2str(round(pdrho0.sigma,4))],...
        'units','normalized','verticalalignment','top')
    ylabel('occurences')


    subplot(4,4,9,'parent',tb(gg));
    histfit(out(gg).bootstat(:,5));
    xlabel('\rho (\omega\rightarrow \infty) [1/\sigma_0]')
    title('high frequency limit resitivity');
    pdrhoinf = fitdist(out(gg).bootstat(:,5),'normal');
    rhoinf(gg,1)=pdrhoinf.mu; 
    rhoinf(gg,2)=pdrhoinf.sigma;
    text(.01,.99,[num2str(round(pdrhoinf.mu,4)) '\pm' num2str(round(pdrhoinf.sigma,4))],...
        'units','normalized','verticalalignment','top')
    ylabel('occurences')

        subplot(4,4,10,'parent',tb(gg));
    histfit(out(gg).bootstat(:,6)/t);
    xlabel('Tx size/t')
    title('temperature size x');
    pdTx = fitdist(out(gg).bootstat(:,6)/t,'normal');

    Tx(gg,1)=pdTx.mu; 
    Tx(gg,2)=pdTx.sigma;
    text(.01,.99,[num2str(round(pdTx.mu,2)) '\pm' num2str(round(pdTx.sigma,2))],...
        'units','normalized','verticalalignment','top')
    ylabel('occurences')

        subplot(4,4,11,'parent',tb(gg));
    histfit(out(gg).bootstat(:,7)/t);
    xlabel('Ty size/t')
    title('temperature size y');
    pdTy = fitdist(out(gg).bootstat(:,7)/t,'normal');
    Ty(gg,1)=pdTy.mu; 
    Ty(gg,2)=pdTy.sigma;
    text(.01,.99,[num2str(round(pdTy.mu,2)) '\pm' num2str(round(pdTy.sigma,2))],...
        'units','normalized','verticalalignment','top')
    ylabel('occurences')

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

% errorbar(Tx(:,1),gamma(:,1),gamma(:,2),gamma(:,2),Tx(:,2),Tx(:,2),...
    % 'o','color',co(1,:),'markerfacecolor',co(1,:))


xlabel('spectral T/t [t]')
ylabel('\Gamma [1/s]')
set(gca,'fontsize',14,'fontname','times')
title('current dissipation');
xlim([0 4])


subplot(122);
p0=errorbar(temp(:,1),rho0(:,1),rho0(:,2),rho0(:,2),temp(:,2),temp(:,2),...
    'o','color',.5*co(3,:),'markerfacecolor',co(3,:));
% errorbar(Tx(:,1),rho0(:,1),rho0(:,2),rho0(:,2),Tx(:,2),Tx(:,2),...
    % 'o','color',co(1,:),'markerfacecolor',co(1,:))
hold on
% pinf=errorbar(temp(:,1),rhoinf(:,1),rhoinf(:,2),rhoinf(:,2),temp(:,2),temp(:,2),...
    % 'o','color',.5*co(4,:),'markerfacecolor',co(4,:));
xlabel('spectral T/t [t]')
ylabel('\rho [1/\sigma_0]')
% legend([p0 pinf],{'$\rho(\mathrm{Im}[\sigma]=0)$','$\rho(\omega\rightarrow \infty)$'},'interpreter','latex',...
%     'location','northwest')

set(gca,'fontsize',14,'fontname','times')
title('resistivity');
xlim([0 4])
