function output_best = conductivity_fit(freq,sigma)
% freq      : frequency data
% sigma     : complex conductivity data

% Make sure they are a colummn vector
freq = freq(:);     
sigma = sigma(:);

% Separate real and imaginary parts
y = [real(sigma); imag(sigma)];

%% Lattice Properties
lattice                     = constants;
lattice.depth               = [2.5]; 
lattice.WannierBands        = [1];

lattice.numStates           = 101;       % must be odd
lattice.numK                = 301;      % must be odd    

wannier_opts                = struct;
wannier_opts.Bands          = [1];

%% Calculate Lattice Properties

lattice = calculateBandStructure(lattice);      % calculate band structure
lattice = calculateTunneling(lattice);          % calculate tunneling elements
lattice = wannier(lattice,wannier_opts);        % Calculate wannier function
lattice = calculateWannierMoments(lattice);     % Dipole matrix elements in wannier basis


% First tunneling in Hz
tunnelings = lattice.Tunneling(1,:)*lattice.fr;% 1st tunneling
t= tunnelings(1);
%% LHO
% Number of eigenstates to include in fit
N = 101;  
% Numerical Settings
Nsites          = 601;
TunnelOrder     = 21;
HarmonicBands   = 1;    

%% Display Stuff

figNum1 = 1990;
figNum2 = 1991;

%% Drude Fit
% Fit the drude peak in a combined sort of way

% Complex drude wrapper function
drude_wrapper = @(P,f) [real(drude(P(1),P(2),P(3),f)); 
    imag(drude(P(1),P(2),P(3),f))];

% Complex drude function
function y = drude(amp,f0,G,f)
    w   = 2*pi*f;
    w0  = 2*pi*f0;
    y   = 1i*amp*G*w./((w.^2-w0.^2)+1i*w*G);
end

% Contstruct guess
[amp,ind]       = max(real(sigma)); % Peak response
f0              = freq(ind);        % Frequency of peak response
P_drude_complex = [amp f0 50];      % Guess object

[fout,resnorm,residual,exitflag,output0,lambda,jacobian]= ...
    lsqcurvefit(drude_wrapper,P_drude_complex,freq,y);
conf = nlparci(fout,residual,'jacobian',jacobian);

% Rsquared of fit
SS_res = resnorm;
SS_tot = sum((real(sigma)-mean(real(sigma))).^2+(imag(sigma)-mean(imag(sigma))).^2);
R2 = 1 - SS_res/SS_tot;

amp         = fout(1);                      % Peak              [sigma_0]
amperr      = (conf(1,2)-conf(1,1))/2;      % Peak Error        [sigma_0]
f0          = fout(2);                      % Peak Freq         [Hz]
f0err       = (conf(2,2)-conf(2,1))/2;      % Peak Freq Error   [Hz]
G           = fout(3);                      % Gamma             [1/s]
Gerr        = (conf(3,2)-conf(3,1))/2;      % Gamma Err         [1/s]

ft = linspace(0,2*max(freq),2000);
drude_complex               = struct;
drude_complex.s0            = amp;
drude_complex.s0err         = amperr;
drude_complex.f0            = f0;
drude_complex.f0err         = f0err;
drude_complex.G             = G;
drude_complex.Gerr          = Gerr;
drude_complex.Rsquared      = R2;
drude_complex.FreqFit       = ft;
drude_complex.SigmaFit      = drude(amp,f0,G,ft);
% drude_sum                   = trapz(ft,real(drude(amp,f0,G,ft)));


%% Calculate LHO States


omega           = 2*pi*[(f0+2):1:(f0+12)];

HO_opts = struct;
HO_opts.NumSites = Nsites;
HO_opts.MaxTunnelingOrder = TunnelOrder;
HO_opts.HarmonicBands = HarmonicBands;
for nn=1:length(omega)
    HO_opts.omega = omega(nn);
    HO_opts.Omega = 0.5*lattice.m*omega(nn)^2*(lattice.lambda/2)^2/lattice.h;
    [lattice,LHO(nn)] = calculateLHOSpectrum_sband(lattice,HO_opts);
end    

%% Conductivity Helper Functions

% Total conductivity
    function y = sigma_func(T,G,f)
        y = arrayfun(@( f) sigma_helper(T,G,f),f);
    end

    function y = sigma_helper(T,G,f)
        Z = sum(exp(-eng/T),'all');
        A = -1i*f*((exp(-EE1/T)-exp(-EE2/T))/Z).*d2./((f-dEE)+1i*G/2/(2*pi));
        y = sum(A,'all');
    end

% Real conductivity
    function y = sigma_real(T,G,f)        
        y = real(sigma_func(T,G,f));
    end

% Imaginary conductivity
    function y = sigma_imag(T,G,f)
        y = imag(sigma_func(T,G,f));
    end

%%
output = struct;
hF1 = figure(figNum1);
hF1.Color='w';
hF1.Position=[500 50 1200 300];
clf
ax1 = subplot(1,4,1,'parent',hF1);
ax2 = subplot(1,4,2,'parent',hF1);
ax3 = subplot(1,4,3,'parent',hF1);
ax4 = subplot(1,4,4,'parent',hF1);
cc=jet(length(LHO));


pDR=plot(freq,real(sigma),'ko','markerfacecolor','k','parent',ax1);
hold(ax1,'on');
xlabel(ax1,'Frequency (Hz)')
ylabel(ax1,'Re[\sigma/\sigma_0]')

pDI=plot(freq,imag(sigma),'ko','markerfacecolor','k','parent',ax2);
xlabel(ax2,'Frequency (Hz)')
ylabel(ax2,'Im[\sigma/\sigma_0]')
hold(ax2,'on');

for jj=1:length(LHO)
    fprintf(['Fitting ' num2str(jj) ' of ' num2str(length(LHO))]);
    src = LHO(jj);
    
    % Load Eigenvalues and Dipole Operator
    d2          = abs(src.DipoleOperator).^2;
    eng         = src.EigenValues;
    eng         = eng-eng(1);    
    eng         = eng(1:N);
    d2          = d2(1:N,1:N);
    
    % Meshgrid data
    [EE1,EE2]   = meshgrid(eng,eng);
    dEE         = EE1-EE2;
    
    P = [850 drude_complex.G];
    % P = [1100 51];
    
    %% Constrained Fit
    tdpt_wrapper = @(P,f) [sigma_real(P(1),P(2),f); sigma_imag(P(1),P(2),f)];
    options = optimset('Display','off');    

    [fout,resnorm,residual,exitflag,output0,lambda,jacobian]=lsqcurvefit(tdpt_wrapper,P,freq,y,[],[],options);
    conf = nlparci(fout,residual,'jacobian',jacobian);
    
    SS_res = resnorm;
    SS_tot = sum((real(sigma)-mean(real(sigma))).^2+(imag(sigma)-mean(imag(sigma))).^2);
    R2 = 1 - SS_res/SS_tot;
    
    T = fout(1);
    Terr = (conf(1,2)-conf(1,1))/2;
    
    G = fout(2);
    Gerr = (conf(2,2)-conf(2,1))/2;
    output(jj).TrapFrequency_Hz     = src.omega/(2*pi);
    output(jj).Temperature_Hz       = T;
    output(jj).TemperatureErr_Hz    = Terr;
    output(jj).Gamma_invSec         = G;
    output(jj).GammaErr_invSec      = Gerr;
    output(jj).Rsquared             = R2;
    

    yt = sigma_func(T,G,ft);      
    plot(ft,real(yt),'-','color',cc(jj,:),'parent',ax1);
    hold(ax1,'on')
    plot(ft,imag(yt),'-','color',cc(jj,:),'parent',ax2);
    hold(ax2,'on')
    disp(' done');
    drawnow;
end


plot(omega/(2*pi),[output.Rsquared],'o','parent',ax3)
xlabel(ax3,'trap freq (Hz)')
ylabel(ax3,'R squared')

uistack(pDI,'top')
uistack(pDR,'top')

axes(ax4);
yyaxis left
errorbar(omega/(2*pi),[output.Gamma_invSec],[output.GammaErr_invSec],'o','parent',ax4)
xlabel(ax4,'trap freq (Hz)')
ylabel(ax4,'\Gamma (1/s)')
yyaxis right
errorbar(omega/(2*pi),[output.Temperature_Hz]/t,[output.TemperatureErr_Hz]/563,'o','parent',ax4)
ylabel('temp (t)')
ylim([0 3])

%% Summary Figure
[val,ind]           = max([output.Rsquared]);
output_best         = output(ind);
T                   = output_best.Temperature_Hz;
Terr                = output_best.TemperatureErr_Hz;
G                   = output_best.Gamma_invSec;
Gerr                = output_best.GammaErr_invSec;
f0_best             = output_best.TrapFrequency_Hz;

sR = @(omega) real(sigma_func(T,G,omega/(2*pi)));

Fsum = 2/pi*integral(sR,0,2000);



%% Best Summary
hF2 = figure(figNum2);
hF2.Color='w';
clf(hF2);
axes('parent',hF2)
hF2.Position=[5 50 500 300];
co=get(gca,'colororder');
% axa = subplot(2,1,1,'parent',hF2);

pDrude_R=plot(ft,real(drude(amp,f0,G,ft)),'--','color',co(1,:));
hold on
pDrude_I=plot(ft,imag(drude(amp,f0,G,ft)),'--','color',co(2,:));
ylabel('$\sigma/\sigma_0$','interpreter','latex')
xlabel('drive frequency (Hz)')
set(gca,'box','on','fontname','times')

pTDPT_R=plot(ft,real(sigma_func(T,G,ft)),'-','color',co(1,:),'linewidth',1);
hold on
pTDPT_I=plot(ft,imag(sigma_func(T,G,ft)),'-','color',co(2,:),'linewidth',1);

pData_R=plot(freq,real(sigma),'o','markerfacecolor',co(1,:),'color',co(1,:)*.5);
pData_I=plot(freq,imag(sigma),'s','markerfacecolor',co(2,:),'color',co(2,:)*.5);

str1=['Re[drude]'];
str2=['Im[drude]'];

str3=['Re[tdpt]'];
str4=['Im[tdpt]'];

legStr={str1, str2, str3, str4};
legend([pDrude_R pDrude_I pTDPT_R pTDPT_I],legStr)

str_drude = ['drude $(' ...
    'f_0:' num2str(round(drude_complex.f0,1)) '~\mathrm{Hz}' ...
    ',A:' num2str(round(drude_complex.s0,1))  ...    
    ',\Gamma:' num2str(round(drude_complex.G,1)) '\mathrm{s}^{-1})$'];

str_tdpt = ['tdpt $(f_0:' num2str(round(f0_best,1)) '~\mathrm{Hz}' ...
    ',T:' num2str(round(T/t,2)) 't'...
    ',\Gamma:' num2str(round(G,1)) '\mathrm{s}^{-1})$'];

str_tdpt_2 = ['tdpt err $\pm' num2str(round(Terr/t,1)) 't'...
    ',\pm' num2str(round(Gerr,0)) '\mathrm{s}^{-1})$'];

str = [str_drude newline str_tdpt newline str_tdpt_2];

text(.99,.01,str,'units','normalized','horizontalalignment','right',...
    'verticalalignment','bottom','interpreter','latex')
xlim([0 150]);

end

